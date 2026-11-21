using ..EvolvingDomains: Cut, Grid, CUT_IN, CUT_CUT

function _cutCellVertices(grid::Grid{D}, cell) where {D}
    index = CartesianIndices(ntuple(d -> length(grid.axes[d]) - 1, D))[cell]
    # Vertex v is the corner lifted by one cell along every axis whose bit is set in v - 1.
    ntuple(v -> ntuple(d -> grid.axes[d][index[d] + (((v - 1) >> (d - 1)) & 1)], D), 2^D)
end

function _cutCellDistance(grid::Grid{D}, cell, root) where {D}
    cell_vertices = _cutCellVertices(grid, cell)
    root_vertices = _cutCellVertices(grid, root)
    maximum(Float64(hypot(ntuple(d -> p[d] - q[d], D)...))
            for p in cell_vertices for q in root_vertices)
end

function _cutCellNeighbor(cell, side, dims::NTuple{D,Int}) where {D}
    1 <= side <= 2D || error("invalid Cartesian side index $side")
    axis = (side - 1) ÷ 2 + 1
    stride = prod(dims[1:axis - 1])
    index = (cell - 1) ÷ stride % dims[axis] + 1
    isodd(side) && return index > 1 ? cell - stride : 0
    index < dims[axis] ? cell + stride : 0
end

function _touchCutAggregatedCells!(touched, roots)
    @inbounds for cell in eachindex(roots)
        roots[cell] > 0 && (touched[cell] = true)
    end
    return touched
end

# Gridap local-face conventions: QUAD walks bottom, top, left, right; HEX uses axis order.
_cutFacetOrder(::Val{2}) = ((3, 4, 1, 2), "(left,right,bottom,top)")
_cutFacetOrder(::Val{3}) = ((1, 2, 3, 4, 5, 6), "(-x,+x,-y,+y,-z,+z)")

# Arity-error wording, kept verbatim per dimension.
_cutSubcellKind(::Val{2}) = "triangular"
_cutSubcellKind(::Val{3}) = "tetrahedral"

# Determinant of the N x N edge matrix `m` (row r, column c), expanding along the
# first column: the N-simplex measure is `abs(_cutSimplexDet(m)) / factorial(N)`.
_cutSimplexMinor(m::NTuple{N,NTuple{N}}, row) where {N} =
    ntuple(r -> ntuple(c -> m[r + (r >= row)][c + 1], N - 1), N - 1)

_cutSimplexDet(m::NTuple{1,NTuple{1}}) = m[1][1]

function _cutSimplexDet(m::NTuple{N,NTuple{N}}) where {N}
    det = m[1][1] * _cutSimplexDet(_cutSimplexMinor(m, 1))
    for row in 2:N
        term = m[row][1] * _cutSimplexDet(_cutSimplexMinor(m, row))
        det = iseven(row) ? det - term : det + term
    end
    det
end

function _cutPhysicalCellMeasures(cut::Cut{T,N}, n_cells) where {T,N}
    measures = zeros(Float64, n_cells)
    for subcell in eachindex(cut.subcell_states)
        cut.subcell_states[subcell] != CUT_IN && continue
        first_point = Int(cut.point_pointers[subcell])
        next_point = Int(cut.point_pointers[subcell + 1])
        next_point == first_point + N + 1 || error(
            "cut-cell aggregation requires $(_cutSubcellKind(Val(N))) cut subcells")
        points = ntuple(q -> cut.physical_points[cut.point_ids[first_point + q - 1]], N + 1)
        edges = ntuple(r -> ntuple(c -> points[c + 1][r] - points[1][r], N), N)
        measures[Int(cut.cell_to_bgcell[subcell])] += abs(_cutSimplexDet(edges)) / factorial(N)
    end
    measures
end

"""
    aggregateCutCells(cut, grid)

Provided a `Cut` and its associated grid, use the classification of the cells and the state of
their facets to map every cut cell to a neighbouring in-phase cell used as root. The aggregation
algorithm is a simple greedy search; the result is meant for AgFEM: Badia, Martín & Verdugo,
"The aggregated unfitted finite element method for elliptic problems", Computer Methods in
Applied Mechanics and Engineering 336 (2018) 480–533 (unverified). For that purpose, any aggregation of cut
cells onto neighbouring in-phase roots is admissible. A zero in the returned mapping means the
cell belongs to no aggregate, i.e. it lies outside the phase.
"""
function aggregateCutCells(cut::Cut{T,N}, grid::Grid{D}) where {T,N,D}
    _sameid(cut, grid) || throw(ArgumentError("Cut and Grid mesh identities differ"))
    D in (2, 3) || throw(ArgumentError("Cut aggregation supports only 2D or 3D grids"))
    sides, facet_order = _cutFacetOrder(Val(D))
    dims = ntuple(d -> length(grid.axes[d]) - 1, D)
    n_cells = prod(dims)
    length(cut.bg_states) == n_cells || throw(ArgumentError(
        "Cut background-cell state count does not match Grid dimensions"))
    facet_states = cut.cell_facet_states
    isnothing(facet_states) && throw(ArgumentError(
        "Cut aggregation requires cell_facet_states in $facet_order order"))
    length(facet_states) == n_cells || throw(ArgumentError(
        "Cut cell_facet_states count does not match Grid dimensions"))

    roots = zeros(Int32, n_cells)
    touched = falses(n_cells)
    for cell in 1:n_cells
        if cut.bg_states[cell] == CUT_IN
            roots[cell] = cell
            touched[cell] = true
        end
    end

    fallback_roots = Int[]
    measures = nothing
    for _ in 1:n_cells
        unaggregated = 0
        made_progress = false
        for cell in 1:n_cells
            if !touched[cell] && cut.bg_states[cell] == CUT_CUT
                best_neighbor = 0
                dmin = Inf
                for side in sides
                    facet_states[cell][side] in (CUT_CUT, CUT_IN) || continue
                    neighbor = _cutCellNeighbor(cell, side, dims)
                    if neighbor > 0 && touched[neighbor]
                        distance = _cutCellDistance(grid, cell, roots[neighbor])
                        if (1 + 1e-9) * distance < dmin
                            dmin = distance
                            best_neighbor = neighbor
                        end
                    end
                end
                if best_neighbor > 0
                    roots[cell] = roots[best_neighbor]
                    made_progress = true
                else
                    unaggregated += 1
                end
            end
        end
        unaggregated == 0 && break
        _touchCutAggregatedCells!(touched, roots)

        if !made_progress
            isnothing(measures) && (measures = _cutPhysicalCellMeasures(cut, n_cells))
            root = 0
            best_measure = 0.0
            for cell in 1:n_cells
                if !touched[cell] && cut.bg_states[cell] == CUT_CUT && measures[cell] > best_measure
                    root = cell
                    best_measure = measures[cell]
                end
            end
            root > 0 || break
            roots[root] = root
            touched[root] = true
            push!(fallback_roots, root)
        end
    end

    unresolved = findall(i -> roots[i] == 0 && cut.bg_states[i] == CUT_CUT, 1:n_cells)
    isempty(unresolved) || error(
        "cut-cell aggregation failed for $(length(unresolved)) cells: " *
        join(unresolved[1:min(10, length(unresolved))], ", "))
    isempty(fallback_roots) || @warn(
        "AgFEM component has no interior cell; using its largest cut cell as root",
        roots=fallback_roots, maxlog=3)
    # Single-signed cut: nothing is cut, so there is no aggregate to name. Say so, rather
    # than handing back a full-length vector in which every entry is the "unassigned" 0.
    all(iszero, roots) && @warn(
        "aggregateCutCells: this cut has no interface, so no cell belongs to an " *
        "aggregate; every root is 0 (unassigned)", maxlog=1)
    roots
end
