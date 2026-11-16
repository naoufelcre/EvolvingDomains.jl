"""
    CutSurface{T,N}

Piecewise-linear discretization of the interface Γ = {x ∈ Ω : φ(x) = 0}, using segments in 2D
and triangles in 3D.

- `point_ids`, `point_pointers`: facet-to-vertex connectivity, stored as an incidence matrix
  in the CSR format.
- `physical_points`: interface vertex coordinates.
- `reference_points`: vertex coordinates in the reference cell.
- `facet_to_bgcell`: mapping from facets to background grid cells.
- `normals`: outward normal vector of each facet, from the `φ < 0` phase.
"""
struct CutSurface{T,N}
    facet_to_bgcell::Vector{Int32}
    point_ids::Vector{Int32}
    point_pointers::Vector{Int32}
    physical_points::Vector{NTuple{N,T}}
    reference_points::Vector{NTuple{N,T}}
    normals::Vector{NTuple{N,T}}
end

"""
    Cut(grid, level_id; kwargs...)

Initialize a `Cut` from a grid and the arrays of an external cut. The cut is expected to come
from an external cutter such as [CutCells](https://github.com/sclaus2/CutCells): for every
background cell it returns the cell's classification, its volume subcells with their vertices,
and its interface facets. [`cut(field, level_id)`](@ref) is the normal entry point and fills
all of this from a level-set `Field`.
"""
struct Cut{T,N}
    mesh_id::UInt64
    level_id::Symbol
    bg_states::Vector{Int8}
    subcell_states::Vector{Int8}
    cell_to_bgcell::Vector{Int32}
    cell_subcell_ptrs::Vector{Int32}
    point_ids::Vector{Int32}
    point_pointers::Vector{Int32}
    physical_points::Vector{NTuple{N,T}}
    reference_points::Vector{NTuple{N,T}}
    boundary::CutSurface{T,N}
    # Declared abstract on purpose: `NTuple{2N,Int8}` is not expressible while N is a
    # TypeVar, and making it concrete would add a third type parameter to every signature.
    cell_facet_states::Union{Nothing,AbstractVector{<:Tuple}}

    Cut{T,N}(::Val{:validated}, mesh_id, level_id, bg_states, subcell_states,
           cell_to_bgcell, cell_subcell_ptrs, point_ids, point_pointers,
           physical_points, reference_points, boundary, cell_facet_states) where {T,N} =
        new{T,N}(mesh_id, level_id, bg_states, subcell_states, cell_to_bgcell,
               cell_subcell_ptrs, point_ids, point_pointers, physical_points,
               reference_points, boundary, cell_facet_states)
end

"""
    _cutSubcellPtrs(parents::Vector{Int32}, ncells::Int)

Build one-based CSR row pointers for subcells grouped by ascending parent-cell index.
"""
function _cutSubcellPtrs(parents::Vector{Int32}, ncells::Int)
    ptrs = zeros(Int32, ncells + 1)
    previous = 1
    for cell in parents
        current = Int(cell)
        current >= previous || throw(ArgumentError(
            "subcells must be grouped by increasing background-cell order"))
        previous = current
        ptrs[current + 1] += 1
    end
    for cell in 1:ncells                    # exclusive prefix sum, shifted to one-based rows
        ptrs[cell + 1] += ptrs[cell]
    end
    ptrs .+= 1
    return ptrs
end

"""    subcells(cut::Cut, cell::Integer) -> UnitRange{Int}

Ids of the subcells owned by background `cell`, empty when it owns none. Cells are identified
by a linear index; see the conventions."""
subcells(cut::Cut, cell::Integer) =
    Int(cut.cell_subcell_ptrs[cell]):(Int(cut.cell_subcell_ptrs[cell + 1]) - 1)

"""    cut_cell_ids(cut::Cut) -> Vector{Int32}

Ascending one-based linear (first-axis-fastest) indices of the background cells classified
`CUT_CUT`."""
cut_cell_ids(cut::Cut) = Int32.(findall(==(CUT_CUT), cut.bg_states))

"""`Int8(-1)`: background cell entirely in `Ω⁻` (all corners `φ < 0`)."""
const CUT_IN = Int8(-1)

"""`Int8(1)`: background cell entirely in `Ω⁺` (all corners `φ > 0`)."""
const CUT_OUT = Int8(1)

"""`Int8(0)`: background cell straddling `Γ`; see `Cut.bg_states`."""
const CUT_CUT = Int8(0)

# Internal aliases kept while the remaining call sites migrate to the public names.
const _CUT_IN = CUT_IN
const _CUT_OUT = CUT_OUT
const _CUT_CUT = CUT_CUT

function _cutState(state, name)
    state isa Integer && !(state isa Bool) &&
        (state == _CUT_IN || state == _CUT_OUT || state == _CUT_CUT) ||
        throw(ArgumentError("$name contains a state other than IN (-1), OUT (1), or CUT (0)"))
    return Int8(state)
end

_cutStates(values, name) = Int8[_cutState(state, name) for state in values]

function _cutConnectivity(ids, pointers, npoints, ncells, arity, name)
    data = Int32.(collect(ids))
    ptrs = Int32.(collect(pointers))
    !isempty(ptrs) && first(ptrs) == 1 ||
        throw(ArgumentError("$name pointers must start at one"))
    length(ptrs) == ncells + 1 ||
        throw(ArgumentError("$name pointers must have one entry per item plus a terminal entry"))
    last(ptrs) == length(data) + 1 ||
        throw(ArgumentError("$name terminal pointer must equal length(point_ids)+1"))
    all(i -> 1 <= i <= npoints, data) ||
        throw(ArgumentError("$name point indices are out of range"))
    all(i -> ptrs[i] < ptrs[i + 1] && ptrs[i + 1] - ptrs[i] == arity, 1:ncells) ||
        throw(ArgumentError("each $name must have exactly $arity vertices"))
    return data, ptrs
end

function _validateCutPointArity(points, name, ndim)
    for (index, point) in enumerate(points)
        components = try
            length(point)
        catch
            throw(ArgumentError("$name coordinate $index must have exactly $ndim components"))
        end
        components == ndim ||
            throw(ArgumentError("$name coordinate $index must have exactly $ndim components"))
    end
    return nothing
end

function _cutPoints(points, ::Type{T}, name, ndim) where {T}
    _validateCutPointArity(points, name, ndim)
    result = NTuple{ndim,T}[ntuple(i -> T(p[i]), ndim) for p in points]
    all(p -> all(isfinite, p), result) ||
        throw(ArgumentError("$name must contain only finite coordinates"))
    result
end

function Cut(grid::Grid, level_id::Symbol;
             bg_states, subcell_states=Int8[], cell_to_bgcell=Int32[],
             cell_facet_states=nothing,
             point_ids=Int32[], point_pointers=Int32[1],
             physical_points=NTuple{2,Float64}[], reference_points=NTuple{2,Float64}[],
             boundary_facet_to_bgcell=Int32[], boundary_point_ids=Int32[],
             boundary_point_pointers=Int32[1], boundary_physical_points=NTuple{2,Float64}[],
             boundary_reference_points=NTuple{2,Float64}[], boundary_normals=NTuple{2,Float64}[])
    mesh_id = grid.mesh_id
    isempty(String(level_id)) && throw(ArgumentError("level_id must be nonempty"))
    axes = grid.axes
    ndim = length(axes)
    ndim in (2, 3) || throw(ArgumentError("Cut supports only two- or three-dimensional grids"))
    sizes = ntuple(d -> length(axes[d]) - 1, ndim)
    ncells = prod(sizes)
    backgrounds = _cutStates(bg_states, "bg_states")
    length(backgrounds) == ncells ||
        throw(ArgumentError("bg_states length must equal the number of background cells"))
    facet_states = if isnothing(cell_facet_states)
        nothing
    else
        length(cell_facet_states) == ncells ||
            throw(ArgumentError("cell_facet_states length must equal the number of background cells"))
        states = Vector{NTuple{2ndim,Int8}}(undef, ncells)
        for (cell, sides) in enumerate(cell_facet_states)
            length(sides) == 2ndim || throw(ArgumentError(
                "cell_facet_states[$cell] must contain 2 states per dimension in (-x,+x,-y,+y,-z,+z) order"))
            states[cell] = ntuple(i -> _cutState(sides[i], "cell_facet_states[$cell]"), 2ndim)
        end
        for cell in 1:ncells, d in 1:ndim
            indices = ntuple(a -> ((cell - 1) ÷ prod(sizes[1:a-1]) % sizes[a]) + 1, ndim)
            if indices[d] < sizes[d]
                neighbor = cell + prod(sizes[1:d-1])
                states[cell][2d] == states[neighbor][2d-1] ||
                    throw(ArgumentError("neighboring cells disagree on their shared facet state"))
            end
        end
        states
    end
    sub_states = _cutStates(subcell_states, "subcell_states")
    parents = Int32.(collect(cell_to_bgcell))
    length(parents) == length(sub_states) ||
        throw(ArgumentError("cell_to_bgcell and subcell_states lengths must match"))
    all(c -> 1 <= c <= ncells, parents) ||
        throw(ArgumentError("subcell parent-cell index is out of range"))

    _validateCutPointArity(physical_points, "physical_points", ndim)
    _validateCutPointArity(reference_points, "reference_points", ndim)
    _validateCutPointArity(boundary_physical_points, "boundary_physical_points", ndim)
    _validateCutPointArity(boundary_reference_points, "boundary_reference_points", ndim)
    _validateCutPointArity(boundary_normals, "boundary_normals", ndim)
    coordinateType(points) = isempty(points) ? eltype(first(axes)) :
        foldl(promote_type, (typeof(first(points)[d]) for d in 1:ndim))
    T = foldl(promote_type, (eltype(axis) for axis in axes))
    T = promote_type(T,
                     coordinateType(physical_points), coordinateType(reference_points),
                     coordinateType(boundary_physical_points),
                     coordinateType(boundary_reference_points), coordinateType(boundary_normals))
    isconcretetype(T) || (T = Float64)
    physical = _cutPoints(physical_points, T, "physical_points", ndim)
    reference = _cutPoints(reference_points, T, "reference_points", ndim)
    length(physical) == length(reference) ||
        throw(ArgumentError("physical_points and reference_points lengths must match"))
    ids, ptrs = _cutConnectivity(point_ids, point_pointers, length(physical), length(sub_states), ndim + 1, "subcell")

    facet_parents = Int32.(collect(boundary_facet_to_bgcell))
    normals = _cutPoints(boundary_normals, T, "boundary_normals", ndim)
    bphysical = _cutPoints(boundary_physical_points, T, "boundary_physical_points", ndim)
    breference = _cutPoints(boundary_reference_points, T, "boundary_reference_points", ndim)
    length(normals) == length(facet_parents) ||
        throw(ArgumentError("boundary normals and facet parent-cell lengths must match"))
    length(bphysical) == length(breference) ||
        throw(ArgumentError("boundary physical/reference point lengths must match"))
    all(c -> 1 <= c <= ncells, facet_parents) ||
        throw(ArgumentError("boundary parent-cell index is out of range"))
    issorted(facet_parents) || throw(ArgumentError("boundary facets must be grouped by increasing background-cell order"))
    bids, bptrs = _cutConnectivity(boundary_point_ids, boundary_point_pointers,
                                    length(bphysical), length(facet_parents), ndim, "boundary facet")
    # Coherence: the arrays must describe one geometric situation, not merely have valid
    # lengths. A split cell owns volume pieces; no piece belongs to a cell declared outside.
    # The one cross-array rule the producer guarantees: a split cell owns volume pieces.
    # Every background cell may own them (an entirely outside cell emits its own outside
    # volume), so the converse is not an invariant. A boundary-only Cut carries no volume
    # subcells at all, which is legal and used by reinitialization and transport masks.
    if !isempty(sub_states)
        owners = Set(parents)
        for cell in 1:ncells
            if backgrounds[cell] == CUT_CUT && !(cell in owners)
                throw(ArgumentError("background cell $cell is CUT but owns no subcell"))
            end
        end
    end
    index = _cutSubcellPtrs(parents, ncells)
    return Cut{T,ndim}(Val(:validated), mesh_id, level_id, backgrounds, sub_states,
                  parents, index, ids, ptrs, physical, reference,
                  CutSurface{T,ndim}(facet_parents, bids, bptrs, bphysical, breference, normals),
                  facet_states)
end

"""    _cutFromNative(grid, level_id; kwargs...)

Return, as a `Cut` snapshot, the cut produced by the bound CutCells C++ library.
"""
function _cutFromNative(grid::Grid{D,T}, level_id::Symbol;
                        bg_states, subcell_states, cell_to_bgcell, cell_facet_states,
                        point_ids, point_pointers, physical_points, reference_points,
                        boundary_facet_to_bgcell, boundary_point_ids,
                        boundary_point_pointers, boundary_physical_points,
                        boundary_reference_points, boundary_normals) where {D,T}
    boundary = CutSurface{T,D}(boundary_facet_to_bgcell, boundary_point_ids,
        boundary_point_pointers, boundary_physical_points, boundary_reference_points,
        boundary_normals)
    index = _cutSubcellPtrs(cell_to_bgcell, length(bg_states))
    return Cut{T,D}(Val(:validated), grid.mesh_id, level_id, bg_states, subcell_states,
        cell_to_bgcell, index, point_ids, point_pointers, physical_points,
        reference_points, boundary, cell_facet_states)
end
