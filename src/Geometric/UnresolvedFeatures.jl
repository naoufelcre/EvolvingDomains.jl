# Luo, Hu & Adams (2016) scale separation for level-set interfaces, nodal adaptation.
#
# Published method (J. Comput. Phys. 308 (2016) 411–420, §3):
#   S   = cells cut by phi = 0                                            Eq. (7)
#   P,N = cells cut by phi = +eps and phi = -eps                          Eq. (9)
#   S1  = cut cells whose 3^N neighbourhood misses P or misses N          Eq. (10)
#   S1P = cells of S1 whose neighbourhood misses P                        Eq. (12)
#   S1N = cells of S1 whose neighbourhood misses N                        Eq. (13)
#   a cell of S1P is re-distanced from the cells of P inside a 9^N search region:
#       q_p   = x_p + (phi_p - eps) n_p                                   Eq. (14)
#       phi_s = eps - min_p |x_s - q_p|                                   Eq. (15)
#   a cell of S1N is re-distanced from N with -eps and the opposite sign.
#
# Adaptations:
#   * nodal storage: a cell's corners are grid nodes and its value is their average; the
#     correction of Eq. (15) is evaluated at the cell's nodes, not at a cell centre.
#   * when the required auxiliary side is absent from the whole search region, Eq. (15) has
#     no source. An isolated island — connected same-sign region that never reaches eps —
#     is then replaced by the surrounding phase. A structure that reaches the auxiliary
#     level set elsewhere is never handled this way.

module UnresolvedFeatures

using StaticArrays: SVector
using LinearAlgebra: norm

using ..Geometric: _gridMetrics, _GridMetrics, _checklevel
using ....EvolvingDomains: Field, _levelvalues

const _EPS_FACTOR = 0.75        # eps = 0.75 h, Luo et al. (2016), Eq. (9)
const _NEIGHBOURHOOD = 1        # 3x3 in 2D, 3x3x3 in 3D
const _SEARCH = 4               # 9x9 in 2D, 9x9x9 in 3D

_cellDims(dims::NTuple{N,Int}) where {N} = ntuple(d -> dims[d] - 1, N)
_nodeDims(cell_dims::NTuple{N,Int}) where {N} = ntuple(d -> cell_dims[d] + 1, N)

@inline _corners(::Val{N}) where {N} = CartesianIndices(ntuple(_ -> 1:2, N))

@inline _corner(c::CartesianIndex{N}, k::CartesianIndex{N}) where {N} =
    CartesianIndex(ntuple(d -> c[d] + k[d] - 1, N))

# The paper's cell test: corners on both sides of `level` (Eq. 7 and 9).
@inline function _straddles(phi, c::CartesianIndex{N}, level::Float64) where {N}
    below = above = false
    for k in _corners(Val(N))
        v = phi[_corner(c, k)]
        below |= v < level
        above |= v > level
    end
    return below && above
end

@inline function _cellValue(phi, c::CartesianIndex{N}) where {N}
    s = 0.0
    for k in _corners(Val(N))
        s += phi[_corner(c, k)]
    end
    return s / 2^N
end

function _cellsWhere(phi, cell_dims::NTuple{N,Int}, level::Float64) where {N}
    mask = falses(cell_dims)
    for c in CartesianIndices(cell_dims)
        mask[c] = _straddles(phi, c, level)
    end
    return mask
end

# Any cell of `mask` within Chebyshev radius `r` of `c`, clipped to the grid.
function _anyInBox(mask, c::CartesianIndex{N}, r::Int) where {N}
    ranges = ntuple(d -> max(1, c[d] - r):min(size(mask, d), c[d] + r), N)
    for I in CartesianIndices(ranges)
        mask[I] && return true
    end
    return false
end

# Central differences on nodes, clamped at the box edges.
function _nodalGradient(phi::Array{Float64,N}, info::_GridMetrics{N}) where {N}
    dims, spacing = info.dims, info.spacing
    grad = Array{SVector{N,Float64}}(undef, dims)
    for I in CartesianIndices(dims)
        comps = ntuple(N) do d
            lo = max(I[d] - 1, 1)
            hi = min(I[d] + 1, dims[d])
            jlo = CartesianIndex(ntuple(k -> k == d ? lo : I[k], N))
            jhi = CartesianIndex(ntuple(k -> k == d ? hi : I[k], N))
            (phi[jhi] - phi[jlo]) / ((hi - lo) * spacing[d])
        end
        grad[I] = SVector{N,Float64}(comps)
    end
    return grad
end

# Unit normal at a cell, from its corner nodes; `nothing` when the gradient vanishes.
function _cellNormal(grad, c::CartesianIndex{N}) where {N}
    acc = zero(SVector{N,Float64})
    for k in _corners(Val(N))
        acc += grad[_corner(c, k)]
    end
    len = norm(acc)
    return len > 0 ? acc / len : nothing
end

@inline function _nodePosition(info::_GridMetrics{N}, I::CartesianIndex{N}) where {N}
    SVector{N,Float64}(ntuple(d -> info.origin[d] + (I[d] - 1) * info.spacing[d], N))
end

@inline function _cellCenter(info::_GridMetrics{N}, c::CartesianIndex{N}) where {N}
    SVector{N,Float64}(ntuple(d -> info.origin[d] + (c[d] - 0.5) * info.spacing[d], N))
end

# Points of the auxiliary level set implied by the cells of `mask` near `c`, Eq. (14).
function _auxPointsNear(phi, grad, info::_GridMetrics{N}, mask, c::CartesianIndex{N},
                        side::Float64, eps::Float64) where {N}
    cell_dims = size(mask)
    ranges = ntuple(d -> max(1, c[d] - _SEARCH):min(cell_dims[d], c[d] + _SEARCH), N)
    points = SVector{N,Float64}[]
    for p in CartesianIndices(ranges)
        mask[p] || continue
        n = _cellNormal(grad, p)
        n === nothing && continue
        # Projection of the auxiliary cell centre onto phi = side*eps, along +grad phi:
        # the signed distance to that isosurface is side*eps - phi, not phi - side*eps.
        push!(points, _cellCenter(info, p) + (side * eps - _cellValue(phi, p)) * n)
    end
    return points
end

@inline function _nearestDistance(x::SVector{N,Float64}, points) where {N}
    best = Inf
    for q in points
        d = hypot(ntuple(e -> x[e] - q[e], N)...)
        d < best && (best = d)
    end
    return best
end

# Nodes of the cells of `mask`, as a node-grid mask.
function _nodeMaskOfCells(mask, cell_dims::NTuple{N,Int}) where {N}
    nodes = falses(_nodeDims(cell_dims))
    for c in CartesianIndices(cell_dims)
        mask[c] || continue
        for k in _corners(Val(N))
            nodes[_corner(c, k)] = true
        end
    end
    return nodes
end

# Connected component of `seed` nodes over `allowed`, under face adjacency.
function _component(seed, allowed, ::Val{N}) where {N}
    member = falses(size(allowed))
    queue = CartesianIndex{N}[]
    for I in CartesianIndices(seed)
        (seed[I] && allowed[I]) || continue
        member[I] = true
        push!(queue, I)
    end
    head = 1
    while head <= length(queue)
        I = queue[head]
        head += 1
        for d in 1:N, s in (-1, 1)
            J = CartesianIndex(ntuple(k -> k == d ? I[k] + s : I[k], N))
            checkbounds(Bool, member, J) || continue
            (member[J] || !allowed[J]) && continue
            member[J] = true
            push!(queue, J)
        end
    end
    return member
end

# Value of the nearest node outside `island`, in face-adjacency hops, ties by traversal
# order (deterministic). One pass over the node grid per island.
function _nearestExteriorValues(phi, island, ::Val{N}) where {N}
    dims = size(island)
    value = Array{Float64}(undef, dims)
    seen = falses(dims)
    queue = CartesianIndex{N}[]
    for I in CartesianIndices(dims)
        island[I] && continue
        seen[I] = true
        value[I] = phi[I]
        push!(queue, I)
    end
    head = 1
    while head <= length(queue)
        I = queue[head]
        head += 1
        for d in 1:N, s in (-1, 1)
            J = CartesianIndex(ntuple(k -> k == d ? I[k] + s : I[k], N))
            checkbounds(Bool, seen, J) || continue
            seen[J] && continue
            seen[J] = true
            value[J] = value[I]
            push!(queue, J)
        end
    end
    return value
end

"""
    filterUnresolvedFeatures!(field, level)

Let a 2D or 3D field store an approximate signed-distance function. The method removes
geometrically unresolved structures of the interface defined by its zero level set on one
level of `field`. The nodal values are mutated in place and the field is returned.

The algorithm is an adaptation of Luo, Hu & Adams (2016), "Efficient formulation of scale
separation for multi-scale modeling of interfacial flows", Journal of Computational Physics
308, 411–420, https://doi.org/10.1016/j.jcp.2015.11.044.
"""
function filterUnresolvedFeatures!(field::Field, level::Symbol)
    phi = _levelvalues(field, level)
    grid = getproperty(field, :discretization).levels[level]
    info = _gridMetrics(grid)
    phi isa Array{Float64,length(info.dims)} || throw(ArgumentError(
        "filterUnresolvedFeatures! requires a contiguous Array{Float64,N}; materialize views first."))
    _checklevel(phi, nothing, info, "Level-set Field values")
    cell_dims = _cellDims(info.dims)
    _filter!(phi, info, cell_dims)
    return field
end

function _filter!(phi::Array{Float64,N}, info::_GridMetrics{N},
                  cell_dims::NTuple{N,Int}) where {N}
    eps = _EPS_FACTOR * minimum(info.spacing)
    S = _cellsWhere(phi, cell_dims, 0.0)
    any(S) || return phi                                   # single-signed: nothing to filter
    maximum(abs, phi) >= eps || throw(ArgumentError(
        "filterUnresolvedFeatures! requires an approximately signed-distance level set: " *
        "|phi| never reaches eps = $eps, so no auxiliary level set exists and the whole " *
        "phase would be classified as unresolved. Reinitialize the field first."))
    P = _cellsWhere(phi, cell_dims, +eps)
    Nset = _cellsWhere(phi, cell_dims, -eps)
    grad = _nodalGradient(phi, info)

    # Classification. `side[c]` is the auxiliary side cell c is re-distanced from
    # (Eq. 12/13: the side its neighbourhood misses); `missing` are the cells for which no
    # auxiliary cell of that side exists in the search region.
    resolved = falses(cell_dims)
    side = zeros(Int8, cell_dims)
    missing = falses(cell_dims)
    missing_side = zeros(Int8, cell_dims)
    for c in CartesianIndices(cell_dims)
        S[c] || continue
        hasP = _anyInBox(P, c, _NEIGHBOURHOOD)
        hasN = _anyInBox(Nset, c, _NEIGHBOURHOOD)
        if hasP && hasN
            resolved[c] = true
        elseif hasP != hasN
            side[c] = hasP ? Int8(-1) : Int8(1)
        else
            # Missing both: the cell's own phase says which auxiliary level set it cannot
            # reach. This case is undefined by the paper.
            missing[c] = true
            missing_side[c] = _cellValue(phi, c) >= 0 ? Int8(1) : Int8(-1)
        end
    end
    protected = _nodeMaskOfCells(resolved, cell_dims)

    values = copy(phi)
    touched = falses(size(phi))
    for c in CartesianIndices(cell_dims)
        side[c] == 0 && continue
        s = Float64(side[c])
        points = _auxPointsNear(phi, grad, info, s > 0 ? P : Nset, c, s, eps)
        if isempty(points)
            missing[c] = true
            missing_side[c] = side[c]
            continue
        end
        for k in _corners(Val(N))
            I = _corner(c, k)
            protected[I] && continue             # never move a resolved cell's node
            d = _nearestDistance(_nodePosition(info, I), points)
            values[I] = s * (eps - d)
            touched[I] = true
        end
    end

    for s in (Int8(1), Int8(-1))
        seeds = falses(cell_dims)
        for c in CartesianIndices(cell_dims)
            seeds[c] = missing[c] && missing_side[c] == s
        end
        any(seeds) || continue
        _fallback!(values, touched, phi, info, cell_dims, seeds, Float64(s),
                   protected, eps)
    end

    for I in CartesianIndices(phi)
        touched[I] && (phi[I] = values[I])
    end
    return phi
end

# Interpolation fallback for an island that no auxiliary level set reaches (project
# adaptation, not part of Luo et al.). The whole connected same-sign region is replaced by
# the nearest exterior sample, which removes the structure without creating a crossing.
function _fallback!(values, touched, phi, info::_GridMetrics{N}, cell_dims, seeds,
                    s::Float64, protected, eps) where {N}
    inphase = BitArray(map(v -> s * v > 0, phi))
    island = _component(_nodeMaskOfCells(seeds, cell_dims), inphase, Val(N))
    any(island) || return nothing
    for I in CartesianIndices(island)
        island[I] || continue
        protected[I] && return nothing        # touches resolved geometry: not ours to delete
        s * phi[I] >= eps && return nothing   # reaches the auxiliary set elsewhere
    end
    exterior = _nearestExteriorValues(phi, island, Val(N))
    for I in CartesianIndices(island)
        island[I] || continue
        values[I] = exterior[I]
        touched[I] = true
    end
    return nothing
end

end # module UnresolvedFeatures
