module Extension

using StaticArrays: SVector
using LinearAlgebra: dot, cross
using ..Geometric: _cutGridInfo, _getInterpolator, _sameid
using ....EvolvingDomains: Cut, Field, Grid, _levelgrid, _levelvalues

# =============================================================================
# Closest-point geometry on the reconstructed interface Γ_h.
#
# `_closestPointOnFacet` returns the point of the facet (segment in 2D, triangle in 3D) that
# minimises the Euclidean distance to `p`. Distances are derived from it, so the search and
# the evaluation can never disagree about which point is closest.

@inline function _closestPointOnSegment(p::SVector{N,Float64}, a::SVector{N,Float64},
                                        b::SVector{N,Float64}) where {N}
    e = b - a
    len2 = dot(e, e)
    len2 > 0 || return a                     # a degenerate facet is its first vertex
    return a + clamp(dot(p - a, e) / len2, 0.0, 1.0) * e
end

@inline _closestPointOnFacet(p::SVector{2,Float64}, v::NTuple{2,SVector{2,Float64}}) =
    _closestPointOnSegment(p, v[1], v[2])

@inline function _closestPointOnFacet(p::SVector{3,Float64},
                                      v::NTuple{3,SVector{3,Float64}})
    a, b, c = v
    e1, e2 = b - a, c - a
    n = cross(e1, e2)
    n2 = dot(n, n)                           # also d00 d11 - d01², the barycentric denominator
    if n2 > 0
        q = p - (dot(p - a, n) / n2) * n                                  # plane projection
        w = q - a
        d00, d01, d11 = dot(e1, e1), dot(e1, e2), dot(e2, e2)
        s = (d11 * dot(w, e1) - d01 * dot(w, e2)) / n2
        t = (d00 * dot(w, e2) - d01 * dot(w, e1)) / n2
        s >= 0 && t >= 0 && s + t <= 1 && return q                        # interior projection
    end
    qab = _closestPointOnSegment(p, a, b)
    qbc = _closestPointOnSegment(p, b, c)
    qca = _closestPointOnSegment(p, c, a)
    dab, dbc, dca = dot(p - qab, p - qab), dot(p - qbc, p - qbc), dot(p - qca, p - qca)
    return dab <= dbc ? (dab <= dca ? qab : qca) : (dbc <= dca ? qbc : qca)
end

@inline function _facetDistanceSquared(p, v)
    q = p - _closestPointOnFacet(p, v)
    return dot(q, q)
end

# Node index just outside the dilated facet box on one axis, clamped into the grid (`±Inf`
# band keeps the whole axis). `rounding` is `floor` for the low bound and `ceil` for the high.
@inline function _nodeBound(x, origin::Float64, h::Float64, dims::Int, rounding)
    t = (Float64(x) - origin) / h
    isfinite(t) || return signbit(t) ? 1 : dims
    # Saturate before converting: a band far larger than the box must clamp, not overflow Int.
    return clamp(rounding(Int, clamp(t, -2.0 * dims, 2.0 * dims)) + 1, 1, dims)
end

# Axis-aligned box of a facet group's vertices, and the node index ranges covering it dilated
# by `band`, with one node of slack on each side (the exact distance test filters the corners).
function _facetBox(vertices::AbstractVector{<:SVector{N,Float64}}, band::Float64,
                   origin, spacing, dims) where {N}
    box = ntuple(N) do d
        (minimum(v[d] for v in vertices), maximum(v[d] for v in vertices))
    end
    ranges = ntuple(N) do d
        o, h, n = origin[d], spacing[d], dims[d]
        low = _nodeBound(box[d][1] - band, o, h, n, floor)
        high = _nodeBound(box[d][2] + band, o, h, n, ceil)
        low:high
    end
    return box, ranges
end

# Squared distance from `x` to the facet's axis-aligned box: a lower bound on the distance to
# any point of the facet, so a node that fails this test cannot be nearest to it.
@inline function _boxDistanceSquared(x::SVector{N,Float64}, box::NTuple{N}) where {N}
    d2 = 0.0
    for d in 1:N
        low, high = box[d]
        gap = x[d] < low ? low - x[d] : x[d] > high ? x[d] - high : 0.0
        d2 += gap * gap
    end
    return d2
end

@inline function _facetVertices(surface, facet::Int, ::Val{N}) where {N}
    first = Int(surface.point_pointers[facet])
    return ntuple(q -> SVector{N,Float64}(
        surface.physical_points[surface.point_ids[first + q - 1]]), N)
end

@inline _nodePoint(I, origin, spacing) =
    SVector(ntuple(d -> origin[d] + (I[d] - 1) * spacing[d], length(origin)))

# =============================================================================
# Banded nearest-facet search.
#
# Facets are grouped by parent background cell (`facet_to_bgcell` is ascending by contract),
# so a group shares one node neighbourhood. Each group is visited as a box: a node outside
# the box cannot be nearest to any of its facets. Two passes, `b/2` then `b`, pin the answer
# for the nodes closest to the interface first, so the second pass rejects nearly every facet
# with one box test. A node's winner only ever changes to a strictly closer facet, visited in
# ascending facet order, so the result is deterministic and independent of the pass split.

# `accept(node, facet, first_touch)` is called whenever a node's nearest facet improves, so a
# caller that only needs the winning facet writes its output inline and allocates nothing but
# the distance buffer; a caller that must evaluate the source at the projected point records
# the winners and evaluates afterwards.
function _nearestFacetSearch(accept, surface, info, b::Float64)
    dims, origin, spacing = info.dims, info.origin, info.spacing
    N = length(dims)
    nearest = fill(Inf, dims)
    coords = ntuple(d -> Float64[origin[d] + (i - 1) * spacing[d] for i in 1:dims[d]], N)
    cell_of = surface.facet_to_bgcell
    nfacets = length(cell_of)
    group = Vector{NTuple{N,SVector{N,Float64}}}()
    flat = SVector{N,Float64}[]
    boxes = Vector{NTuple{N,Tuple{Float64,Float64}}}()

    for radius in (isfinite(b) ? (0.5 * b, b) : (b,))
        r2 = radius * radius
        facet = 1
        while facet <= nfacets
            last = facet
            while last < nfacets && cell_of[last + 1] == cell_of[facet]
                last += 1
            end
            resize!(group, last - facet + 1)
            resize!(boxes, last - facet + 1)
            empty!(flat)
            for q in eachindex(group)
                group[q] = _facetVertices(surface, facet + q - 1, Val(N))
                boxes[q] = ntuple(d -> (minimum(v[d] for v in group[q]),
                                        maximum(v[d] for v in group[q])), N)
                append!(flat, group[q])
            end
            box, ranges = _facetBox(flat, radius, origin, spacing, dims)
            low, high = box[N]
            oN, hN, nN = origin[N], spacing[N], dims[N]
            for I in CartesianIndices(ranges[1:N - 1])
                gap2 = 0.0
                for d in 1:N - 1
                    xd = coords[d][I[d]]
                    lo, hi = box[d]
                    g = xd < lo ? lo - xd : xd > hi ? xd - hi : 0.0
                    gap2 += g * g
                end
                gap2 <= r2 || continue
                remaining = sqrt(r2 - gap2)
                k1 = _nodeBound(low - remaining, oN, hN, nN, floor)
                k2 = _nodeBound(high + remaining, oN, hN, nN, ceil)
                @inbounds for k in k1:k2
                    J = CartesianIndex(ntuple(d -> d == N ? k : I[d], N))
                    x = SVector{N,Float64}(ntuple(d -> coords[d][J[d]], N))
                    _boxDistanceSquared(x, box) < nearest[J] || continue
                    for q in eachindex(group)
                        _boxDistanceSquared(x, boxes[q]) < nearest[J] || continue
                        distance = _facetDistanceSquared(x, group[q])
                        distance <= r2 || continue
                        if distance < nearest[J]
                            accept(J, facet + q - 1, isinf(nearest[J]))
                            nearest[J] = distance
                        end
                    end
                end
            end
            facet = last + 1
        end
    end
    return nothing
end

_defaultBand(info) = 3 * maximum(info.spacing)

"""
    extend!(destination, source, cut::Cut, grid::Grid; band = nothing)
    extend!(destination::Field, source::Field, cut::Cut; level_id = cut.level_id, band = nothing)

Extend interface quantities off the reconstructed interface `Γ_h` of `cut` by Euclidean
closest-point projection: for a grid node `x`,

    extend!(…)[x] = u(P(x)),    P(x) = argmin_{p ∈ Γ_h} ‖x − p‖₂,

where `P(x)` is the closest point on the actual segments (2D) or triangles (3D) of the cut.
Supports per-facet and nodal source values.

# References

- D. Adalsteinsson and J. A. Sethian, *The Fast Construction of Extension Velocities in Level
  Set Methods*, J. Comput. Phys. 148 (1999) 2–22,
  [doi:10.1006/jcph.1998.6090](https://doi.org/10.1006/jcph.1998.6090) — the normal-constant
  extension principle that this operator realises; the fast-marching algorithm of the paper is
  not what is implemented.
- R. I. Saye and J. A. Sethian, *A Review of Level Set Methods to Model Interfaces Moving under
  Complex Physics*, Handbook of Numerical Analysis 21 (2020), §4.1.3,
  [doi:10.1016/bs.hna.2019.07.003](https://doi.org/10.1016/bs.hna.2019.07.003) — closest-point
  evaluation of interface quantities, which is the geometric construction implemented here on
  the piecewise-linear cut.
"""
function extend!(dest::AbstractArray, facet_values::AbstractVector, cut::Cut, grid::Grid;
                 band=nothing)
    info = _cutGridInfo(cut, grid)
    dims = info.dims
    Base.require_one_based_indexing(dest, facet_values)
    size(dest) == dims || throw(DimensionMismatch(
        "extend!: destination has size $(size(dest)); grid requires $(dims)."))
    values = Base.mightalias(dest, facet_values) ? copy(facet_values) : facet_values
    nfacets = length(cut.boundary.facet_to_bgcell)
    length(values) == nfacets || throw(DimensionMismatch(
        "extend!: source holds $(length(values)) values but the cut boundary has $nfacets facets."))
    b = band === nothing ? _defaultBand(info) : Float64(band)
    b >= 0 || throw(ArgumentError("extend!: band must be non-negative; got $band."))
    nfacets == 0 && return dest
    _nearestFacetSearch(cut.boundary, info, b) do J, facet, _
        dest[J] = values[facet]
    end
    return dest
end

function extend!(dest::AbstractArray, source::AbstractArray{S,N}, cut::Cut{T,N},
                 grid::Grid{N}; band=nothing) where {S,T,N}
    info = _cutGridInfo(cut, grid)
    dims = info.dims
    Base.require_one_based_indexing(dest, source)
    size(dest) == dims || throw(DimensionMismatch(
        "extend!: destination has size $(size(dest)); grid requires $(dims)."))
    size(source) == dims || throw(DimensionMismatch(
        "extend!: nodal source has size $(size(source)); grid requires $(dims)."))
    b = band === nothing ? _defaultBand(info) : Float64(band)
    b >= 0 || throw(ArgumentError("extend!: band must be non-negative; got $band."))
    surface = cut.boundary
    isempty(surface.facet_to_bgcell) && return dest
    winner = zeros(Int32, dims)
    touched = CartesianIndex{N}[]
    _nearestFacetSearch(surface, info, b) do J, facet, first
        first && push!(touched, J)
        winner[J] = Int32(facet)
    end
    values = Base.mightalias(dest, source) ? copy(source) : source
    itp = _getInterpolator(vec(values), info)    # also copies: writes cannot poison reads
    for I in touched
        q = _closestPointOnFacet(_nodePoint(I, info.origin, info.spacing),
                                 _facetVertices(surface, Int(winner[I]), Val(N)))
        dest[I] = itp(q...)
    end
    return dest
end

function extend!(dest::Field, source::Field, cut::Cut;
                 level_id::Symbol=cut.level_id, band=nothing)
    dd, sd = getproperty(dest, :discretization), getproperty(source, :discretization)
    grid = _levelgrid(dd, level_id)
    source_grid = _levelgrid(sd, level_id)
    _sameid(grid, source_grid) || throw(ArgumentError(
        "extend!: source and destination Fields belong to different grids."))
    grid.axes == source_grid.axes || throw(ArgumentError(
        "extend!: source and destination Fields do not share their axes."))
    dest_values = _levelvalues(dest, level_id)
    source_values = _levelvalues(source, level_id)
    extend!(dest_values, source_values, cut, grid; band=band)
    return dest
end

end # module Extension
