module ExtendCutValues

using StaticArrays: SVector
using LinearAlgebra: dot, cross
using ..Geometric: _cutGridInfo
using ....EvolvingDomains: Cut, Grid

# Squared distance from `p` to the segment `[a, b]`. A zero-length segment is a point, and
# a hand-built snapshot is only structurally valid, so the degenerate case is handled.
@inline function _segmentDistanceSquared(p::SVector{N,Float64}, a::SVector{N,Float64},
                                        b::SVector{N,Float64}) where {N}
    e = b - a
    len2 = dot(e, e)
    t = len2 > 0 ? clamp(dot(p - a, e) / len2, 0.0, 1.0) : 0.0
    q = p - (a + t * e)
    return dot(q, q)
end

# Squared distance from `p` to the facet: the segment `[a, b]` in 2D, the triangle `[a, b, c]`
# in 3D. The closest point of a triangle is its plane projection when that projection is
# interior and an edge point otherwise; a degenerate triangle has no plane and is all edge.
@inline _facetDistanceSquared(p::SVector{2,Float64}, v::NTuple{2,SVector{2,Float64}}) =
    _segmentDistanceSquared(p, v[1], v[2])

@inline function _facetDistanceSquared(p::SVector{3,Float64}, v::NTuple{3,SVector{3,Float64}})
    a, b, c = v
    e1, e2 = b - a, c - a
    n = cross(e1, e2)
    n2 = dot(n, n)                       # also d00 d11 - d01², the barycentric denominator
    if n2 > 0
        q = p - (dot(p - a, n) / n2) * n                                    # plane projection
        w = q - a
        d00, d01, d11 = dot(e1, e1), dot(e1, e2), dot(e2, e2)
        s = (d11 * dot(w, e1) - d01 * dot(w, e2)) / n2
        t = (d00 * dot(w, e2) - d01 * dot(w, e1)) / n2
        s >= 0 && t >= 0 && s + t <= 1 && return dot(p - q, p - q)
    end
    return min(_segmentDistanceSquared(p, a, b), _segmentDistanceSquared(p, b, c),
               _segmentDistanceSquared(p, c, a))
end

# Node index just outside the dilated facet box on one axis, clamped into the grid (`±Inf`
# band keeps the whole axis). `rounding` is `floor` for the low bound and `ceil` for the high.
@inline function _nodeBound(x, origin::Float64, h::Float64, dims::Int, rounding)
    isfinite(x) || return signbit(x) ? 1 : dims
    return clamp(rounding(Int, (x - origin) / h) + 1, 1, dims)
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

"""
    extendCutValues(values::AbstractVector, cut::Cut, grid::Grid; band=nothing, fillvalue=zero(eltype(values)))

Inputs:
- `values`: one value per facet of `cut.boundary`, in facet order.
- `cut`: discrete interface snapshot; `grid` must have the same mesh identity.
- `grid`: uniform 2D or 3D grid on which node values are produced.
- `band`: nonnegative physical distance from the interface to extend; defaults to three times the largest node spacing, and `Inf` covers the grid.
- `fillvalue`: value assigned to nodes outside the band; defaults to `zero(eltype(values))`.

Returns:
- A node-shaped array with `u_ext(x) = values[f*(x)]`, where `f* = argmin d(x,f)` uses Euclidean distance to the facet as a set. Uncovered nodes contain `fillvalue`; equal-distance ties select the lowest facet index.

Requires: `cut` and `grid` share mesh identity; the grid is uniform and 2D or 3D; `values` has one entry per facet.
Use: Extend per-facet interface data as a piecewise-constant nearest-facet field, which is discontinuous across the medial axis when the facet values differ.

Any element type is extended alike (scalars, `SVector`); `vec(dest)` is the linear view transport operators expect.
"""
function extendCutValues(values::AbstractVector, cut::Cut, grid::Grid;
                         band=nothing, fillvalue=zero(eltype(values)))
    dims = _cutGridInfo(cut, grid).dims
    dest = Array{eltype(values)}(undef, dims)
    return extendCutValues!(dest, values, cut, grid; band=band, fillvalue=fillvalue)
end

"""
    extendCutValues!(dest::AbstractArray, values::AbstractVector, cut::Cut{T,N}, grid::Grid; band=nothing, fillvalue=zero(eltype(values))) where {T,N}

Inputs:
- `dest`: node array to mutate; its shape must match `grid`.
- `values`: one value per facet of `cut.boundary`, in facet order.
- `cut`: discrete interface snapshot; `grid` must have the same mesh identity.
- `grid`: uniform 2D or 3D grid.
- `band`: nonnegative physical distance from the interface; defaults to three times the largest node spacing, and `Inf` covers the grid.
- `fillvalue`: value assigned outside the band; defaults to `zero(eltype(values))`.

Returns:
- Mutated `dest`, with `u_ext(x) = values[f*(x)]` for the nearest facet as a geometric set; uncovered nodes contain `fillvalue`, and ties select the lowest facet index.

Requires: matching `cut` and `grid` mesh identity, uniform 2D/3D grid, one `values` entry per facet, and destination shape equal to the grid node shape.
Use: Extend per-facet interface data as a piecewise-constant nearest-facet field.
"""
function extendCutValues!(dest::AbstractArray, values::AbstractVector, cut::Cut{T,N},
                          grid::Grid; band=nothing,
                          fillvalue=zero(eltype(values))) where {T,N}
    info = _cutGridInfo(cut, grid)          # same mesh identity, uniform, 2D or 3D
    dims, origin, spacing = info.dims, info.origin, info.spacing
    size(dest) == dims || throw(DimensionMismatch(
        "extendCutValues!: destination has size $(size(dest)); grid requires $(dims)."))
    surface = cut.boundary
    nfacets = length(surface.facet_to_bgcell)
    length(values) == nfacets || throw(DimensionMismatch(
        "extendCutValues!: values hold $(length(values)) entries but the cut boundary has " *
        "$nfacets facets."))
    b = band === nothing ? 3 * maximum(spacing) : Float64(band)
    b >= 0 || throw(ArgumentError("extendCutValues!: band must be non-negative."))
    b2 = b * b
    fill!(dest, fillvalue)
    nfacets == 0 && return dest

    # Squared distance to the nearest facet found so far, per node. Facets are visited in
    # ascending order and comparisons are strict, which is what breaks a tie towards the
    # lowest facet index.
    # shortcut: one Float64 per node per call; cache it only if extension becomes hot.
    nearest = fill(Inf, dims)
    coords = ntuple(d -> Float64[origin[d] + (i - 1) * spacing[d] for i in 1:dims[d]], N)
    # Facets are grouped by parent cell (`facet_to_bgcell` is ascending, a `CutSurface`
    # guarantee), and a facet lies inside its parent cell, so a group's facets share one node
    # neighbourhood. A hand-built cut that breaks that geometric guarantee stays correct, but
    # pays for the wider box of its group.
    cell_of = surface.facet_to_bgcell
    group = Vector{NTuple{N,SVector{N,Float64}}}()
    flat = SVector{N,Float64}[]
    boxes = Vector{NTuple{N,Tuple{Float64,Float64}}}()
    # Two passes: the near shell first, which pins the nearest facet for the nodes closest to
    # the interface, then the full band, where that answer rejects nearly every facet with one
    # box test. The result does not depend on the split: a node's value is only ever replaced
    # by a strictly closer facet, and facets are visited in ascending index order within and
    # across the groups each pass visits.
    for radius in (isfinite(b) ? (0.5 * b, b) : (b,))
        r2 = radius * radius
        facet = 1
        while facet <= nfacets
            last = facet
            while last < nfacets && cell_of[last + 1] == cell_of[facet]
                last += 1
            end
            # One cell's facets share their neighbourhood: the node box, the cheap rejection
            # tests and each node's coordinate are computed once for the group.
            resize!(group, last - facet + 1)
            resize!(boxes, last - facet + 1)
            empty!(flat)
            for q in eachindex(group)
                first = Int(surface.point_pointers[facet + q - 1])
                group[q] = ntuple(d -> SVector{N,Float64}(
                    surface.physical_points[surface.point_ids[first + d - 1]]), N)
                boxes[q] = ntuple(d -> (minimum(v[d] for v in group[q]),
                                        maximum(v[d] for v in group[q])), N)
                append!(flat, group[q])
            end
            box, ranges = _facetBox(flat, radius, origin, spacing, dims)
            # On the last axis only the nodes within the rest of the band are visited, which
            # turns the dilated box into a ball: most of a box is outside the ball.
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
                    _boxDistanceSquared(x, box) < nearest[J] || continue  # cheap cull: same winner
                    for q in eachindex(group)
                        _boxDistanceSquared(x, boxes[q]) < nearest[J] || continue
                        distance = _facetDistanceSquared(x, group[q])
                        distance <= b2 || continue
                        if distance < nearest[J]
                            nearest[J] = distance
                            dest[J] = values[facet + q - 1]
                        end
                    end
                end
            end
            facet = last + 1
        end
    end
    return dest
end

end # module
