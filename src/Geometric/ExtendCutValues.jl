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
    extendCutValues(values, cut::Cut, grid::Grid; band, fillvalue) -> Array
    extendCutValues!(dest, values, cut, grid; band, fillvalue) -> dest

Extend data defined **per interface facet** of `cut` onto the nodes of `grid` by
nearest-facet projection, writing `dest` (the first form allocates it at the node shape of
`grid`):

    u_ext(x) = values[f*(x)],   f*(x) = argmin_{f ∈ Γ_h} d(x, f)

`d` is the Euclidean distance to the facet **as a set** — the segment in 2D, the triangle in
3D, not the centroid — and `Γ_h = cut.boundary` is the piecewise-linear interface of the
snapshot.

`values` is indexed in the cut's own facet order: `values[f]` belongs to facet `f` of
`cut.boundary`, the facet whose parent cell is `cut.boundary.facet_to_bgcell[f]`. A producer
of per-facet data emits it in that order and needs no mapping of its own. Scalars,
`SVector`s and any other element type are extended alike; `dest` has the shape of the grid
node array, and `vec(dest)` is the linear view that transport operators expect.

`band` is the distance around the interface where the extension is defined (default three
times the largest node spacing; `Inf` extends over the whole grid). Nodes farther than
`band` from every facet keep `fillvalue` (default `zero(eltype(values))`). Where two facets
are exactly equidistant, the lowest facet index wins.

**This is the piecewise-constant nearest-facet extension, not a smooth normal extension.**
Each node carries one facet's value, so on a curved interface the nearest-facet regions meet
in jumps of the facet values and `u_ext` is discontinuous across the medial axis, even when
the interface field is smooth. That is the chosen approximation, not an implementation
defect: refine, or use a normal-projection method such as `extendClosestPoint!`, when a
continuous extension is required.

Cost is `O(nfacets × nodes within band of a facet)`: a facet group visits only its own
neighbourhood, and two cheap bounds (the group's box, then the facet's own box) reject nearly
every candidate before the exact distance is computed, so the whole grid may be large and the
band narrow. Measured at `band = 3h`: 1.7 ms with 2048 facets on 33³, 7 ms with 8192 on 65³,
32 ms with 32 768 on 129³ (tiled plane); 12/18/42 ms for real sphere cuts with 6288/9384/
24 264 facets. In the mean-curvature-flow validation the extension is ~2 % of an evolution
step, the facet curvature it transfers being an order of magnitude more expensive. The cost
grows with the **facet count**, so a much finer interface is the case to watch: a node-driven
search over a cell→facet index would make it `O(band nodes × facets of a few cells)` and
share the index across several fields; it is not offered yet.

`Cut` and `Grid` must carry the same `mesh_id`, and `grid` must be uniformly spaced and 2D or
3D, the same requirements as `reinitialize!`. `values` must hold one entry per facet of
`cut.boundary` and `dest` must have the grid's node shape. A cut with no interface
(`isempty(cut.boundary.facet_to_bgcell)`) leaves `dest` uniformly `fillvalue`.

# Example

```julia
snapshot = cut(geometry, :fine)
κ = interfaceCurvature(snapshot.boundary; ...)      # one value per facet, e.g. Curvature.jl
velocity = extendCutValues(-μ .* κ .* snapshot.boundary.normals, snapshot, grid; band=3h)
eulerianTransport!(geometry, vec(velocity), dt)
```
"""
function extendCutValues(values::AbstractVector, cut::Cut, grid::Grid;
                         band=nothing, fillvalue=zero(eltype(values)))
    dims = _cutGridInfo(cut, grid).dims
    dest = Array{eltype(values)}(undef, dims)
    return extendCutValues!(dest, values, cut, grid; band=band, fillvalue=fillvalue)
end

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
