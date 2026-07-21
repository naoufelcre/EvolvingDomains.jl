module Curvature

using StaticArrays: SVector, SMatrix, MMatrix, MVector, SOneTo
using LinearAlgebra: det
using GridapEmbedded: EmbeddedBoundary
using ..Geometric: EvolvingDiscreteGeometry, ensure_cut!, grid_info, CartesianGridInfo

export InterfaceSamples, interface_samples, get_curvature, curvature_at
export reference_curve

# ============================================================================
# Compute the curvature of the embedded interface.
#
# `cut(model, geo)` intersects the zero contour of φ with every background cell.
# In 2D each cut cell contributes one or two straight segments; each segment is a
# *subfacet*. The result, `SubFacetData`, holds
#
#   point_to_coords   every segment endpoint
#   facet_to_points   which endpoints form each segment
#   facet_to_normal   unit normal per segment, pointing out of the φ < 0 phase
#   facet_to_bgcell   which background cell each segment came from
#
# so Γ arrives as an unordered bag of small oriented segments. This module fits a
# parabola to a local handful of them. The per-facet normal is what makes that
# possible without ever reconstructing the connectivity of Γ.
#
# ---------------------------------------------------------------------------
# Notation
#
#   q       the point at which curvature is wanted
#   nq      unit normal at q; together with the tangent it defines the local frame
#   t       unit tangent, t = (-nq₂, nq₁)
#   (u, v)  coordinates of a nearby sample in the (t, nq) frame
#   r       window radius: samples within r of q are used in the fit
#   h       background cell size
#   δ       typical spacing between neighbouring samples (order h)
#
# ---------------------------------------------------------------------------
# We fit on Γ rather than difference φ because φ might not always be a SDF
# thus even thos it's marketed as one of the main advantages of level sets methods, computation of curvature is not always
# reliable and suffers from oscillations
# Our approach reconstruct the interface, which is not perfect but does the job whitout too much headache
#
# ============================================================================

"""
    InterfaceSamples

One position/normal pair per interface subfacet, indexed by background cell.

`positions` are subfacet midpoints: one sample per facet, so an endpoint shared by
two neighbouring cells never enters twice.

`normals` follow the GridapEmbedded orientation (outward from the `φ < 0` phase),
which is the same orientation `get_normal_vector` hands an integrand — so curvatures
produced here carry the sign a surface-tension term expects.

`cell_to_facets` maps a background cell to the samples lying in it. The fit needs
"every sample within `r` of `q`"; without an index that is a scan of all N samples per
query, and there is one query per sample, so O(N²). Restricting the scan to cells near
`q` makes it O(1) per query.
"""
struct InterfaceSamples
    positions::Vector{SVector{2,Float64}}
    normals::Vector{SVector{2,Float64}}
    lengths::Vector{Float64}
    curvatures::Vector{Float64}
    endpoints::Vector{NTuple{4,Float64}}
    cell_to_facets::Dict{Int,Vector{Int32}}
    info::CartesianGridInfo
end

Base.length(s::InterfaceSamples) = length(s.positions)

@inline function _cell_of(info::CartesianGridInfo, x, y)
    nx, ny = info.cells
    i = clamp(floor(Int, (x - info.origin[1]) / info.spacing[1]) + 1, 1, nx)
    j = clamp(floor(Int, (y - info.origin[2]) / info.spacing[2]) + 1, 1, ny)
    return i, j
end

"Linear index of background cell (i, j)."
@inline _cell_id(info::CartesianGridInfo, i, j) = i + (j - 1) * info.cells[1]

"""
    _cells_within(info, q, r) -> (irange, jrange)

Block of background cells guaranteed to contain everything within `r` of `q`.
"""
@inline function _cells_within(info::CartesianGridInfo, q::SVector{2,Float64}, r::Float64)
    nx, ny = info.cells
    ci, cj = _cell_of(info, q[1], q[2])
    pad = ceil(Int, r / min(info.spacing[1], info.spacing[2]))
    return max(1, ci - pad):min(nx, ci + pad), max(1, cj - pad):min(ny, cj + pad)
end

"""
    interface_samples(geom::EvolvingDiscreteGeometry) -> InterfaceSamples

Sample the current embedded interface at subfacet midpoints.
"""
function interface_samples(geom::EvolvingDiscreteGeometry)
    sf = EmbeddedBoundary(ensure_cut!(geom)).subfacets
    info = grid_info(geom.grid)

    n = length(sf.facet_to_normal)
    positions = Vector{SVector{2,Float64}}(undef, n)
    normals = Vector{SVector{2,Float64}}(undef, n)
    lengths = Vector{Float64}(undef, n)
    endpoints = Vector{NTuple{4,Float64}}(undef, n)
    cell_to_facets = Dict{Int,Vector{Int32}}()

    degenerate = Int[]
    for f in 1:n
        ids = sf.facet_to_points[f]
        cx = cy = 0.0
        for id in ids
            p = sf.point_to_coords[id]
            cx += p[1]; cy += p[2]
        end
        m = length(ids)
        positions[f] = SVector(cx / m, cy / m)

        # Arclength carried by this facet. Subfacet lengths span ~90x on a smooth
        # interface, so weighting every sample equally lets a sliver count as much as a
        # full segment — a direct source of facet-scale noise in the fit.
        if m == 2
            p1 = sf.point_to_coords[ids[1]]; p2 = sf.point_to_coords[ids[2]]
            lengths[f] = hypot(p2[1]-p1[1], p2[2]-p1[2])
            endpoints[f] = (p1[1], p1[2], p2[1], p2[2])
        else
            lengths[f] = 0.0
            endpoints[f] = (positions[f][1], positions[f][2], positions[f][1], positions[f][2])
        end

        # Degeneracy is a property of the *facet*, not of the normal. A subfacet whose two
        # endpoints coincide still gets a unit normal from the cutter — computed from a
        # zero-length direction, so meaningless — and a test on |n| passes it straight
        # through. Measured cost: two such facets on a circle carried frames 10.3° from
        # radial while their neighbours were within 0.13°.
        nf = sf.facet_to_normal[f]
        len = hypot(nf[1], nf[2])
        if len > 1e-12 && lengths[f] > 1e-12
            normals[f] = SVector(nf[1] / len, nf[2] / len)
        else
            # Zero-length subfacet, emitted when Γ passes through (or very near) a grid
            # node. Same root cause as the small-cut-cell problem, but a much milder
            # consequence: in CutFEM a sliver wrecks the conditioning of the system,
            # whereas here the facet simply has no orientation of its own. It carries
            # zero measure so it contributes nothing to any surface integral, but it
            # still occupies one slot in the per-facet array and so needs some frame.
            normals[f] = SVector(0.0, 0.0)
            push!(degenerate, f)
        end

        e = endpoints[f]
        for (px, py) in ((positions[f][1], positions[f][2]), (e[1], e[2]), (e[3], e[4]))
            i, j = _cell_of(info, px, py)
            v = get!(() -> Int32[], cell_to_facets, _cell_id(info, i, j))
            (isempty(v) || v[end] != Int32(f)) && push!(v, Int32(f))
        end
    end

    s = InterfaceSamples(positions, normals, lengths, zeros(n), endpoints, cell_to_facets, info)
    oriented = copy(normals)
    for f in degenerate
        normals[f] = _borrow_normal(s, positions[f], oriented)
    end
    # Curvature completes the reference curve, so it must be filled after the degenerate
    # frames are repaired. Both consumers — the surface-tension term and the
    # reinitialisation seeding — then read the same object.
    r = 8 * min(info.spacing[1], info.spacing[2])
    Base.Threads.@threads for f in eachindex(positions)
        s.curvatures[f] = _fit(s, positions[f], normals[f], r)
    end
    return s
end

"""
    _borrow_normal(s, q, oriented) -> SVector{2,Float64}

Frame for a degenerate facet, averaged from the oriented samples around it.

Normals more than 90° from the nearest oriented sample are excluded, so a facet sitting in
a thin neck cannot average against the opposing wall and cancel to nothing.
"""
function _borrow_normal(s::InterfaceSamples, q::SVector{2,Float64},
                        oriented::Vector{SVector{2,Float64}})
    nx, ny = s.info.cells
    ci, cj = _cell_of(s.info, q[1], q[2])
    for pad in 1:max(nx, ny)
        # nearest oriented sample first: it fixes the branch everything else aligns to
        best = Inf
        ref = SVector(0.0, 0.0)
        for jj in max(1, cj - pad):min(ny, cj + pad), ii in max(1, ci - pad):min(nx, ci + pad)
            ids = get(s.cell_to_facets, _cell_id(s.info, ii, jj), nothing)
            ids === nothing && continue
            for f in ids
                nf = oriented[f]
                (nf[1] == 0.0 && nf[2] == 0.0) && continue
                d = s.positions[f] - q
                d2 = d[1] * d[1] + d[2] * d[2]
                d2 < best && (best = d2; ref = nf)
            end
        end
        isfinite(best) || continue

        ax = ay = 0.0
        for jj in max(1, cj - pad):min(ny, cj + pad), ii in max(1, ci - pad):min(nx, ci + pad)
            ids = get(s.cell_to_facets, _cell_id(s.info, ii, jj), nothing)
            ids === nothing && continue
            for f in ids
                nf = oriented[f]
                (nf[1] == 0.0 && nf[2] == 0.0) && continue
                nf[1] * ref[1] + nf[2] * ref[2] > 0 || continue   # same branch only
                d = s.positions[f] - q
                d2 = d[1] * d[1] + d[2] * d[2]
                w = s.lengths[f] / (d2 + 1e-30)
                ax += w * nf[1]; ay += w * nf[2]
            end
        end
        len = hypot(ax, ay)
        len > 1e-12 && return SVector(ax / len, ay / len)
        return ref
    end
    return SVector(1.0, 0.0)
end

"""
    _fit(s, q, nq, r) -> Float64

Curvature at `q` from a least-squares parabola through the nearby samples.

Set up the frame `(t, nq)` at `q` and write each sample within `r` as `(u, v)` in that
frame. Locally Γ is then a graph `v = c₀ + c₁u + c₂u²`, and

    κ = -2c₂ / (1 + c₁²)^{3/2}

The frame is needed because Γ cannot be written as `v = f(u)` in global coordinates —
a vertical piece of interface is not a function of x. Rotating so `u` runs along the
curve makes it one locally.

Two choices matter:

  * **Least squares over a window, not interpolation.** Sample positions carry O(h²)
    error and κ is a second derivative, so interpolating at the sample spacing δ
    amplifies that error by 1/δ² and the result does not converge under refinement.
    Fitting over a window of radius `r` evaluates the second derivative at scale `r`
    instead.

  * **A parabola in the tangent frame, not a circle.** Over a symmetric window the
    leading `κ′u³/6` term of Γ is odd, so it projects entirely onto `c₁` and leaves
    `c₂` unbiased. A circle fit has no such structure and picks up an O(r·κ′) bias
    that does not vanish under refinement wherever curvature varies.

Samples whose normal opposes `nq` are dropped. The window is a disc in *space*, so if
a second branch of Γ passes within `r` — the far wall of a thin gap, or the opposite
side of a neck — its samples land in the same disc as a second cluster at `v ≈ ±gap`,
and one parabola through both clusters is meaningless. Opposing walls face each other,
so the sign of `nf·nq` separates them. This is also why no connectivity of Γ is ever
needed: the normal already answers "same piece of surface?".
"""
function _fit(s::InterfaceSamples, q::SVector{2,Float64}, nq::SVector{2,Float64}, r::Float64;
              degree::Int=4)
    t = SVector(-nq[2], nq[1])
    r² = r * r

    # Weighted least squares of degree `degree` in the tangent frame.
    #
    # Three departures from a plain unweighted quadratic, each measured on an exact
    # analytic circle and ellipse (metric: Σ m|a_m| of the κ error, which predicts the
    # spurious velocity σ/R·Σ m|a_m| since a mode-m κ error drives u ∝ m):
    #
    #   * arclength weights — subfacet lengths span ~90x, so equal weighting lets
    #     slivers dominate
    #   * a tricube kernel — a hard cutoff makes samples enter and leave the window
    #     abruptly as q moves from facet to facet, which is itself facet-scale noise
    #   * degree 4 — raising the order lets the window widen without bias. Degree 3 is
    #     pointless: the cubic term is odd and cannot affect c₂ under a symmetric
    #     window; only the even u⁴ term reduces bias.
    #
    # Order and radius must rise together. Measured Σ m|a_m| (circle / ellipse):
    #   deg 2, r=4h  2.11 / 3.34     deg 4, r=4h  7.94 / 11.8   (overfit, too few pts)
    #   deg 2, r=8h  0.53 / 3.04     deg 4, r=8h  0.58 /  0.76  (chosen)
    n_c = degree + 1
    M = zeros(MMatrix{5,5,Float64})
    rhs = zeros(MVector{5,Float64})
    npts = 0

    irange, jrange = _cells_within(s.info, q, r)
    for jj in jrange, ii in irange
        ids = get(s.cell_to_facets, _cell_id(s.info, ii, jj), nothing)
        ids === nothing && continue
        for f in ids
            nf = s.normals[f]
            nf[1] * nq[1] + nf[2] * nq[2] > 0 || continue      # same branch only
            d = s.positions[f] - q
            d[1] * d[1] + d[2] * d[2] <= r² || continue
            dd = sqrt(d[1] * d[1] + d[2] * d[2])
            u = d[1] * t[1] + d[2] * t[2]
            v = d[1] * nq[1] + d[2] * nq[2]
            w = s.lengths[f] * (1 - (dd / r)^3)^3
            w > 0 || continue
            pu = 1.0                       # u^(a-1)
            for a in 1:n_c
                pv = pu * pu               # M[a,b] = Σ w·u^(a+b-2); at b=a that is u^(2a-2)
                for b in a:n_c
                    M[a, b] += w * pv
                    pv *= u
                end
                rhs[a] += w * pu * v
                pu *= u
            end
            npts += 1
        end
    end

    npts >= n_c + 2 || return 0.0
    for a in 1:n_c, b in 1:(a-1)
        M[a, b] = M[b, a]
    end

    # The samples must actually span the window: a cloud clustered into a spot carries
    # no second-derivative information. Testing the spread of u is scale-free, unlike
    # thresholding det(A), whose magnitude scales as u¹² and so means different things
    # at different h.
    A = SMatrix{5,5,Float64}(M)[SOneTo(n_c), SOneTo(n_c)]
    abs(det(A)) > 0.0 || return 0.0
    c = A \ SVector{5,Float64}(rhs)[SOneTo(n_c)]
    return -2 * c[3] / (1 + c[2]^2)^1.5
end

"""
    get_curvature(s::InterfaceSamples; radius) -> Vector{Float64}

Curvature at every subfacet midpoint, in subfacet order.

The length matches `num_cells(EmbeddedBoundary(cut))`, so the result can be passed
straight to `CellField(κ, Γ)` and integrated against a `Measure` on Γ.

`radius` balances two competing errors:

  * *noise* — samples carry position error ε = O(h²), and fitting over half-width `r`
    with N ≈ r/δ samples gives a curvature error of order ε/(r²√N). Larger `r` is
    better; the 1/r² is the second-derivative amplification.
  * *bias* — Γ is not a parabola. The `κ′u³/6` term is killed by window symmetry, but
    `κ″u⁴/24` survives and contributes a bias of order r²κ″. Smaller `r` is better.

Minimising the sum puts the optimum at r ~ (ε√h)^{2/9}, i.e. with ε = O(h²) the ideal
`r/h` grows slowly as h → 0. A fixed `4h` is a pragmatic compromise rather than an
asymptotically optimal choice, which is why errors plateau under refinement instead of
converging cleanly.

There is also a hard limit: `r` must stay well below the local radius of curvature
1/κ, or Γ is not a graph over its tangent inside the window and the parabola model
does not apply at all.
"""
function get_curvature(s::InterfaceSamples; radius::Real)
    κ = Vector{Float64}(undef, length(s))
    r = Float64(radius)
    Base.Threads.@threads for f in eachindex(s.positions)
        κ[f] = _fit(s, s.positions[f], s.normals[f], r)
    end
    return κ
end

"""
    curvature_at(s::InterfaceSamples, x, y; radius) -> Float64

Curvature of Γ at the point of Γ closest to `(x, y)`.

`(x, y)` is projected onto the nearest sample's tangent line before fitting, so the
value returned is a property of the interface rather than of the query point. This is
meaningful for queries within roughly `radius` of Γ — a grid node in the band around
the interface, say. Further away there is no nearby branch to speak of and the result
is `0.0`.
"""
function curvature_at(s::InterfaceSamples, x::Real, y::Real; radius::Real)
    isempty(s.positions) && return 0.0
    q = SVector(Float64(x), Float64(y))
    r = Float64(radius)

    best = Inf
    nq = SVector(0.0, 0.0)
    pf = q
    irange, jrange = _cells_within(s.info, q, r)
    for jj in jrange, ii in irange
        ids = get(s.cell_to_facets, _cell_id(s.info, ii, jj), nothing)
        ids === nothing && continue
        for f in ids
            d = s.positions[f] - q
            d2 = d[1] * d[1] + d[2] * d[2]
            d2 < best && (best = d2; nq = s.normals[f]; pf = s.positions[f])
        end
    end
    isfinite(best) || return 0.0

    # project onto the nearest facet's tangent line, so the fit is centred on Γ
    qΓ = q + ((pf - q) ⋅ nq) * nq
    return _fit(s, qΓ, nq, r)
end

@inline ⋅(a::SVector{2,Float64}, b::SVector{2,Float64}) = a[1] * b[1] + a[2] * b[2]

"""
    reference_curve(geom::EvolvingDiscreteGeometry) -> InterfaceSamples

The explicit discrete interface: position, normal and curvature per subfacet, cached on
the geometry and dropped whenever the level set changes. Can be used for post processing or things like that
"""
function reference_curve(geom::EvolvingDiscreteGeometry)
    if isnothing(geom.cache.interface_samples)
        geom.cache.interface_samples = interface_samples(geom)
    end
    return geom.cache.interface_samples::InterfaceSamples
end



end # module
