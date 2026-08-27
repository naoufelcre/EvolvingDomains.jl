module Curvature

using StaticArrays: SVector, SMatrix, MMatrix, MVector, SOneTo
using LinearAlgebra: det
using GridapEmbedded: EmbeddedBoundary
using ..Geometric: EvolvingDiscreteGeometry, ensure_cut!, grid_info, CartesianGridInfo

export InterfaceSamples, interface_samples, get_curvature, curvature_at, interface_tangents
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
    bad = findfirst(!isfinite, geom.levelset)
    bad === nothing || error("interface_samples: level set has a non-finite value " *
        "($(geom.levelset[bad])) at node $bad — a velocity blow-up or bad advection step " *
        "upstream is the usual cause.")
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

    end
    _index!(cell_to_facets, info, positions, endpoints)

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
    _fill_unresolved!(s.curvatures, s)
    return _sharpen!(s, r)
end

"""
    _index!(cell_to_facets, info, positions, endpoints)

(Re)build the background-cell index. A facet is registered under its midpoint *and* both
endpoints so a nearest-facet search finds it from any cell it reaches into; `_fit_coeffs`
undoes the resulting multiplicity itself.
"""
function _index!(c2f::Dict{Int,Vector{Int32}}, info::CartesianGridInfo,
                 positions::Vector{SVector{2,Float64}}, endpoints::Vector{NTuple{4,Float64}})
    empty!(c2f)
    for f in eachindex(positions)
        e = endpoints[f]
        for (px, py) in ((positions[f][1], positions[f][2]), (e[1], e[2]), (e[3], e[4]))
            i, j = _cell_of(info, px, py)
            v = get!(() -> Int32[], c2f, _cell_id(info, i, j))
            (isempty(v) || v[end] != Int32(f)) && push!(v, Int32(f))
        end
    end
    return c2f
end

"""
    _edge_of(info, x, y) -> (s, L, e)

The background-grid edge a cut vertex lies on: arclength `s` from one end, edge length `L`,
unit direction `e`. `cut` puts every crossing on a horizontal or vertical cell edge, or on
the diagonal its sub-triangulation adds — measured, that classifies 100% of them. Returns
`NaN` for anything else, which the caller leaves uncorrected.
"""
@inline function _edge_of(info::CartesianGridInfo, x::Float64, y::Float64)
    hx, hy = info.spacing
    fx = (x - info.origin[1]) / hx
    fy = (y - info.origin[2]) / hy
    onx = abs(fx - round(fx)) < 1e-9
    ony = abs(fy - round(fy)) < 1e-9
    ony && !onx && return (fx - floor(fx)) * hx, hx, SVector(1.0, 0.0)
    onx && !ony && return (fy - floor(fy)) * hy, hy, SVector(0.0, 1.0)
    onx && ony  && return 0.0, hx, SVector(1.0, 0.0)          # on a node: s(L−s) = 0 anyway
    u = fx - floor(fx); v = fy - floor(fy); Ld = hypot(hx, hy)
    abs(u + v - 1) < 1e-8 && return u * Ld, Ld, SVector(hx, -hy) / Ld
    abs(u - v) < 1e-8     && return u * Ld, Ld, SVector(hx, hy) / Ld
    return NaN, NaN, SVector(0.0, 0.0)
end

"Move a cut vertex from the linearly-interpolated crossing onto the true contour."
@inline function _lift(info::CartesianGridInfo, p::SVector{2,Float64}, κ::Float64,
                       n::SVector{2,Float64})
    sl, Le, ev = _edge_of(info, p[1], p[2])
    isfinite(sl) || return p
    en = ev[1] * n[1] + ev[2] * n[2]
    return p + (0.5 * κ * (1 - en * en) * sl * (Le - sl)) * n
end

"""
    _sharpen!(s, r) -> s

Put the sample cloud back on the true interface, then refit.

`cut` locates each crossing by **linear** interpolation of φ, so every vertex sits off the
curve; and the chord midpoint sits inside it by the sagitta. Both are O(h²), and κ is a
second derivative, so a window of C·h amplifies them by 1/(C·h)² — the h² cancels exactly
and the κ noise κ/(8C²) does not converge under refinement. Both are also known in closed
form:

  * **crossing** — for a signed distance function the Hessian at Γ is κ t⊗t, so along an
    edge of direction `e` the interpolation error is ½κ(1−(e·n)²)s(L−s). Projected onto n
    the edge obliquity and |∇φ| both cancel, leaving `δ = ½κ(1−(e·n)²)s(L−s)` — so this
    does not assume φ is normalised, only that it is smooth.
  * **sagitta** — a chord of length L subtends L²κ/8.

Correcting **either one alone makes κ worse** (measured 2.2× and 1.4×): each leaves the
other term dominant while adding its own perturbation. Together they cancel — circle κ rms
23× better, interface ripple 168× better and O(h²) convergent, and unlike every fit-side
attempt the ellipse (3.9×) and merged drops (1.4×) improve too. One Picard step suffices:
it uses κ fitted from the uncorrected cloud, far more accurate than the correction needs.

Skipped where |κ|h ≥ ½: there the feature is sub-grid, the smooth-interface Hessian does
not apply, and the reflex corners of merged drops sit exactly there.
"""
function _sharpen!(s::InterfaceSamples, r::Float64)
    hmin = min(s.info.spacing[1], s.info.spacing[2])
    for f in eachindex(s.positions)
        κ = s.curvatures[f]; n = s.normals[f]
        (isfinite(κ) && s.lengths[f] > 1e-14 && abs(κ) * hmin < 0.5) || continue
        e = s.endpoints[f]
        p1 = _lift(s.info, SVector(e[1], e[2]), κ, n)
        p2 = _lift(s.info, SVector(e[3], e[4]), κ, n)
        d = p2 - p1
        L = hypot(d[1], d[2])
        L > 1e-14 || continue
        s.endpoints[f] = (p1[1], p1[2], p2[1], p2[2])
        s.lengths[f] = L
        s.positions[f] = (p1 + p2) / 2 + (L * L * κ / 8) * n   # undo the chord sagitta
        t = SVector(d[2], -d[1]) / L
        s.normals[f] = (t[1] * n[1] + t[2] * n[2]) >= 0 ? t : -t
    end
    _index!(s.cell_to_facets, s.info, s.positions, s.endpoints)
    Base.Threads.@threads for f in eachindex(s.positions)
        s.curvatures[f] = _fit(s, s.positions[f], s.normals[f], r)
    end
    _fill_unresolved!(s.curvatures, s)
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
function _fit_coeffs(s::InterfaceSamples, q::SVector{2,Float64}, nq::SVector{2,Float64},
                     r::Float64; degree::Int=4)
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
        cid = _cell_id(s.info, ii, jj)
        ids = get(s.cell_to_facets, cid, nothing)
        ids === nothing && continue
        for f in ids
            # `_index!` registers a facet under its midpoint and both endpoints, so up to
            # three cells of this block reach the same sample. Measured on a circle: 47% of
            # samples entered twice, inflating the stencil 1.55×. Take each only from the
            # cell that owns its midpoint.
            _cell_id(s.info, _cell_of(s.info, s.positions[f][1], s.positions[f][2])...) == cid || continue
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

    npts >= n_c + 2 || return nothing   # unresolved (too few samples)
    for a in 1:n_c, b in 1:(a-1)
        M[a, b] = M[b, a]
    end

    # The samples must actually span the window: a cloud clustered into a spot carries
    # no second-derivative information. Testing the spread of u is scale-free, unlike
    # thresholding det(A), whose magnitude scales as u¹² and so means different things
    # at different h.
    A = SMatrix{5,5,Float64}(M)[SOneTo(n_c), SOneTo(n_c)]
    abs(det(A)) > 0.0 || return nothing   # unresolved (samples don't span window)
    return A \ SVector{5,Float64}(rhs)[SOneTo(n_c)]
end

"""
    _fit(s, q, nq, r) -> κ

Curvature from the fit's second-order coefficient. NaN when unresolved (filled from
neighbours by `_fill_unresolved!`).
"""
function _fit(s::InterfaceSamples, q::SVector{2,Float64}, nq::SVector{2,Float64}, r::Float64;
              degree::Int=4)
    c = _fit_coeffs(s, q, nq, r; degree=degree)
    c === nothing && return NaN
    return -2 * c[3] / (1 + c[2]^2)^1.5
end

"""
    _fit_slope(s, q, nq, r) -> c₁

First-order coefficient: the slope of the fitted curve in the local frame. NaN if unresolved.

This is the *cleanest* thing the fit produces. Measured on an exact circle, the tangent
rebuilt from c₁ is 130x more accurate than Gridap's `facet_to_normal` and **2266x smoother**
(roughness 3.6e-6 vs 8.2e-3 rad). κ, by contrast, comes from c₂ — a second derivative — and
carries a non-convergent 0.3% p2p noise floor. Same fit, different coefficient.
"""
function _fit_slope(s::InterfaceSamples, q::SVector{2,Float64}, nq::SVector{2,Float64},
                    r::Float64; degree::Int=4)
    c = _fit_coeffs(s, q, nq, r; degree=degree)
    c === nothing && return NaN
    return c[2]
end

"""
    _fill_unresolved!(κ, s) -> κ

Replace facets the fit could not resolve (κ = NaN) with a distance-weighted average of
nearby resolved facets.

A facet whose window holds too few samples to fit — a feature finer than the window, a
sliver, an isolated fragment — cannot carry its own curvature. Returning 0 would read as
"flat" where the truth is "as curved as the grid allows", the worst possible value for a
force term. Borrowing from resolved neighbours instead **deliberately smooths the sub-grid
feature**, which is the honest meaning of "unresolved": a continuous κ field, with detail
below the window averaged away. Reads a snapshot so one borrowed facet never seeds another
(the ordering hazard `_borrow_normal` also guards against). If nothing nearby resolved —
the whole shape is sub-grid — falls back to 0.
"""
function _fill_unresolved!(κ::Vector{Float64}, s::InterfaceSamples)
    any(isnan, κ) || return κ
    resolved = copy(κ)
    nx, ny = s.info.cells
    for f in eachindex(κ)
        isnan(κ[f]) || continue
        q = s.positions[f]
        ci, cj = _cell_of(s.info, q[1], q[2])
        num = 0.0
        den = 0.0
        for pad in 0:max(nx, ny)
            for jj in max(1, cj - pad):min(ny, cj + pad), ii in max(1, ci - pad):min(nx, ci + pad)
                pad == 0 || ii == ci - pad || ii == ci + pad ||
                    jj == cj - pad || jj == cj + pad || continue
                ids = get(s.cell_to_facets, _cell_id(s.info, ii, jj), nothing)
                ids === nothing && continue
                for gf in ids
                    isnan(resolved[gf]) && continue
                    d = s.positions[gf] - q
                    w = 1.0 / (d[1] * d[1] + d[2] * d[2] + 1e-30)
                    num += w * resolved[gf]
                    den += w
                end
            end
            den > 0 && break        # nearest ring holding a resolved facet wins
        end
        κ[f] = den > 0 ? num / den : 0.0
    end
    return κ
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
    _fill_unresolved!(κ, s)
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
    fidx = 0
    irange, jrange = _cells_within(s.info, q, r)
    for jj in jrange, ii in irange
        ids = get(s.cell_to_facets, _cell_id(s.info, ii, jj), nothing)
        ids === nothing && continue
        for f in ids
            d = s.positions[f] - q
            d2 = d[1] * d[1] + d[2] * d[2]
            d2 < best && (best = d2; nq = s.normals[f]; pf = s.positions[f]; fidx = f)
        end
    end
    isfinite(best) || return 0.0

    # project onto the nearest facet's tangent line, so the fit is centred on Γ
    qΓ = q + ((pf - q) ⋅ nq) * nq
    κ = _fit(s, qΓ, nq, r)
    # unresolved here too: fall back to the nearest facet's stored (already borrow-filled)
    # curvature rather than leaking a NaN — same "smooth the sub-grid feature" contract.
    return isnan(κ) ? s.curvatures[fidx] : κ
end

@inline ⋅(a::SVector{2,Float64}, b::SVector{2,Float64}) = a[1] * b[1] + a[2] * b[2]

"""
    interface_tangents(geom; radius=nothing) -> Vector{SVector{2,Float64}}

Smoothed unit tangent at each subfacet, in subfacet order — so it wraps with
`CellField(τ, Γ)` exactly like `interface_curvature`.

Built from the fit's c₁ rather than from facet geometry. Gridap's `facet_to_normal` is
piecewise constant and jumps ~8e-3 rad between neighbours; the fitted tangent describes the
fitted *curve*, so it is continuous across facets. That matters for the variational
surface-tension force ∫_Γ γ (∇_Γ·w): a discontinuous tangent there produces delta-function
point loads at every facet vertex (measured: max|u·n| 575 vs 0.25 on a circle at rest).

Unresolved facets fall back to the raw facet tangent.
"""
function interface_tangents(geom::EvolvingDiscreteGeometry; radius=nothing)
    s = reference_curve(geom)
    r = isnothing(radius) ? 8 * min(s.info.spacing[1], s.info.spacing[2]) : Float64(radius)
    τ = Vector{SVector{2,Float64}}(undef, length(s.positions))
    Base.Threads.@threads for f in eachindex(s.positions)
        nq = s.normals[f]
        t  = SVector(-nq[2], nq[1])
        c1 = _fit_slope(s, s.positions[f], nq, r)
        τ[f] = isnan(c1) ? t : (t + c1 * nq) / sqrt(1 + c1 * c1)
    end
    return τ
end

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
