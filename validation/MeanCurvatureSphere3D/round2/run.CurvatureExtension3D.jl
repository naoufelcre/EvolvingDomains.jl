#!/usr/bin/env julia
"""
Historic round. Superseded by `round3/run.CurvatureExtension3D.jl`, which replaces this
file's local facet-to-node adapter with the package operator `extend!` and is the
maintained form of this experiment; the outputs here are kept for comparison.

CurvatureExtension3D — mean-curvature flow driven by Curvature.jl.

    EvolvingDomains CutSurface  ->  Curvature.jl (one kappa per facet)
                                ->  facet normal speed  V_f = -mu * kappa_f
                                ->  closest-point extension of V_f * n_f to background nodes
                                ->  EvolvingDomains eulerianTransport! (WENO5 + SSP-RK3)

    julia --project=../env run.CurvatureExtension3D.jl            # 32^3, both geometries
    julia --project=../env run.CurvatureExtension3D.jl 40 0.12    # n, tmax

Two independent level sets are advanced with the same grid, same dt sequence, same transport
and same reinitialisation schedule:

  * `cv`  — velocity from Curvature.jl facet curvature, extended to the grid;
  * `fd`  — velocity from a nodal finite-difference curvature (div of nodal normals), the
            comparison method of round 1, kept unchanged.

Curvature.jl drives the `cv` interface at every step; it is not a diagnostic there.

Sign convention: phi < 0 inside, n outward, kappa = div n > 0 on a convex surface, so
V_n = -mu*kappa < 0 makes the sphere contract. Geometry B (axisymmetric P2 perturbation) is
driven only by its numerically computed curvature; no analytic trajectory is imposed.
"""

using Libdl: dlext
using Printf
using LinearAlgebra: norm, dot, cross
using StaticArrays: SVector

# ---------------------------------------------------------------------------------------
# Environment: native cut library before `using EvolvingDomains`; sibling Curvature.jl
# checkout on LOAD_PATH (no project, manifest or package source is modified).
# ---------------------------------------------------------------------------------------
function locateCutLibrary()
    override = get(ENV, "EVOLVINGDOMAINS_CUTCELLS_LIBRARY", "")
    if !isempty(override)
        isfile(override) || error("EVOLVINGDOMAINS_CUTCELLS_LIBRARY=$override does not exist.")
        return override
    end
    name = "libevolvingdomains_cut_cells." * dlext
    for dir in ("evolvingdomains-cutcells-install", "evolvingdomains-cutcells-install-mc")
        candidate = joinpath(tempdir(), dir, "lib", name)
        isfile(candidate) && return candidate
    end
    error("no native cut library: build one with scripts/build_cxx_bindings.sh and set " *
          "EVOLVINGDOMAINS_CUTCELLS_LIBRARY.")
end

function locateCurvatureSource()
    override = get(ENV, "CURVATURE_JL_DIR", "")
    isempty(override) || return override
    dir = @__DIR__
    for _ in 1:8
        candidate = joinpath(dirname(dir), "Curvature.jl", "src", "Curvature.jl")
        isfile(candidate) && return dirname(dirname(candidate))
        dir = dirname(dir)
    end
    error("Curvature.jl source checkout not found; set CURVATURE_JL_DIR.")
end

ENV["EVOLVINGDOMAINS_CUTCELLS_LIBRARY"] = locateCutLibrary()
push!(LOAD_PATH, locateCurvatureSource())

using EvolvingDomains
using Curvature: interfaceCurvature

# ---------------------------------------------------------------------------------------
# Configuration
# ---------------------------------------------------------------------------------------
Base.@kwdef struct Config
    n::Int = 32
    R0::Float64 = 0.23
    mu::Float64 = 0.01
    tmax::Float64 = 0.12
    epsilon::Float64 = 0.15                  # geometry B perturbation amplitude
    # Explicit-diffusion scale: curvature flow with speed -mu*kappa is stable for
    # mu*dt/h^2 <= 1/6 as a diffusion operator; C = 1/12 is a factor-two margin.
    # Not a rigorous bound for the full WENO5 + extension + reinit algorithm.
    dt_coefficient::Float64 = 1 / 12
    reinit_every::Int = 10
    band_in::Float64 = 1.5                   # full speed for |phi| <= band_in*h
    band_out::Float64 = 3.0                  # velocity zero for |phi| >= band_out*h
    search_radius::Int = 2                   # cell block for the nearest-facet search
    eps_gradient::Float64 = 0.5              # |grad phi| clamp for the FD method
    curvature_radius::Float64 = 4.0          # Curvature.jl fit window, in cells
    snapshot_times::Vector{Float64} = Float64[]
end

exactRadius(cfg, t) = sqrt(max(cfg.R0^2 - 4 * cfg.mu * t, 0.0))
exactVolume(cfg, t) = 4π * exactRadius(cfg, t)^3 / 3

# ---------------------------------------------------------------------------------------
# Initial level sets
# ---------------------------------------------------------------------------------------

"""
    sphereField(x, R0) — geometry A: phi0 = |x| - R0, an exact signed-distance function.

    perturbedField(x, R0, eps) — geometry B: phi0 = r - R0*(1 + eps*P2(z/r)).

The perturbed level set is **not** a signed-distance function: |grad phi0| differs from 1 by
O(eps) and depends on direction, so the curvature of the *field* and the curvature of the
*interface* are not interchangeable there either. It is a valid level set (correct sign, single
closed zero crossing), so the cut geometry is right; only its metric is not Euclidean. The run
therefore reinitialises it, measuring the cut geometry before and after.

At r = 0 the direction z/r is undefined. Every directional limit of phi0 at the centre is
negative (between -R0*(1-eps/2) and -R0*(1+eps)), so the sign, and hence the interface, is
unaffected; the value stored there is the angular average -R0*(1-eps/2), which keeps the
interior nodal field finite and scalar without dividing by zero.
"""
function sphereField(x, R0)
    return [norm(SVector(a, b, c)) - R0 for a in x, b in x, c in x]
end

function perturbedField(x, R0, eps)
    phi = Array{Float64}(undef, length(x), length(x), length(x))
    for k in eachindex(x), j in eachindex(x), i in eachindex(x)
        px, py, pz = x[i], x[j], x[k]
        r = sqrt(px * px + py * py + pz * pz)
        s = r > 1e-12 ? pz / r : 0.0
        phi[i, j, k] = r - R0 * (1 + eps * (3s^2 - 1) / 2)
    end
    return phi
end

# ---------------------------------------------------------------------------------------
# Interface extraction from the public CutSurface
# ---------------------------------------------------------------------------------------

"""
    interfaceData(field, x) -> NamedTuple

Public `Cut` -> public `CutSurface` -> plain geometry. Per facet `q`: vertices from the CSR
row `point_pointers[q]:point_pointers[q+1]-1` of `point_ids` into `physical_points`, and
`boundary.normals[q]`, already oriented out of the phi < 0 phase. `cellof[q]` is the linear
background-cell index of the facet's parent cell (first axis fastest, cell counts dims-1),
used to index facets spatially without an all-pairs search.

Also returns area-weighted radius, sphericity (relative area-weighted std of the centroid
radius) and the divergence-theorem volume of the triangulated zero set.
"""
function interfaceData(field, x)
    surface = cut(field, :init).boundary
    nf = length(surface.facet_to_bgcell)
    facets = Vector{NTuple{3,SVector{3,Float64}}}(undef, nf)
    normals = Vector{SVector{3,Float64}}(undef, nf)
    centroids = Vector{SVector{3,Float64}}(undef, nf)
    areas = Vector{Float64}(undef, nf)
    @inbounds for q in 1:nf
        a1 = Int(surface.point_pointers[q])
        a2 = Int(surface.point_pointers[q + 1]) - 1
        a2 - a1 + 1 == 3 || error("facet $q has $(a2 - a1 + 1) vertices, expected 3")
        p1 = SVector{3,Float64}(surface.physical_points[Int(surface.point_ids[a1])])
        p2 = SVector{3,Float64}(surface.physical_points[Int(surface.point_ids[a1 + 1])])
        p3 = SVector{3,Float64}(surface.physical_points[Int(surface.point_ids[a1 + 2])])
        facets[q] = (p1, p2, p3)
        normals[q] = SVector{3,Float64}(surface.normals[q])
        centroids[q] = (p1 + p2 + p3) / 3
        areas[q] = norm(cross(p2 - p1, p3 - p1)) / 2
    end
    total = sum(areas)
    radii = norm.(centroids)
    radius = total > 0 ? sum(areas .* radii) / total : NaN
    spread = total > 0 ? sqrt(max(sum(areas .* (radii .- radius) .^ 2) / total, 0.0)) : NaN
    volume = sum(areas[q] * dot(centroids[q], normals[q]) for q in 1:nf) / 3
    return (facets=facets, normals=normals, centroids=centroids, areas=areas,
            cellof=surface.facet_to_bgcell, radius=radius, radius_std=spread,
            sphericity=radius > 0 ? spread / radius : NaN,
            volume=volume, area=total, nfacets=nf)
end

# ---------------------------------------------------------------------------------------
# Curvature.jl: one mean curvature per interface facet
# ---------------------------------------------------------------------------------------

"""
    facetCurvature(geo, x, cfg) -> Vector{Float64}

The required cross-package call: plain facet geometry from EvolvingDomains in, one curvature
value per facet out, in facet order (NaN where Curvature.jl cannot resolve the facet). The
lattice (`spacing`, `origin`) supplies the fit window and the neighbour index; it is the same
grid the cut was taken from.
"""
function facetCurvature(geo, x, cfg)
    h = x[2] - x[1]
    return interfaceCurvature(geo.facets, geo.normals; spacing=h, origin=x[1],
                              radius=cfg.curvature_radius * h)
end

# ---------------------------------------------------------------------------------------
# Facet -> node extension (application-local; not added to either package)
# ---------------------------------------------------------------------------------------

"""
    FacetIndex(cellof, cells) — CSR of facets grouped by background cell.

`cellof[q]` is the parent cell of facet `q`; the arrays let a node look up the facets of a
small block of cells instead of scanning all facets.
"""
struct FacetIndex
    start::Vector{Int32}
    facets::Vector{Int32}
end

function FacetIndex(cellof::Vector{Int32}, ncells::Int)
    start = zeros(Int32, ncells + 1)
    for c in cellof
        1 <= c <= ncells || error("facet parent cell $c outside 1:$ncells")
        start[c + 1] += 1
    end
    for c in 1:ncells
        start[c + 1] += start[c]
    end
    facets = Vector{Int32}(undef, length(cellof))
    fill = copy(start)
    for q in eachindex(cellof)
        c = cellof[q]
        facets[fill[c] + 1] = Int32(q)
        fill[c] += 1
    end
    return FacetIndex(start, facets)
end

"""
    closestPointOnTriangle(a, b, c, p) -> SVector{3,Float64}

Exact closest point of `p` on the closed triangle (Ericson, *Real-Time Collision Detection*,
§5.1.5): seven Voronoi regions, no iteration, no centroid approximation.
"""
function closestPointOnTriangle(a, b, c, p)
    ab, ac, ap = b - a, c - a, p - a
    d1, d2 = dot(ab, ap), dot(ac, ap)
    if d1 <= 0 && d2 <= 0
        return a
    end
    bp = p - b
    d3, d4 = dot(ab, bp), dot(ac, bp)
    if d3 >= 0 && d4 <= d3
        return b
    end
    vc = d1 * d4 - d3 * d2
    if vc <= 0 && d1 >= 0 && d3 <= 0
        return a + (d1 / (d1 - d3)) * ab
    end
    cp = p - c
    d5, d6 = dot(ab, cp), dot(ac, cp)
    if d6 >= 0 && d5 <= d6
        return c
    end
    vb = d5 * d2 - d1 * d6
    if vb <= 0 && d2 >= 0 && d6 <= 0
        return a + (d2 / (d2 - d6)) * ac
    end
    va = d3 * d6 - d5 * d4
    if va <= 0 && (d4 - d3) >= 0 && (d5 - d6) >= 0
        return b + ((d4 - d3) / ((d4 - d3) + (d5 - d6))) * (c - b)
    end
    denom = 1 / (va + vb + vc)
    return a + ab * (vb * denom) + ac * (vc * denom)
end

"""
    assignNearestFacets!(facet_of, phi, x, geo, index, cfg) -> (band, unresolved)

Closest-point extension index. For every node of the narrow band |phi| <= band_out*h, store
the facet minimizing the true point-to-triangle distance (`0` outside the band, `-1` inside
the band with no facet found in the searched block). Candidate facets come from the CSR block
of background cells within `search_radius` cells of the node; the search widens once if the
block is empty, so the cost is proportional to band nodes times facets of a few cells, never
to nodes times facets of the whole interface.

This is a **closest-point extension of interface values**, not a normal (PDE) extension: it
gives v_f(x) = value of the facet nearest x, which is exactly the map the facet data supports.
It is constant along the direction normal to that facet (piecewise, facet by facet), which is
why the resulting velocity field is only O(h)-smooth across facet boundaries.
"""
function assignNearestFacets!(facet_of, phi, x, geo, index::FacetIndex, cfg)
    n = length(x)
    h = x[2] - x[1]
    nc = n - 1
    rmax = cfg.band_out * h
    band = 0
    unresolved = 0
    @inbounds for k in 1:n, j in 1:n, i in 1:n
        d = abs(phi[i, j, k])
        if d > rmax
            facet_of[i, j, k] = 0
            continue
        end
        band += 1
        p = SVector{3,Float64}(x[i], x[j], x[k])
        best = Inf
        bf = 0
        for radius in (cfg.search_radius, 6)
            i1, i2 = max(1, i - 1 - radius), min(nc, i + radius)
            j1, j2 = max(1, j - 1 - radius), min(nc, j + radius)
            k1, k2 = max(1, k - 1 - radius), min(nc, k + radius)
            for ck in k1:k2, cj in j1:j2, ci in i1:i2
                cell = ci + (cj - 1) * nc + (ck - 1) * nc * nc
                for t in (index.start[cell] + 1):index.start[cell + 1]
                    q = Int(index.facets[t])
                    f = geo.facets[q]
                    cp = closestPointOnTriangle(f[1], f[2], f[3], p)
                    dd = (cp - p)[1]^2 + (cp - p)[2]^2 + (cp - p)[3]^2
                    if dd < best
                        best = dd
                        bf = q
                    end
                end
            end
            bf != 0 && break
        end
        facet_of[i, j, k] = Int32(bf)
        bf == 0 && (unresolved += 1)
    end
    return (band=band, unresolved=unresolved)
end

"""
    facetVelocity!(velocity, phi, x, geo, kappa, facet_of, cfg) -> (nodes, nonfinite)

Background velocity from the per-facet curvature:

    v_h(x) = chi(|phi|) * (-mu * kappa_{f(x)}) * n_{f(x)}

i.e. the facet's own outward normal and the normal speed V_f = -mu*kappa_f, extended to the
node by the nearest facet f(x). `chi` is a C^1 cosine taper: 1 for |phi| <= band_in*h, 0 for
|phi| >= band_out*h. The zero level set sits in the full-speed zone, so on the interface the
transport sees v.nabla(phi) = -mu*kappa_f*(n_f.nabla(phi)) = -mu*kappa|nabla(phi)| up to the
angle between the facet normal and the grid gradient.

A facet with non-finite curvature contributes zero velocity at its nodes (counted, never
hidden): the alternative, carrying NaN into the transport, is rejected by
`eulerianTransport!` for good reason.

Away from the band the velocity is exactly zero. That is documented behaviour, not an
extension: no interface value is propagated to the far field, and the level set there is only
transported passively.

shortcut: velocity is tapered to zero at band_out*h instead of being extended along normals to
the whole box; adequate while the band exceeds the per-step interface displacement (0.02 h
here). Upgrade to a normal extension if the band ever lags the front.
"""
function facetVelocity!(velocity, phi, x, geo, kappa, facet_of, cfg)
    n = length(x)
    h = x[2] - x[1]
    rin, rout = cfg.band_in * h, cfg.band_out * h
    nodes = 0
    nonfinite = 0
    @inbounds for k in 1:n, j in 1:n, i in 1:n
        q = facet_of[i, j, k]
        v = SVector(0.0, 0.0, 0.0)
        if q > 0
            rate = -cfg.mu * kappa[q]
            if isfinite(rate)
                d = abs(phi[i, j, k])
                chi = d <= rin ? 1.0 :
                      d >= rout ? 0.0 :
                      (c = cos(π / 2 * (d - rin) / (rout - rin)); c * c)
                nf = geo.normals[q]
                v = SVector(chi * rate * nf[1], chi * rate * nf[2], chi * rate * nf[3])
                nodes += 1
            else
                nonfinite += 1
            end
        end
        velocity[i + (j - 1) * n + (k - 1) * n * n] = v
    end
    return (nodes=nodes, nonfinite=nonfinite)
end

"""
    facetValuesToNodes(geo, values, facet_of, phi) -> Vector{Float64}

The scalar form of the same extension (used by the Stage C known-surface test): nodal value
= value of the nearest facet, zero outside the band and where the facet holds no value.
"""
function facetValuesToNodes(geo, values, facet_of, phi)
    out = zeros(Float64, size(phi))
    @inbounds for index in CartesianIndices(phi)
        q = facet_of[index]
        q > 0 && (out[index] = values[q])
    end
    return out
end

# ---------------------------------------------------------------------------------------
# Independent finite-difference curvature baseline (round-1 operator, unchanged)
# ---------------------------------------------------------------------------------------

computeNormals!(normals, phi, h, eps) = begin
    n = size(phi, 1)
    inv2h = 1 / (2h)
    clamped = 0
    @inbounds for k in 1:n, j in 1:n, i in 1:n
        gx = (phi[clamp(i + 1, 1, n), j, k] - phi[clamp(i - 1, 1, n), j, k]) * inv2h
        gy = (phi[i, clamp(j + 1, 1, n), k] - phi[i, clamp(j - 1, 1, n), k]) * inv2h
        gz = (phi[i, j, clamp(k + 1, 1, n)] - phi[i, j, clamp(k - 1, 1, n)]) * inv2h
        m = sqrt(gx * gx + gy * gy + gz * gz)
        s = m > eps ? 1 / m : 1 / eps
        m > eps || (clamped += 1)
        normals[i, j, k] = SVector(gx * s, gy * s, gz * s)
    end
    clamped
end

computeCurvature!(kappa, normals, h) = begin
    n = size(kappa, 1)
    inv2h = 1 / (2h)
    @inbounds for k in 1:n, j in 1:n, i in 1:n
        a = normals[clamp(i + 1, 1, n), j, k]
        b = normals[clamp(i - 1, 1, n), j, k]
        c = normals[i, clamp(j + 1, 1, n), k]
        d = normals[i, clamp(j - 1, 1, n), k]
        e = normals[i, j, clamp(k + 1, 1, n)]
        f = normals[i, j, clamp(k - 1, 1, n)]
        kappa[i, j, k] = ((a[1] - b[1]) + (c[2] - d[2]) + (e[3] - f[3])) * inv2h
    end
    count(!isfinite, kappa)
end

"""
    nodalVelocity!(velocity, phi, kappa, normals, h, cfg) -> full

Nodal finite-difference variant: v = chi(|phi|) * (-mu*kappa_h) * n_h with the same taper.
Kept identical in band treatment to `facetVelocity!` so the comparison isolates the curvature
operator, not the velocity assembly.
"""
function nodalVelocity!(velocity, phi, kappa, normals, h, cfg)
    n = size(phi, 1)
    rin, rout = cfg.band_in * h, cfg.band_out * h
    full = 0
    @inbounds for k in 1:n, j in 1:n, i in 1:n
        d = abs(phi[i, j, k])
        chi = d <= rin ? 1.0 : d >= rout ? 0.0 :
              (c = cos(π / 2 * (d - rin) / (rout - rin)); c * c)
        if d <= rin && chi > 0
            full += 1
        end
        rate = -cfg.mu * kappa[i, j, k]
        nv = normals[i, j, k]
        s = isfinite(rate) ? chi * rate : 0.0
        velocity[i + (j - 1) * n + (k - 1) * n * n] = SVector(s * nv[1], s * nv[2], s * nv[3])
    end
    return full
end

# ---------------------------------------------------------------------------------------
# Diagnostics
# ---------------------------------------------------------------------------------------

"""
    gradientBandStats(phi, x, h) -> (mean, min, max)

Magnitude of the central-difference gradient over the narrow band |phi| <= 0.8 h. For a
level set that behaves as a distance function this stays near 1; a growing value means the
band is being compressed by the velocity taper (all band level sets advected at nearly the
same normal speed, the tapered ones lagging), which is the mechanism that degrades WENO5
reconstruction if `reinitialize!` is not keeping it in check.
"""
function gradientBandStats(phi, x, h)
    n = length(x)
    inv2h = 1 / (2h)
    total, count_, lo, hi = 0.0, 0, Inf, -Inf
    @inbounds for k in 1:n, j in 1:n, i in 1:n
        abs(phi[i, j, k]) <= 0.8 * h || continue
        gx = (phi[clamp(i + 1, 1, n), j, k] - phi[clamp(i - 1, 1, n), j, k]) * inv2h
        gy = (phi[i, clamp(j + 1, 1, n), k] - phi[i, clamp(j - 1, 1, n), k]) * inv2h
        gz = (phi[i, j, clamp(k + 1, 1, n)] - phi[i, j, clamp(k - 1, 1, n)]) * inv2h
        g = sqrt(gx * gx + gy * gy + gz * gz)
        total += g
        count_ += 1
        lo = min(lo, g)
        hi = max(hi, g)
    end
    return count_ == 0 ? (NaN, NaN, NaN) : (total / count_, lo, hi)
end

"""
    sampleTrilinear(phi, x, p)

Trilinear interpolation of a nodal field on the uniform cubic grid `x`; used to read the
nodal finite-difference curvature at Curvature.jl's facet centroids.
"""
function sampleTrilinear(phi, x, p)
    n = length(x)
    h = x[2] - x[1]
    s = ntuple(d -> clamp((p[d] - x[1]) / h + 1, 1.0, float(n)), 3)
    lo = ntuple(d -> clamp(floor(Int, s[d]), 1, n - 1), 3)
    f = ntuple(d -> s[d] - lo[d], 3)
    i, j, k = lo
    fx, fy, fz = f
    c00 = phi[i, j, k] * (1 - fx) + phi[i + 1, j, k] * fx
    c10 = phi[i, j + 1, k] * (1 - fx) + phi[i + 1, j + 1, k] * fx
    c01 = phi[i, j, k + 1] * (1 - fx) + phi[i + 1, j, k + 1] * fx
    c11 = phi[i, j + 1, k + 1] * (1 - fx) + phi[i + 1, j + 1, k + 1] * fx
    return (c00 * (1 - fy) + c10 * fy) * (1 - fz) + (c01 * (1 - fy) + c11 * fy) * fz
end

"Distance from the centre to the zero crossing of the trilinear interpolant along a coordinate axis."
function axisCrossing(phi, x, unit::Int, side::Int)
    at = s -> sampleTrilinear(phi, x, ntuple(d -> d == unit ? side * s : 0.0, 3))
    at(0.0) < 0 || return NaN
    at(0.5) > 0 || return NaN
    lo, hi = 0.0, 0.5
    for _ in 1:60
        mid = (lo + hi) / 2
        at(mid) < 0 ? (lo = mid) : (hi = mid)
    end
    return (lo + hi) / 2
end

"""
    curvatureSummary(geo, kappa_cv, kappa_fd)

Area-weighted mean/std of the Curvature.jl facet curvatures, the band mean of the nodal
finite-difference curvature, and the area-weighted rms mismatch between the two at the facet
centroids (relative to |kappa| = 2/R0, the natural scale here).
"""
function curvatureSummary(geo, kappa_cv, kappa_fd, phi, x, h)
    ok = findall(isfinite, kappa_cv)
    w = geo.areas[ok]
    mean_cv = sum(w .* kappa_cv[ok]) / sum(w)
    std_cv = sqrt(max(sum(w .* (kappa_cv[ok] .- mean_cv) .^ 2) / sum(w), 0.0))
    band = abs.(phi) .<= 0.8 * h
    kfd = kappa_fd[band]
    mean_fd = sum(kfd) / length(kfd)
    onfacets = [sampleTrilinear(kappa_fd, x, geo.centroids[q]) for q in ok]
    rms = sqrt(sum(w .* (onfacets .- kappa_cv[ok]) .^ 2) / sum(w))
    return (kappa_cv=mean_cv, kappa_cv_std=std_cv, kappa_cv_nan=length(kappa_cv) - length(ok),
            kappa_fd=mean_fd, kappa_cv_min=minimum(kappa_cv[ok]),
            kappa_cv_max=maximum(kappa_cv[ok]), rms_cv_fd=rms)
end

# ---------------------------------------------------------------------------------------
# Self-checks (Stage C)
# ---------------------------------------------------------------------------------------

"""
    verifyTransportOrdering(n) -> max error

A planar level set phi = x - x0 advected by a constant v = (vx, 0, 0) satisfies
phi(t+dt) = x - x0 - vx*dt exactly (a linear field is in WENO5's reconstruction space), so a
transposed velocity array or a wrong fastest axis shows up as a large error.
"""
function verifyTransportOrdering()
    n = 16
    x = collect(range(-0.5, 0.5; length=n))
    h = x[2] - x[1]
    field = Field(Discretization(Dict(:init => Grid((x, x, x)))),
                  Dict(:init => [a - 0.1 for a in x, _ in x, _ in x]))
    vx = 0.04
    velocity = fill(SVector(vx, 0.0, 0.0), n^3)
    dt = 0.2 * h / vx
    eulerianTransport!(field, velocity, dt; level_id=:init, cache=WENO5Cache())
    c = n ÷ 2
    return abs(field.data[:init][c, c, c] - (x[c] - 0.1 - vx * dt))
end

"""
    verifyExtension(cfg, results) -> rows

Independent test of the facet-to-node extension on a **known surface function**, using the
analytic sphere as ground truth:

  * scalar g_f = 1 on every facet. The exact extension of the constant 1 is 1 at every point
    whose nearest interface point exists, so the nodal field must be 1 in the band; the test
    reports max |value - 1|.
  * scalar g_f = x-coordinate of the facet centroid. Its exact closest-point extension at a
    node p is the x-coordinate of the closest point of the *smooth* sphere, 0.23 * p_x/|p|;
    the test reports the max and rms deviation from that closed form. The measured error
    therefore combines the O(h^2) distance between the faceted and the smooth surface with
    the O(h^2) centroid-value difference -- it is a property of the discretisation, not of
    the search, so it must shrink with n.
  * vector g_f = n_f (via `facetVelocity!` with kappa = -1/mu, so the rate is +1): the nodal
    velocity must be a unit vector and must agree with the analytic normal p/|p| up to the
    facet normals' O(h) error, which is why that error is expected to be much larger than the
    scalar one and to converge only at first order.
"""
function verifyExtension(cfg)
    rows = NamedTuple[]
    for n in (28, 32, 40)
        x = collect(range(-0.5, 0.5; length=n))
        h = x[2] - x[1]
        phi = sphereField(x, cfg.R0)
        field = Field(Discretization(Dict(:init => Grid((x, x, x)))), Dict(:init => phi))
        geo = interfaceData(field, x)
        nc = n - 1
        index = FacetIndex(geo.cellof, nc * nc * nc)
        facet_of = zeros(Int32, size(phi))
        cfgband = Config(n=n, band_out=cfg.band_out, search_radius=cfg.search_radius)
        info = assignNearestFacets!(facet_of, phi, x, geo, index, cfgband)
        bandmask = facet_of .> 0

        ones_f = ones(Float64, geo.nfacets)
        nodal_const = facetValuesToNodes(geo, ones_f, facet_of, phi)
        err_const = maximum(abs.(nodal_const[bandmask] .- 1))

        xval = [c[1] for c in geo.centroids]
        nodal_x = facetValuesToNodes(geo, xval, facet_of, phi)
        exact_x = [cfg.R0 * x[i] / norm(SVector(x[i], x[j], x[k]))
                   for i in 1:n, j in 1:n, k in 1:n]
        dev = abs.(nodal_x .- exact_x)
        err_x_max = maximum(dev[bandmask])
        err_x_rms = sqrt(sum(dev[bandmask] .^ 2) / count(bandmask))

        velocity = fill(SVector(0.0, 0.0, 0.0), n^3)
        kappa_unit = fill(-1 / cfg.mu, geo.nfacets)     # rate = -mu*kappa = +1
        vinfo = facetVelocity!(velocity, phi, x, geo, kappa_unit, facet_of, cfgband)
        vnorm_err = 0.0
        vdir_err = 0.0
        for k in 1:n, j in 1:n, i in 1:n
            facet_of[i, j, k] > 0 || continue
            v = velocity[i + (j - 1) * n + (k - 1) * n * n]
            d = abs(phi[i, j, k]) <= cfg.band_in * h ? 1.0 : 0.0
            d == 1.0 || continue
            vnorm_err = max(vnorm_err, abs(norm(v) - 1))
            p = SVector(Float64(x[i]), Float64(x[j]), Float64(x[k]))
            vdir_err = max(vdir_err, norm(v .- p / norm(p)))
        end
        push!(rows, (n=n, h=h, nfacets=geo.nfacets, band=info.band,
                     unresolved=info.unresolved, err_const=err_const,
                     err_x_max=err_x_max, err_x_rms=err_x_rms,
                     vnorm_err=vnorm_err, vdir_err=vdir_err,
                     velocity_nodes=vinfo.nodes))
        @printf("Stage C  n=%d  band=%d  |g-1|max=%.2e  |g_x-x*|max=%.2e rms=%.2e  |v-1|max=%.2e  |v-n|max=%.2e\n",
                n, info.band, err_const, err_x_max, err_x_rms, vnorm_err, vdir_err)
    end
    return rows
end

# ---------------------------------------------------------------------------------------
# Driver
# ---------------------------------------------------------------------------------------

"""
    runCase(cfg, geometry, method, outdir) -> summary

Advance one level set from t = 0 to tmax with `method = :cv` (Curvature.jl) or `:fd` (nodal
finite differences) and write the time series, snapshots and cost breakdown to `outdir`.

Timings are accumulated per operator and measured with `time_ns` around each call, after a
warm-up cut and curvature call so the first-touch compilation is not charged to the loop.
"""
function runCase(cfg, geometry::Symbol, method::Symbol, outdir)
    n, mu = cfg.n, cfg.mu
    x = collect(range(-0.5, 0.5; length=n))
    h = x[2] - x[1]
    phi0 = geometry === :A ? sphereField(x, cfg.R0) : perturbedField(x, cfg.R0, cfg.epsilon)
    field = Field(Discretization(Dict(:init => Grid((x, x, x)))), Dict(:init => copy(phi0)))
    phi = field.data[:init]
    nc = n - 1
    dt = cfg.dt_coefficient * h^2 / mu
    nsteps = max(1, ceil(Int, cfg.tmax / dt))
    dt = cfg.tmax / nsteps                      # exact final time

    velocity = fill(SVector(0.0, 0.0, 0.0), n^3)
    normals = Array{SVector{3,Float64}}(undef, n, n, n)
    kappa_fd = zeros(Float64, n, n, n)
    facet_of = zeros(Int32, n, n, n)
    cache = WENO5Cache()

    mkpath(outdir)
    cost = Dict{Symbol,Float64}(:cut => 0.0, :curvature => 0.0, :extension => 0.0,
                                :fd => 0.0, :transport => 0.0, :reinit => 0.0)
    clock(f) = (t0 = time_ns(); r = f(); (r, (time_ns() - t0) / 1e9))

    # Warm-up: one cut, one curvature call, one transport step at t = 0 (not timed).
    warm = interfaceData(field, x)
    if method === :cv
        facetCurvature(warm, x, cfg)
    end
    computeNormals!(normals, phi, h, cfg.eps_gradient)
    computeCurvature!(kappa_fd, normals, h)
    nodalVelocity!(velocity, phi, kappa_fd, normals, h, cfg)
    eulerianTransport!(field, velocity, 0.0; level_id=:init, cache=cache)

    rows = NamedTuple[]
    snapshots = NamedTuple[]
    info_unreached = 0        # band nodes with no facet in the searched block (cv only)
    reinit_log = NamedTuple[]
    t = 0.0
    step = 0
    targets = cfg.snapshot_times
    next_target = 1

    function snapshot(tag)
        geo, t_cut = clock(() -> interfaceData(field, x))
        cost[:cut] += t_cut
        kcv, t_kcv = clock(() -> facetCurvature(geo, x, cfg))
        cost[:curvature] += t_kcv
        # The nodal finite-difference curvature is recomputed here for both methods: in the
        # `fd` run it is the driver, in the `cv` run it is only this comparison field.
        kfd = kappa_fd
        _, t_fd = clock(() -> begin
            computeNormals!(normals, phi, h, cfg.eps_gradient)
            computeCurvature!(kfd, normals, h)
        end)
        cost[:fd] += t_fd
        cs = curvatureSummary(geo, kcv, kfd, phi, x, h)
        gstat = gradientBandStats(phi, x, h)
        crossings = [axisCrossing(phi, x, u, s) for u in 1:3 for s in (-1, 1)]
        row = (tag=tag, t=t, nfacets=geo.nfacets, radius=geo.radius,
               radius_exact=geometry === :A ? exactRadius(cfg, t) : NaN,
               volume=geo.volume,
               volume_exact=geometry === :A ? exactVolume(cfg, t) : NaN,
               area=geo.area, sphericity=geo.sphericity,
               crossing_z=0.5 * (crossings[5] + crossings[6]),
               crossing_x=0.5 * (crossings[1] + crossings[2]),
               crossing_y=0.5 * (crossings[3] + crossings[4]),
               kappa_cv=cs.kappa_cv, kappa_cv_std=cs.kappa_cv_std,
               kappa_fd=cs.kappa_fd, rms_cv_fd=cs.rms_cv_fd,
               grad_mean=gstat[1], grad_min=gstat[2], grad_max=gstat[3],
               kappa_cv_nan=cs.kappa_cv_nan,
               # band nodes the last extension could not resolve; 0 in every run of this report
               unreached=info_unreached)
        push!(snapshots, row)
        # geometry for the visualiser (facets + curvatures) and the field for cross-sections
        open(joinpath(outdir, "facets_$tag.csv"), "w") do io
            println(io, "# x,y,z,area,kappa_cv,kappa_fd")
            for q in 1:geo.nfacets
                c = geo.centroids[q]
                @printf(io, "%.10e,%.10e,%.10e,%.10e,%.10e,%.10e\n", c[1], c[2], c[3],
                        geo.areas[q], kcv[q], sampleTrilinear(kfd, x, c))
            end
        end
        open(joinpath(outdir, "phi_$tag.bin"), "w") do io
            write(io, vec(Float32.(phi)))
        end
        return row
    end

    push!(rows, snapshot("t0000"))
    # One cut per step: the cut taken at the end of a step is the geometry the next step's
    # velocity is built from, so the recorded row and the advected state always agree.
    geo, t_cut0 = clock(() -> interfaceData(field, x))
    cost[:cut] += t_cut0
    @printf("%s %s n=%d h=%.4f dt=%.5f steps=%d\n", geometry, method, n, h, dt, nsteps)

    for _ in 1:nsteps
        dt_step = min(dt, cfg.tmax - t)
        t0 = time_ns()
        if method === :cv
            kappa, t_k = clock(() -> facetCurvature(geo, x, cfg))
            cost[:curvature] += t_k
            index = FacetIndex(geo.cellof, nc * nc * nc)
            info, t_e = clock(() -> begin
                a = assignNearestFacets!(facet_of, phi, x, geo, index, cfg)
                b = facetVelocity!(velocity, phi, x, geo, kappa, facet_of, cfg)
                (a, b)
            end)
            cost[:extension] += t_e
            vec_info = info[2]
            info_unreached = info[1].unresolved
        else
            _, t_fd = clock(() -> begin
                computeNormals!(normals, phi, h, cfg.eps_gradient)
                computeCurvature!(kappa_fd, normals, h)
            end)
            cost[:fd] += t_fd
            vec_info = (nodes=nodalVelocity!(velocity, phi, kappa_fd, normals, h, cfg),
                        nonfinite=0)
        end
        _, t_tr = clock(() -> eulerianTransport!(field, velocity, dt_step;
                                                 level_id=:init, cache=cache))
        cost[:transport] += t_tr
        t += dt_step
        step += 1

        geo, t_c = clock(() -> interfaceData(field, x))
        cost[:cut] += t_c

        if cfg.reinit_every > 0 && step % cfg.reinit_every == 0
            # `geo` is cut #1 (before), the Cut passed in is cut #2 (independent construction of
            # the same field), `after` is cut #3 (after): three separately built snapshots.
            cut_ids, t_b = clock(() -> cut(field, :init))
            _, t_r = clock(() -> reinitialize!(field, cut_ids))
            after, t_c2 = clock(() -> interfaceData(field, x))
            cost[:cut] += t_b + t_c2
            cost[:reinit] += t_r
            push!(reinit_log, (t=t, step=step,
                               radius_before=geo.radius, radius_after=after.radius,
                               volume_before=geo.volume, volume_after=after.volume,
                               area_before=geo.area, area_after=after.area,
                               nfacets_before=geo.nfacets, nfacets_after=after.nfacets))
            geo = after
        end

        push!(rows, (tag="t$(lpad(step, 4, '0'))", t=t, nfacets=geo.nfacets,
                     radius=geo.radius,
                     radius_exact=geometry === :A ? exactRadius(cfg, t) : NaN,
                     volume=geo.volume,
                     volume_exact=geometry === :A ? exactVolume(cfg, t) : NaN,
                     area=geo.area, sphericity=geo.sphericity,
                     crossing_z=NaN, crossing_x=NaN, crossing_y=NaN,
                     kappa_cv=NaN, kappa_cv_std=NaN, kappa_fd=NaN, rms_cv_fd=NaN,
                     grad_mean=gradientBandStats(phi, x, h)[1], grad_min=NaN, grad_max=NaN,
                     kappa_cv_nan=NaN, unreached=Float64(vec_info.nonfinite)))
        if next_target <= length(targets) && t >= targets[next_target] - 1e-12
            snapshot("t" * @sprintf("%04d", next_target))
            next_target += 1
        end
        step_time = (time_ns() - t0) / 1e9
        @printf("  step %3d  t=%.4f  R=%.6f  V=%.6e  band=%d  %.2fs\n", step, t,
                geo.radius, geo.volume, vec_info.nodes, step_time)
    end
    final = snapshot("final")

    # CSV: per-step radius/volume/area from the cut at the end of each step.
    open(joinpath(outdir, "timeseries.csv"), "w") do io
        println(io, "step,t,nfacets,radius,radius_exact,volume,volume_exact,area,sphericity,grad_band,nonfinite_velocity")
        for (i, r) in enumerate(rows)
            @printf(io, "%d,%.10e,%d,%.10e,%.10e,%.10e,%.10e,%.10e,%.10e,%.10e,%.0f\n", i - 1, r.t,
                    r.nfacets, r.radius, r.radius_exact, r.volume, r.volume_exact, r.area,
                    r.sphericity, r.grad_mean, isnan(r.unreached) ? 0.0 : r.unreached)
        end
    end
    open(joinpath(outdir, "snapshots.csv"), "w") do io
        println(io, "tag,t,n,nfacets,radius,radius_exact,volume,volume_exact,area,sphericity," *
                    "crossing_x,crossing_y,crossing_z,kappa_cv,kappa_cv_std,kappa_fd," *
                    "rms_cv_fd,grad_mean,grad_min,grad_max,shrink,shrink_exact," *
                    "kappa_cv_nan,band_nodes_unreached")
        r0 = rows[1].radius
        for r in snapshots
            @printf(io, "%s,%.10e,%d,%d,%.10e,%.10e,%.10e,%.10e,%.10e,%.10e,%.10e,%.10e,%.10e,%.10e,%.10e,%.10e,%.10e,%.10e,%.10e,%.10e,%.10e,%.10e,%d,%d\n",
                    r.tag, r.t, n, r.nfacets, r.radius, r.radius_exact, r.volume,
                    r.volume_exact, r.area, r.sphericity, r.crossing_x, r.crossing_y,
                    r.crossing_z, r.kappa_cv, r.kappa_cv_std, r.kappa_fd, r.rms_cv_fd,
                    r.grad_mean, r.grad_min, r.grad_max, r0 - r.radius,
                    geometry === :A ? cfg.R0 - r.radius_exact : NaN,
                    r.kappa_cv_nan, r.unreached)
        end
    end
    open(joinpath(outdir, "reinit.csv"), "w") do io
        println(io, "step,t,radius_before,radius_after,volume_before,volume_after," *
                    "area_before,area_after,nfacets_before,nfacets_after")
        for r in reinit_log
            @printf(io, "%d,%.10e,%.10e,%.10e,%.10e,%.10e,%.10e,%.10e,%d,%d\n", r.step, r.t,
                    r.radius_before, r.radius_after, r.volume_before, r.volume_after,
                    r.area_before, r.area_after, r.nfacets_before, r.nfacets_after)
        end
    end
    total_cost = sum(values(cost))
    open(joinpath(outdir, "cost.csv"), "w") do io
        println(io, "component,seconds,fraction,per_step_ms")
        for k in (:cut, :curvature, :extension, :fd, :transport, :reinit)
            @printf(io, "%s,%.6f,%.4f,%.3f\n", k, cost[k], cost[k] / total_cost,
                    1e3 * cost[k] / nsteps)
        end
        @printf(io, "total,%.6f,1.0,%.3f\n", total_cost, 1e3 * total_cost / nsteps)
    end

    summary = (geometry=geometry, method=method, n=n, h=h, dt=dt, nsteps=nsteps,
               tmax=cfg.tmax, final=final, cost=cost, total_cost=total_cost,
               reinit=reinit_log, snapshots=snapshots, outdir=outdir)
    writeSummary(summary, cfg)
    return summary
end

"writeSummary(s, cfg) — human-readable summary of one (geometry, method) run."
function writeSummary(s, cfg)
    f = s.final
    open(joinpath(s.outdir, "summary.md"), "w") do io
        println(io, "# CurvatureExtension3D — geometry ", s.geometry, ", method ", s.method)
        println(io, "\n- grid: ", s.n, "^3 nodes, h = ", round(s.h; sigdigits=6))
        println(io, "- dt = ", round(s.dt; sigdigits=6), " (C = ", cfg.dt_coefficient,
                ", dt = C h^2/mu), steps = ", s.nsteps, ", tmax = ", s.tmax)
        println(io, "- reinitialisation every ", cfg.reinit_every, " steps; band |phi| <= ",
                cfg.band_out, "h; Curvature.jl window ", cfg.curvature_radius, "h")
        println(io, "\n## Final interface\n")
        @printf(io, "| quantity | value |\n|---|---|\n")
        @printf(io, "| facets | %d |\n", f.nfacets)
        @printf(io, "| area-weighted radius | %.8f |\n", f.radius)
        if isfinite(f.radius_exact)
            @printf(io, "| exact radius R(t) | %.8f |\n", f.radius_exact)
            @printf(io, "| radius relative error | %.3e |\n", (f.radius - f.radius_exact) / f.radius_exact)
            @printf(io, "| volume | %.8e |\n| exact volume | %.8e |\n", f.volume, f.volume_exact)
            @printf(io, "| volume relative error | %.3e |\n", (f.volume - f.volume_exact) / f.volume_exact)
        else
            @printf(io, "| volume | %.8e |\n", f.volume)
        end
        @printf(io, "| surface area | %.8e |\n", f.area)
        @printf(io, "| sphericity (std R / R) | %.3e |\n", f.sphericity)
        @printf(io, "| crossings x,y,z | %.6f, %.6f, %.6f |\n", f.crossing_x, f.crossing_y, f.crossing_z)
        @printf(io, "| kappa (Curvature.jl, area mean) | %.6f +- %.6f |\n", f.kappa_cv, f.kappa_cv_std)
        @printf(io, "| kappa (nodal FD, band mean) | %.6f |\n", f.kappa_fd)
        @printf(io, "| rms CV-vs-FD at facet centroids | %.6f |\n", f.rms_cv_fd)
        @printf(io, "| |grad phi| in band (mean/min/max) | %.4f / %.4f / %.4f |\n",
                f.grad_mean, f.grad_min, f.grad_max)
        if isfinite(f.radius_exact)
            @printf(io, "| exact mean curvature 2/R(t) | %.6f |\n", 2 / f.radius_exact)
        end
        println(io, "\n## Cost (seconds, warm)\n")
        println(io, "| component | total | fraction | ms/step |\n|---|---|---|---|")
        for k in (:cut, :curvature, :extension, :fd, :transport, :reinit)
            @printf(io, "| %s | %.3f | %.1f%% | %.2f |\n", k, s.cost[k],
                    100 * s.cost[k] / s.total_cost, 1e3 * s.cost[k] / s.nsteps)
        end
        @printf(io, "| total | %.3f | 100%% | %.2f |\n", s.total_cost,
                1e3 * s.total_cost / s.nsteps)
        if !isempty(s.reinit)
            println(io, "\n## Reinitialisation (independent cuts before and after)\n")
            println(io, "| step | t | R before | R after | dR | dV/V |\n|---|---|---|---|---|---|")
            for r in s.reinit
                @printf(io, "| %d | %.4f | %.8f | %.8f | %.2e | %.2e |\n", r.step, r.t,
                        r.radius_before, r.radius_after,
                        abs(r.radius_after - r.radius_before),
                        abs(r.volume_after - r.volume_before) / r.volume_before)
            end
        end
    end
end

# ---------------------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------------------

function main(args)
    n = isempty(args) ? 32 : parse(Int, args[1])
    tmax = length(args) < 2 ? 0.12 : parse(Float64, args[2])
    root = @__DIR__
    cfg = Config(n=n, tmax=tmax,
                 snapshot_times=sort(unique([0.25, 0.5, 0.75, 1.0] .* tmax)))

    @printf("=== CurvatureExtension3D — n=%d, tmax=%.4f, mu=%.3f ===\n", n, tmax, cfg.mu)
    order_err = verifyTransportOrdering()
    @printf("Stage C  transport ordering self-check: max error %.3e (must be ~1e-16)\n", order_err)
    order_err < 1e-12 || error("transport velocity ordering self-check failed")

    stageC = verifyExtension(cfg)
    mkpath(joinpath(root, "output", "stageC"))
    open(joinpath(root, "output", "stageC", "extension_test.csv"), "w") do io
        println(io, "n,h,nfacets,band,unresolved,err_const,err_x_max,err_x_rms,vnorm_err,vdir_err,velocity_nodes")
        for r in stageC
            @printf(io, "%d,%.6e,%d,%d,%d,%.6e,%.6e,%.6e,%.6e,%.6e,%d\n", r.n, r.h,
                    r.nfacets, r.band, r.unresolved, r.err_const, r.err_x_max, r.err_x_rms,
                    r.vnorm_err, r.vdir_err, r.velocity_nodes)
        end
    end

    summaries = []
    for geometry in (:A, :B), method in (:cv, :fd)
        # `fd` does not need Curvature.jl to drive, but its snapshots still record it.
        outdir = joinpath(root, "output", "n$(n)", string(geometry), string(method))
        push!(summaries, runCase(cfg, geometry, method, outdir))
    end

    # Machine-readable cross-method comparison.
    open(joinpath(root, "output", "n$(n)", "comparison.csv"), "w") do io
        println(io, "geometry,method,t,radius,radius_exact,volume,volume_exact,area," *
                    "sphericity,crossing_z,kappa_cv,kappa_fd,rms_cv_fd,shrink,total_cost,ms_per_step")
        for s in summaries, r in s.snapshots
            @printf(io, "%s,%s,%.10e,%.10e,%.10e,%.10e,%.10e,%.10e,%.10e,%.10e,%.10e,%.10e,%.10e,%.10e,%.6f,%.3f\n",
                    s.geometry, s.method, r.t, r.radius, r.radius_exact, r.volume,
                    r.volume_exact, r.area, r.sphericity, r.crossing_z, r.kappa_cv,
                    r.kappa_fd, r.rms_cv_fd, s.snapshots[1].radius - r.radius, s.total_cost,
                    1e3 * s.total_cost / s.nsteps)
        end
    end

    println("\n=== done; outputs under ", joinpath(root, "output", "n$n"))
    return summaries
end

summaries = main(ARGS)

# Figures are produced by a separate script that only reads the outputs above.
if get(ENV, "CURVATUREEXT_SKIP_FIGURES", "0") != "1"
    visualizer = joinpath(@__DIR__, "visualize.CurvatureExtension3D.jl")
    if isfile(visualizer)
        include(visualizer)
        Base.invokelatest(CurvatureExtensionFigures.makeFigures, joinpath(@__DIR__, "output"))
    end
end
