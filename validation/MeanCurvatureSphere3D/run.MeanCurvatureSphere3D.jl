#!/usr/bin/env julia
"""
Mean curvature flow of a sphere in D = [-1/2, 1/2]^3.

    V_n = -mu * kappa,        kappa = div n,   n = grad(phi) / |grad(phi)|

is solved through the level-set transport of EvolvingDomains. The interface is a sphere of
initial radius R0, so the exact solution is R(t) = sqrt(R0^2 - 4 mu t).

    julia --project=../env run.MeanCurvatureSphere3D.jl          # 48^3, the default case
    julia --project=../env run.MeanCurvatureSphere3D.jl 64       # finer case

What is computed and what is not
--------------------------------
At every step the velocity driving the transport is built from the *discrete* curvature of
the current level set. Nothing is read from the analytic radius: R(t) enters only in
`diagnostics`, as the reference the result is compared with.

Curvature: Curvature.jl is used through its public interface for validation (it consumes the
public `CutSurface` geometry of an EvolvingDomains `Cut`). It cannot drive the transport: it
returns one curvature per interface facet, not a nodal field over the grid, and it has no
public extension operator, so a per-step use would need a nearest-facet extension and one
fit per facet, i.e. seconds per step. The transport therefore uses the finite-difference
operator `computeCurvature!` below — a documented fallback, cross-checked against
Curvature.jl at every snapshot (`fd_vs_cv_rms` column, curvature figure).

Data ordering: EvolvingDomains wants one `SVector{3,Float64}` per grid node, first axis
fastest, i.e. `velocity[i + (j-1)*n + (k-1)*n^2]` for `phi[i, j, k]`.
"""

using Libdl: dlext
using Printf

# ---------------------------------------------------------------------------------------
# Native cut library. This checkout has no published artifact, so the documented developer
# fallback (scripts/build_cxx_bindings.sh) has to be pointed at. The environment must be set
# before `using EvolvingDomains`, and LOAD_PATH before `using Curvature`.
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

ENV["EVOLVINGDOMAINS_CUTCELLS_LIBRARY"] = locateCutLibrary()
# Curvature.jl is a sibling source checkout; putting its project on LOAD_PATH loads it with
# its own dependency graph. No project or manifest of either package is modified.
push!(LOAD_PATH, normpath(joinpath(@__DIR__, "..", "..", "..", "Curvature.jl")))

using EvolvingDomains
using StaticArrays: SVector
using LinearAlgebra: norm, dot, cross
using Curvature: interfaceCurvature      # public: facet curvature from plain geometry


# ---------------------------------------------------------------------------------------
# Configuration
# ---------------------------------------------------------------------------------------

Base.@kwdef struct Config
    n::Int = 48                        # nodes per axis on a uniform cubic grid
    R0::Float64 = 0.27
    mu::Float64 = 0.01
    tmax::Float64 = 0.6
    # Explicit diffusion (mu) in 3D is stable for mu*dt/h^2 <= 1/6; 1/12 is a factor-two
    # safety margin: dt <= C*h^2/mu with C = dt_coefficient.
    dt_coefficient::Float64 = 1 / 12
    reinit_every::Int = 10             # steps between reinitialize! calls; 0 disables
    band_in::Float64 = 1.5             # full interface speed for |phi| <= band_in*h
    band_out::Float64 = 4.0            # velocity tapered to zero at |phi| = band_out*h
    eps_gradient::Float64 = 0.5        # |grad phi| clamp: below this, n is not defined
    curvature_radius::Float64 = 4.0    # Curvature.jl fit window, in cells
    coordinates::Vector{Float64} = [0.0, 0.2, 0.4, 0.6]
end

exactRadius(cfg, t) = sqrt(max(cfg.R0^2 - 4 * cfg.mu * t, 0.0))
exactVolume(cfg, t) = 4π * exactRadius(cfg, t)^3 / 3

# ---------------------------------------------------------------------------------------
# Geometry: normals and mean curvature of the level set
# ---------------------------------------------------------------------------------------

"""
    computeNormals!(normals, phi, h, eps) -> clamped

Outward unit normals n = grad(phi)/|grad(phi)| by central differences, clamped indices at
the box boundary. Where |grad(phi)| < eps the direction is kept and the magnitude stored is
|n| = |grad|/eps <= 1; the returned count says how many nodes were clamped. eps = 0.5 is a
length-scale match, not a tuning knob: a reinitialised level set has |grad phi| = 1 near the
interface, and clamping below half of that leaves the interface untouched.
"""
function computeNormals!(normals, phi, h, eps)
    nx, ny, nz = size(phi)
    inv2h = 1 / (2h)
    clamped = 0
    @inbounds for k in 1:nz, j in 1:ny, i in 1:nx
        gx = (phi[clamp(i + 1, 1, nx), j, k] - phi[clamp(i - 1, 1, nx), j, k]) * inv2h
        gy = (phi[i, clamp(j + 1, 1, ny), k] - phi[i, clamp(j - 1, 1, ny), k]) * inv2h
        gz = (phi[i, j, clamp(k + 1, 1, nz)] - phi[i, j, clamp(k - 1, 1, nz)]) * inv2h
        m = sqrt(gx * gx + gy * gy + gz * gz)
        if m > eps
            s = 1 / m
        else
            s = 1 / eps
            clamped += 1
        end
        normals[i, j, k] = SVector(gx * s, gy * s, gz * s)
    end
    return clamped
end

"""
    computeCurvature!(kappa, normals, h) -> nonfinite

Discrete mean curvature kappa = div n by central differences of the nodal normal field. This
is the fallback operator that drives the transport: Curvature.jl gives per-facet values on
the interface, not a field over the grid. Away from the interface the value is meaningless
(the normal field there is not a surface normal), which is why the velocity is only used
inside the band. Returns the number of non-finite entries (expected 0).
"""
function computeCurvature!(kappa, normals, h)
    nx, ny, nz = size(kappa)
    inv2h = 1 / (2h)
    @inbounds for k in 1:nz, j in 1:ny, i in 1:nx
        a = normals[clamp(i + 1, 1, nx), j, k]
        b = normals[clamp(i - 1, 1, nx), j, k]
        c = normals[i, clamp(j + 1, 1, ny), k]
        d = normals[i, clamp(j - 1, 1, ny), k]
        e = normals[i, j, clamp(k + 1, 1, nz)]
        f = normals[i, j, clamp(k - 1, 1, nz)]
        kappa[i, j, k] = ((a[1] - b[1]) + (c[2] - d[2]) + (e[3] - f[3])) * inv2h
    end
    return count(!isfinite, kappa)
end

"""
    dynamics(cfg)

The geometric dynamics alone: the normal rate F = -mu*kappa, i.e. V_n = -mu*kappa with the
outward normal of `phi`. Positive kappa (convex, outward normal) therefore gives an inward
velocity and a shrinking sphere.
"""
dynamics(cfg) = kappa -> -cfg.mu .* kappa

"""
    velocityFromRates!(velocity, phi, rate, normals, h, cfg) -> full

Narrow-band extension of the interface velocity v = rate * n to the whole grid: full speed
for |phi| <= band_in*h, a C^1 cosine taper to zero at |phi| = band_out*h. The zero level set
lies in the full-speed zone (it moves ~0.15 h per step), so on Γ the transport sees exactly
v·∇φ = -mu*kappa*|∇φ|; the taper only holds the inert part of the level set still. Returns
the number of nodes carrying the full interface speed.

shortcut: the velocity is tapered rather than extended along normals beyond the band, which
is adequate while the band is wider than the per-step interface displacement; upgrade to a
normal (PDE) extension if that stops holding.
"""
function velocityFromRates!(velocity, phi, rate, normals, h, cfg)
    nx, ny, nz = size(phi)
    rin, rout = cfg.band_in * h, cfg.band_out * h
    full = 0
    @inbounds for k in 1:nz, j in 1:ny, i in 1:nx
        d = abs(phi[i, j, k])
        if d <= rin
            chi = 1.0
            full += 1
        elseif d >= rout
            chi = 0.0
        else
            c = cos(π / 2 * (d - rin) / (rout - rin))
            chi = c * c
        end
        n = normals[i, j, k]
        s = chi * rate[i, j, k]
        velocity[i + (j - 1) * nx + (k - 1) * nx * ny] = SVector(s * n[1], s * n[2], s * n[3])
    end
    return full
end

"""Advance phi through EvolvingDomains: WENO5 + SSP-RK3 Eulerian advection, first axis fastest."""
advanceGeometry!(field, velocity, dt, cache) =
    eulerianTransport!(field, velocity, dt; level_id=:init, cache=cache)

"""
    verifyTransportOrdering() -> error

Self-check of the one thing the type system cannot check: the node ordering of the velocity
array. A planar level set phi = x - x0 advected by a constant v = (vx, vy, vz) satisfies
phi(x, t+dt) = x - x0 - vx*dt exactly (the linear field is in WENO5's reconstruction space),
so an x-axis component applied to the y index, a transposed array, or a wrong fastest axis
all show up as a large error at one interior node.
"""
function verifyTransportOrdering()
    n = 16
    x = collect(range(-0.5, 0.5; length=n))
    h = x[2] - x[1]
    field = Field(Discretization(Dict(:init => Grid((x, x, x)))),
                  Dict(:init => [a - 0.1 for a in x, _ in x, _ in x]))
    vx, vy, vz = 0.04, -0.03, 0.02
    velocity = fill(SVector(vx, vy, vz), n^3)
    dt = 0.2 * h / vx                     # displacement 0.2 h along x
    eulerianTransport!(field, velocity, dt; level_id=:init, cache=WENO5Cache())
    i, j, k = n ÷ 2, n ÷ 2, n ÷ 2
    return abs(field.data[:init][i, j, k] - (x[i] - 0.1 - vx * dt))
end

# ---------------------------------------------------------------------------------------
# Interface extraction (public CutSurface -> plain geometry)
# ---------------------------------------------------------------------------------------

"""
    interfaceGeometry(field, h)

Triangulated zero level set of `field` from the public `CutSurface` of `cut`, converted to
plain geometry for Curvature.jl:

  * facet `q` owns the CSR row `point_ids[point_pointers[q] : point_pointers[q+1]-1]`, whose
    three entries index `physical_points`;
  * `boundary.normals[q]` is that facet's unit normal, already oriented out of the phi < 0
    phase, i.e. outward for this sphere.

Returns the facets, their normals, centroids and areas, the area-weighted mean radius with
its standard deviation, and the enclosed volume from the divergence theorem,
V = (1/3) * sum(area * centroid . normal). These are zero-crossing estimates of the discrete
interface, not analytic values.
"""
function interfaceGeometry(field, h)
    boundary = cut(field, :init).boundary
    nf = length(boundary.facet_to_bgcell)
    facets = Vector{NTuple{3,SVector{3,Float64}}}(undef, nf)
    normals = Vector{SVector{3,Float64}}(undef, nf)
    centroids = Vector{SVector{3,Float64}}(undef, nf)
    areas = Vector{Float64}(undef, nf)
    @inbounds for q in 1:nf
        a1 = Int(boundary.point_pointers[q])
        a2 = Int(boundary.point_pointers[q + 1]) - 1
        a2 - a1 + 1 == 3 || error("facet $q has $(a2 - a1 + 1) vertices, expected 3")
        p1 = SVector{3,Float64}(boundary.physical_points[Int(boundary.point_ids[a1])])
        p2 = SVector{3,Float64}(boundary.physical_points[Int(boundary.point_ids[a1 + 1])])
        p3 = SVector{3,Float64}(boundary.physical_points[Int(boundary.point_ids[a1 + 2])])
        facets[q] = (p1, p2, p3)
        normals[q] = SVector{3,Float64}(boundary.normals[q])
        centroids[q] = (p1 + p2 + p3) / 3
        areas[q] = norm(cross(p2 - p1, p3 - p1)) / 2
    end
    total = sum(areas)
    mean_radius = sum(areas .* norm.(centroids)) / total
    variance = max(sum(areas .* (norm.(centroids) .- mean_radius) .^ 2) / total, 0.0)
    volume = sum(areas[q] * dot(centroids[q], normals[q]) for q in 1:nf) / 3
    return (facets=facets, normals=normals, centroids=centroids, areas=areas,
            radius=mean_radius, std=sqrt(variance), volume=volume, nfacets=nf)
end

# ---------------------------------------------------------------------------------------
# Numerically justified point values and axis crossings
# ---------------------------------------------------------------------------------------

"""
    sampleTrilinear(phi, x, p)

Trilinear interpolation of a nodal field on the uniform cubic grid `x` (all axes equal, as
here). Used to evaluate the discrete curvature at a facet centroid and to locate
coordinate-axis crossings of phi.
"""
function sampleTrilinear(phi, x, p)
    n = length(x)
    h = x[2] - x[1]
    s = ntuple(d -> clamp((p[d] - x[1]) / h + 1, 1.0, float(n)), 3)
    lo = ntuple(d -> clamp(floor(Int, s[d]), 1, n - 1), 3)
    fx, fy, fz = ntuple(d -> s[d] - lo[d], 3)
    i, j, k = lo
    c00 = phi[i, j, k] * (1 - fx) + phi[i + 1, j, k] * fx
    c10 = phi[i, j + 1, k] * (1 - fx) + phi[i + 1, j + 1, k] * fx
    c01 = phi[i, j, k + 1] * (1 - fx) + phi[i + 1, j, k + 1] * fx
    c11 = phi[i, j + 1, k + 1] * (1 - fx) + phi[i + 1, j + 1, k + 1] * fx
    return (c00 * (1 - fy) + c10 * fy) * (1 - fz) + (c01 * (1 - fy) + c11 * fy) * fz
end

"""
    axisCrossing(phi, x, unit, side)

Zero crossing of the trilinear interpolant of `phi` on the ray from the centre along
`side * e_unit`, by bisection on the sign change. Returns the distance from the centre, or
`NaN` if the ray does not leave the sphere inside the box.
"""
function axisCrossing(phi, x, unit::Int, side::Int)
    at = s -> sampleTrilinear(phi, x, ntuple(d -> d == unit ? side * s : 0.0, 3))
    at(0.0) < 0 || return NaN
    at(0.5) > 0 || return NaN
    lo, hi = 0.0, 0.5
    for _ in 1:80
        mid = (lo + hi) / 2
        at(mid) < 0 ? (lo = mid) : (hi = mid)
    end
    return (lo + hi) / 2
end

# ---------------------------------------------------------------------------------------
# Diagnostics
# ---------------------------------------------------------------------------------------

"""
    diagnostics(field, workspace, cfg, t) -> NamedTuple

Everything compared against the exact sphere at time `t`:

  * estimated radius: area-weighted mean distance of the cut-surface centroids (a zero-level
    crossing estimate, not the analytic radius), with its area-weighted standard deviation
    as the sphericity measure;
  * enclosed volume by the divergence theorem on the same triangulation, against
    4*pi*R(t)^3/3;
  * the six coordinate-axis crossings of the trilinear interpolant;
  * average mean curvature from Curvature.jl (area-weighted) and from the finite-difference
    operator in the band, against 2/R(t);
  * the area-weighted rms difference of the two curvature estimates at the facet centroids.

The returned NamedTuple also carries the facet centroids, areas and per-facet curvatures so
the caller can save the interface without a second fit.
"""
function diagnostics(field, workspace, cfg, t)
    phi, x = field.data[:init], workspace.x
    h = x[2] - x[1]
    geo = interfaceGeometry(field, h)
    rex = exactRadius(cfg, t)

    kappa_cv = interfaceCurvature(geo.facets, geo.normals; spacing=h, origin=x[1],
                                  radius=cfg.curvature_radius * h)
    ok = findall(isfinite, kappa_cv)
    w = geo.areas[ok]
    kcv = sum(w .* kappa_cv[ok]) / sum(w)
    kcv_std = sqrt(max(sum(w .* (kappa_cv[ok] .- kcv) .^ 2) / sum(w), 0.0))

    kappa_fd = workspace.curvature
    band = abs.(phi) .<= 0.8 * h
    kfd_band = kappa_fd[band]
    kfd = sum(kfd_band) / length(kfd_band)
    kfd_std = sqrt(sum((kfd_band .- kfd) .^ 2) / length(kfd_band))

    kfd_on_facets = [sampleTrilinear(kappa_fd, x, geo.centroids[q]) for q in ok]
    rms = sqrt(sum(w .* (kfd_on_facets .- kappa_cv[ok]) .^ 2) / sum(w)) / abs(2 / rex)

    crossings = [axisCrossing(phi, x, unit, side) for unit in 1:3 for side in (-1, 1)]
    axis_radius = [(crossings[2d - 1] + crossings[2d]) / 2 for d in 1:3]
    axis_err = maximum(abs.(crossings .- rex))

    return (t=t, exact_radius=rex, radius=geo.radius, radius_std=geo.std,
            radius_err=geo.radius - rex, radius_rel_err=(geo.radius - rex) / rex,
            volume=geo.volume, exact_volume=exactVolume(cfg, t),
            volume_rel_err=(geo.volume - exactVolume(cfg, t)) / exactVolume(cfg, t),
            axis_radius=axis_radius, axis_crossings=crossings, axis_err_max=axis_err,
            kappa_cv=kcv, kappa_cv_std=kcv_std, kappa_cv_nan=length(kappa_cv) - length(ok),
            kappa_fd=kfd, kappa_fd_std=kfd_std, kappa_fd_vs_cv_rms=rms,
            kappa_exact=2 / rex, nfacets=geo.nfacets,
            centroids=geo.centroids[ok], areas=geo.areas[ok], kappa_cv_per_facet=kappa_cv[ok])
end

"""Mean |grad phi| over the band |phi| <= 1.5h — the distance-function check used for reinit."""
function gradientStats(phi, h)
    nx, ny, nz = size(phi)
    inv2h = 1 / (2h)
    total = 0.0
    n = 0
    @inbounds for k in 1:nz, j in 1:ny, i in 1:nx
        abs(phi[i, j, k]) <= 1.5h || continue
        gx = (phi[clamp(i + 1, 1, nx), j, k] - phi[clamp(i - 1, 1, nx), j, k]) * inv2h
        gy = (phi[i, clamp(j + 1, 1, ny), k] - phi[i, clamp(j - 1, 1, ny), k]) * inv2h
        gz = (phi[i, j, clamp(k + 1, 1, nz)] - phi[i, j, clamp(k - 1, 1, nz)]) * inv2h
        total += sqrt(gx * gx + gy * gy + gz * gz)
        n += 1
    end
    return total / n
end

# ---------------------------------------------------------------------------------------
# Output
# ---------------------------------------------------------------------------------------

const TIMESERIES_COLUMNS = ["t", "radius", "radius_std", "radius_err", "radius_rel_err",
                            "exact_radius", "volume", "exact_volume", "volume_rel_err",
                            "axis_rx", "axis_ry", "axis_rz", "axis_err_max",
                            "kappa_fd", "kappa_fd_std", "kappa_cv", "kappa_cv_std",
                            "kappa_exact", "kappa_cv_nan", "fd_vs_cv_rms", "nfacets"]

_fmt(x) = @sprintf("%.10e", Float64(x))

function writeTimeseries(path, rows)
    open(path, "w") do io
        println(io, "# MeanCurvatureSphere3D timeseries; lengths in [-0.5,0.5]^3, kappa = div n")
        println(io, join(TIMESERIES_COLUMNS, ","))
        for d in rows
            cols = [_fmt(d.t), _fmt(d.radius), _fmt(d.radius_std), _fmt(d.radius_err),
                    _fmt(d.radius_rel_err), _fmt(d.exact_radius), _fmt(d.volume),
                    _fmt(d.exact_volume),
                    _fmt(d.volume_rel_err), _fmt(d.axis_radius[1]), _fmt(d.axis_radius[2]),
                    _fmt(d.axis_radius[3]), _fmt(d.axis_err_max), _fmt(d.kappa_fd),
                    _fmt(d.kappa_fd_std), _fmt(d.kappa_cv), _fmt(d.kappa_cv_std),
                    _fmt(d.kappa_exact), string(d.kappa_cv_nan), _fmt(d.kappa_fd_vs_cv_rms),
                    string(d.nfacets)]
            println(io, join(cols, ","))
        end
    end
end

"""
    saveSnapshot(outputdir, tag, row)

Writes the interface facet centroids, areas and Curvature.jl curvatures, plus the level set
as raw Float64 (first axis fastest). The plotting script needs neither EvolvingDomains nor
the native cut library.
"""
function saveSnapshot(outputdir, tag, row, phi, n)
    open(joinpath(outputdir, "interface_t$tag.csv"), "w") do io
        println(io, "x,y,z,area,kappa_cv")
        for q in eachindex(row.areas)
            c = row.centroids[q]
            @printf(io, "%.10e,%.10e,%.10e,%.10e,%.10e\n", c[1], c[2], c[3],
                    row.areas[q], row.kappa_cv_per_facet[q])
        end
    end
    open(joinpath(outputdir, "phi_n$(n)_t$tag.bin"), "w") do io
        write(io, phi)
    end
    return nothing
end

# ---------------------------------------------------------------------------------------
# Simulation
# ---------------------------------------------------------------------------------------

function makeSetup(cfg)
    x = collect(range(-0.5, 0.5; length=cfg.n))
    h = x[2] - x[1]
    grid = Grid((x, x, x))
    disc = Discretization(Dict(:init => grid))
    phi = [sqrt(a^2 + b^2 + c^2) - cfg.R0 for a in x, b in x, c in x]
    field = Field(disc, Dict(:init => phi))
    return x, h, field
end

"""
    simulate(cfg; outputdir)

Explicit time loop. Per step: normals and curvature of the *current* level set, the rate
F = -mu*kappa, the extended velocity, and one Eulerian transport step. Reinitialisation runs
every `reinit_every` steps; its effect on the discrete interface is measured (the zero
crossings are preserved by construction — what it restores is |grad phi| = 1).
"""
function simulate(cfg; outputdir=joinpath(@__DIR__, "output", "n$(cfg.n)"))
    cfg.coordinates[1] == 0.0 || error("the first snapshot must be the initial state at t = 0")
    cfg.coordinates[end] == cfg.tmax ||
        error("the last snapshot must be the final time; got $(cfg.coordinates[end]) vs $(cfg.tmax)")
    mkpath(outputdir)
    x, h, field = makeSetup(cfg)
    phi = field.data[:init]
    workspace = (x=x, normals=Array{SVector{3,Float64}}(undef, cfg.n, cfg.n, cfg.n),
                 curvature=zeros(cfg.n, cfg.n, cfg.n),
                 velocity=Vector{SVector{3,Float64}}(undef, cfg.n^3))
    cache = WENO5Cache()

    dt = cfg.dt_coefficient * h^2 / cfg.mu
    clamped = computeNormals!(workspace.normals, phi, h, cfg.eps_gradient)
    nonfinite = computeCurvature!(workspace.curvature, workspace.normals, h)

    rows = NamedTuple[]
    tags = String[]
    reinit_shifts = NamedTuple{(:radius, :volume),Tuple{Float64,Float64}}[]
    reinit_grad_before = Float64[]
    reinit_grad_after = Float64[]
    step_times = Float64[]

    row0 = diagnostics(field, workspace, cfg, 0.0)
    saveSnapshot(outputdir, "0.000", row0, phi, cfg.n)
    push!(rows, row0)
    push!(tags, "0.000")

    t = 0.0
    step = 0
    tstart = time()
    next_target = 2
    while t < cfg.tmax - 1e-12
        dt_step = min(dt, cfg.tmax - t)
        steptime = @elapsed begin
            rate = dynamics(cfg)(workspace.curvature)
            velocityFromRates!(workspace.velocity, phi, rate, workspace.normals, h, cfg)
            advanceGeometry!(field, workspace.velocity, dt_step, cache)
            t += dt_step
            step += 1

            if cfg.reinit_every > 0 && step % cfg.reinit_every == 0
                before = interfaceGeometry(field, h)
                gb = gradientStats(phi, h)
                reinitialize!(field, cut(field, :init))
                after = interfaceGeometry(field, h)
                ga = gradientStats(phi, h)
                push!(reinit_shifts, (radius=abs(after.radius - before.radius),
                                      volume=abs(after.volume - before.volume)))
                push!(reinit_grad_before, gb)
                push!(reinit_grad_after, ga)
            end

            clamped += computeNormals!(workspace.normals, phi, h, cfg.eps_gradient)
            nonfinite += computeCurvature!(workspace.curvature, workspace.normals, h)
        end
        push!(step_times, steptime)

        if next_target <= length(cfg.coordinates) && t >= cfg.coordinates[next_target] - 1e-12
            tag = next_target == length(cfg.coordinates) ? "final" :
                  @sprintf("%.3f", cfg.coordinates[next_target])
            row = diagnostics(field, workspace, cfg, t)
            saveSnapshot(outputdir, tag, row, phi, cfg.n)
            push!(rows, row)
            push!(tags, tag)
            next_target += 1
        end
    end
    elapsed = time() - tstart

    writeTimeseries(joinpath(outputdir, "timeseries.csv"), rows)
    open(joinpath(outputdir, "snapshots.csv"), "w") do io
        println(io, "tag,t,n")
        for (i, tag) in enumerate(tags)
            println(io, tag, ",", _fmt(rows[i].t), ",", cfg.n)
        end
    end
    summary = (cfg=cfg, h=h, dt=dt, nsteps=step, elapsed=elapsed, rows=rows,
               final=rows[end], reinit_shifts=reinit_shifts,
               reinit_grad_before=reinit_grad_before, reinit_grad_after=reinit_grad_after,
               clamped=clamped, nonfinite=nonfinite, tags=tags,
               step_times=step_times, outputdir=outputdir)
    writeSummary(summary)
    return summary
end

function writeSummary(s)
    cfg, final = s.cfg, s.final
    path = joinpath(s.outputdir, "summary.md")
    open(path, "w") do io
        println(io, "# MeanCurvatureSphere3D — ", cfg.n, "^3, dt = ", round(s.dt; sigdigits=5),
                ", ", s.nsteps, " steps")
        println(io)
        println(io, "R0 = ", cfg.R0, ", mu = ", cfg.mu, ", t_final = ", final.t,
                ", C = ", cfg.dt_coefficient, ", reinit every ", cfg.reinit_every, " steps, ",
                "Curvature.jl window ", cfg.curvature_radius, " h")
        println(io)
        println(io, "| quantity | numerical | exact | error | relative |")
        println(io, "|---|---|---|---|---|")
        @printf(io, "| radius | %.10f | %.10f | %+.3e | %+.3e |\n",
                final.radius, final.exact_radius, final.radius_err, final.radius_rel_err)
        @printf(io, "| volume | %.10f | %.10f | %+.3e | %+.3e |\n",
                final.volume, final.exact_volume, final.volume - final.exact_volume,
                final.volume_rel_err)
        @printf(io, "| mean curvature | %.6f | %.6f | %+.3e | %+.3e |\n",
                final.kappa_cv, final.kappa_exact, final.kappa_cv - final.kappa_exact,
                (final.kappa_cv - final.kappa_exact) / final.kappa_exact)
        println(io)
        println(io, "Sphericity: area-weighted std of the interface radius = ",
                round(final.radius_std; sigdigits=4), " (",
                round(100 * final.radius_std / final.radius; sigdigits=3), " % of R).")
        println(io, "Axis crossings at the final time: ",
                join(round.(final.axis_crossings; digits=6), ", "),
                " (max |error| = ", round(final.axis_err_max; sigdigits=4), ").")
        println(io, "Curvature: FD band mean ", round(final.kappa_fd; sigdigits=6),
                " (std ", round(final.kappa_fd_std; sigdigits=3), "), Curvature.jl area mean ",
                round(final.kappa_cv; sigdigits=6), " (std ",
                round(final.kappa_cv_std; sigdigits=3), "), 2/R = ",
                round(final.kappa_exact; sigdigits=6), ".")
        println(io, "FD vs Curvature.jl at facet centroids: area-weighted rms ",
                round(100 * final.kappa_fd_vs_cv_rms; sigdigits=3), " % of 2/R.")
        println(io)
        println(io, "Non-finite curvature values: ", s.nonfinite,
                "; gradient-clamped node evaluations (per-step counts summed): ", s.clamped)
        if !isempty(s.reinit_shifts)
            println(io, "\nReinitialisation (", length(s.reinit_shifts), " calls): max |change| ",
                    "of the interface radius = ",
                    round(maximum(x.radius for x in s.reinit_shifts); sigdigits=3),
                    ", of the enclosed volume = ",
                    round(maximum(x.volume for x in s.reinit_shifts); sigdigits=3),
                    "; mean |grad phi| in the band ",
                    round(sum(s.reinit_grad_before) / length(s.reinit_grad_before); sigdigits=5),
                    " -> ",
                    round(sum(s.reinit_grad_after) / length(s.reinit_grad_after); sigdigits=5),
                    ".")
        end
        println(io, @sprintf("\nRuntime %.2f s (incl. snapshots and figures); %.2f ms per transport step; %.1f M node-visits/s.",
                             s.elapsed, 1000 * sum(s.step_times) / s.nsteps,
                             s.nsteps * cfg.n^3 / sum(s.step_times) / 1e6))
    end
    return path
end

function printSummary(s)
    cfg, final = s.cfg, s.final
    println("\n=== MeanCurvatureSphere3D, ", cfg.n, "^3 = ", cfg.n^3, " nodes, h = ",
            round(s.h; sigdigits=5), " ===")
    println("dt = ", round(s.dt; sigdigits=6), " (C = ", cfg.dt_coefficient,
            ", mu*dt/h^2 = ", round(cfg.mu * s.dt / s.h^2; sigdigits=4), "), ", s.nsteps,
            " steps, t = ", round(final.t; digits=4))
    println("radius      : ", round(final.radius; digits=8), "  exact ",
            round(final.exact_radius; digits=8), "  err ", round(final.radius_err; sigdigits=4),
            "  rel ", round(final.radius_rel_err; sigdigits=4))
    println("volume      : ", round(final.volume; digits=8), "  exact ",
            round(final.exact_volume; digits=8), "  rel ", round(final.volume_rel_err; sigdigits=4))
    println("sphericity  : std(R) = ", round(final.radius_std; sigdigits=4), "  axis err max = ",
            round(final.axis_err_max; sigdigits=4))
    println("curvature   : Curvature.jl ", round(final.kappa_cv; digits=6), "  FD band ",
            round(final.kappa_fd; digits=6), "  exact 2/R = ", round(final.kappa_exact; digits=6))
    println("curvature   : FD vs Curvature.jl rms ",
            round(100 * final.kappa_fd_vs_cv_rms; sigdigits=3), " % of 2/R")
    println("integrity   : non-finite curvature ", s.nonfinite, ", reinit calls ",
            length(s.reinit_shifts), ", max reinit radius shift ",
            isempty(s.reinit_shifts) ? 0.0 : round(maximum(x.radius for x in s.reinit_shifts); sigdigits=3))
    @printf("performance : %.2f s total, %.1f ms/step, %.1f M node-visits/s in the time loop\n",
            s.elapsed, 1000 * sum(s.step_times) / s.nsteps,
            s.nsteps * cfg.n^3 / sum(s.step_times) / 1e6)
    println("snapshots   : ", join(s.tags, ", "), " in ", s.outputdir)
end

function main(args=ARGS)
    n = isempty(args) ? 48 : parse(Int, args[1])
    cfg = Config(; n=n)
    ordering_error = verifyTransportOrdering()
    @printf("transport ordering self-check: |error| = %.2e %s\n", ordering_error,
            ordering_error < 1e-10 ? "(ok)" : "(FAILED)")
    ordering_error < 1e-10 || error("velocity/level-set node ordering is inconsistent")
    summary = simulate(cfg)
    printSummary(summary)

    # The entry point saves every figure.
    visualizer = joinpath(@__DIR__, "visualize.MeanCurvatureSphere3D.jl")
    isfile(visualizer) || error("visualizer not found: $visualizer")
    include(visualizer)
    # `include` defines the figure module at run time; call it in the current world age.
    Base.invokelatest(MeanCurvatureSphereFigures.makeFigures, summary.outputdir)
    return summary
end

if abspath(PROGRAM_FILE) == @__FILE__
    main()
end
