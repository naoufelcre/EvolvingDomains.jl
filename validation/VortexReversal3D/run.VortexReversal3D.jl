#!/usr/bin/env julia
# =============================================================================
# VortexReversal3D — three-dimensional reversible vortex benchmark.
#
#   D = [0,1]^3,  0 ≤ t ≤ T = 1
#   v(x,t) = cos(πt/T) V(x)          (divergence free, tangent to ∂D, reverses at T/2)
#   ∂_t φ + v·∇φ = 0                 φ_0 = |x - (½,½,½)| - 0.16
#   ∂_t c + v·∇c = 0                 c_0 = exp(-|x - (0.43,0.52,0.49)|² / 0.08²)
#
# Because v(·,T-t) = -v(·,t), the exact fields at t = T equal the initial ones, so every
# error below is a *return* error and needs no manufactured solution.
#
# EvolvingDomains is used as a data representation (Grid / Discretization / Field) plus two
# transport kernels (eulerianTransport! = WENO5 + SSP-RK3, lagrangianTransport! = CIP).
# Everything else — the velocity, the time loop, the diagnostics — is plain array code.
#
# Run:  julia --project=validation/env validation/VortexReversal3D/run.VortexReversal3D.jl
# =============================================================================

# The native CutCells artifact is not published in this checkout (no Artifacts.toml), so
# `cut` needs the documented environment override. Point it at the library built in this
# repository unless the caller already chose one. Never touch package sources.
let stem = (Sys.iswindows() ? "" : "lib") * "evolvingdomains_cut_cells." *
           (Sys.iswindows() ? "dll" : Sys.isapple() ? "dylib" : "so"),
    local_lib = normpath(joinpath(@__DIR__, "..", "lib", "cutcells", "lib", stem))
    isempty(get(ENV, "EVOLVINGDOMAINS_CUTCELLS_LIBRARY", "")) && isfile(local_lib) &&
        (ENV["EVOLVINGDOMAINS_CUTCELLS_LIBRARY"] = local_lib)
end

using EvolvingDomains
using StaticArrays: SVector
using LinearAlgebra: norm, dot, cross
using Printf

# ----------------------------------------------------------------------------- model
const T_END = 1.0
const CENTRE = SVector(0.5, 0.5, 0.5)
const RADIUS = 0.16
const BUMP = SVector(0.43, 0.52, 0.49)
const BUMP_WIDTH = 0.08
const EXACT_VOLUME = 4π * RADIUS^3 / 3

"""Spatial factor V of the velocity; v(x,t) = cos(πt/T) V(x). Symmetric in the sense that
each component is a product of one squared sine and two full sines, which is what makes
∇·V = 0 and V·n = 0 on ∂D."""
prescribedVelocity(p) =
    SVector(2sinpi(p[1])^2 * sinpi(2p[2]) * sinpi(2p[3]),
            -sinpi(2p[1]) * sinpi(p[2])^2 * sinpi(2p[3]),
            -sinpi(2p[1]) * sinpi(2p[2]) * sinpi(p[3])^2)

velocity(p, t) = cos(π * t / T_END) * prescribedVelocity(p)

# ----------------------------------------------------------------------- configuration
"""
    configuration(nodes; cfl = 0.4)

Everything the run needs, as arrays and named tuples: grid, fields, current time, step and
numerical parameters. `dt` comes from a CFL bound on the *frozen* velocity magnitude
(`|cos| ≤ 1`, so t = 0 is the worst instant) and is then rounded down so that `steps` is
even and `t = T/2` is a step boundary.
"""
function configuration(nodes::Int; cfl::Float64 = 0.4)
    axes = ntuple(_ -> collect(range(0.0, 1.0; length = nodes)), 3)
    grid = Grid(axes)
    discretization = Discretization(Dict(:init => grid))
    # node order, first axis fastest: `vec` of an (nx,ny,nz) array (README, "Velocity")
    points = vec([SVector(x, y, z) for x in axes[1], y in axes[2], z in axes[3]])
    phi0 = [norm(SVector(x, y, z) - CENTRE) - RADIUS
            for x in axes[1], y in axes[2], z in axes[3]]
    c0 = [exp(-norm(SVector(x, y, z) - BUMP)^2 / BUMP_WIDTH^2)
          for x in axes[1], y in axes[2], z in axes[3]]
    field(values) = Field(discretization, Dict(:init => copy(values)))

    h = minimum(minimum(diff(axis)) for axis in axes)
    vmax = maximum(norm ∘ prescribedVelocity, points)
    steps = 2 * ceil(Int, T_END / (2 * (cfl * h / vmax)))
    dt = T_END / steps

    return (; nodes, axes, grid, discretization, points, h, cfl, vmax, dt, steps,
              t_end = T_END, phi0, c0,
              geometry = field(phi0),
              scalarEuler = field(c0),        # passive scalar, WENO5 + SSP-RK3
              scalarCip = field(c0),          # passive scalar, CIP (source buffer)
              scalarCipTarget = field(c0),    # CIP destination; must not alias the source
              velocity = Vector{SVector{3,Float64}}(undef, length(points)),
              time = Ref(0.0),
              weno = WENO5Cache(),
              cip = CIPCache(vec(copy(c0)), grid))
end

"""
    verifyVelocityOrdering(config)

Self-check of the untyped, load-bearing node ordering. A CIP step on an affine field under
a constant velocity is exact in the interior only if `config.velocity[k]` belongs to
`config.points[k]`; a transposed ordering gives O(1) errors instead.
"""
function verifyVelocityOrdering(config)
    affine = p -> 1.0 + 2.0p[1] - 0.5p[2] + 3.0p[3]
    shift, dt = SVector(0.01, -0.008, 0.006), 0.2
    n = config.nodes
    source = Field(config.discretization,
                   Dict(:init => reshape([affine(p) for p in config.points], n, n, n)))
    target = Field(config.discretization, Dict(:init => zeros(n, n, n)))
    v = [all(0 < c < 1 for c in p) ? shift : SVector(0.0, 0.0, 0.0) for p in config.points]
    lagrangianTransport!(target, source, v, dt; level_id = :init)
    exact = reshape([affine(p - dt * shift) for p in config.points], n, n, n)
    r = 4:n-3
    err = maximum(abs.(target.data[:init][r, r, r] .- exact[r, r, r]))
    @printf("velocity ordering self-check (%d³): max CIP error on an affine field %.2e\n",
            n, err)
    err < 1e-9 || error("velocity node ordering does not match the kernels (error $err)")
    return err
end

# ----------------------------------------------------------------------------- dynamics
"""
    dynamics(config, t)

Prescribed nodal velocity at time `t`, in the package ordering. The buffer is reused: a
fresh `Vector{SVector{3,Float64}}` per step would be 6 MB at 64³.
"""
function dynamics(config, t)
    v = config.velocity
    @inbounds for k in eachindex(config.points)
        v[k] = velocity(config.points[k], t)
    end
    return v
end

# ------------------------------------------------------------------------------- steps
"""Advance the level set with WENO5 + SSP-RK3 on the frozen velocity `v`.

No reinitialization: the benchmark measures the transport error of a reversing flow, so
`reinitialize!` is deliberately not called and no φ is ever reset from the initial sphere.
"""
function advanceGeometry!(config, v)
    eulerianTransport!(config.geometry, v, config.dt; level_id = :init, cache = config.weno)
    return config
end

"""
    advanceState!(config, scheme, v)

Advance the passive scalar with `:weno5` (eulerianTransport!) or `:cip`
(lagrangianTransport!). The CIP call keeps source and destination buffers distinct and
copies the result back, since the operator is documented as non-aliasing.
"""
function advanceState!(config, scheme::Symbol, v)
    if scheme === :weno5
        eulerianTransport!(config.scalarEuler, v, config.dt; level_id = :init,
                           cache = config.weno)
    elseif scheme === :cip
        lagrangianTransport!(config.scalarCipTarget, config.scalarCip, v, config.dt;
                             level_id = :init, cache = config.cip)
        copyto!(config.scalarCip.data[:init], config.scalarCipTarget.data[:init])
    else
        throw(ArgumentError("unknown scalar transport scheme $(repr(scheme))"))
    end
    return config
end

snapshot(config, t) = (; t,
                       phi = copy(config.geometry.data[:init]),
                       c_weno5 = copy(config.scalarEuler.data[:init]),
                       c_cip = copy(config.scalarCip.data[:init]))

# ------------------------------------------------------------------------------ driver
"""
    simulate(config)

The explicit time loop: velocity evaluated once per global step at the midpoint
`t + dt/2` and shared by the level set and both scalars, then `t → t + dt`. Midpoint
freezing is second order in `dt` for a smooth `v`; the spatial orders are WENO5/CIP for
the scalars and WENO5 for φ, all limited by SSP-RK3 and the two-level time differencing.
"""
function simulate(config)
    config.time[] = 0.0
    snaps = [snapshot(config, 0.0)]
    half = config.steps ÷ 2

    # compilation, cache warm-up and first allocation, excluded from the reported wall clock
    warm = (level = Field(config.discretization, Dict(:init => copy(config.phi0))),
            scalar = Field(config.discretization, Dict(:init => copy(config.c0))),
            target = Field(config.discretization, Dict(:init => zeros(config.nodes,
                                                                     config.nodes,
                                                                     config.nodes))))
    v0 = dynamics(config, 0.0)
    eulerianTransport!(warm.level, v0, config.dt; level_id = :init, cache = config.weno)
    eulerianTransport!(warm.scalar, v0, config.dt; level_id = :init, cache = config.weno)
    lagrangianTransport!(warm.target, warm.scalar, v0, config.dt; level_id = :init,
                         cache = config.cip)

    stats = @timed begin
        for step in 1:config.steps
            v = dynamics(config, config.time[] + config.dt / 2)
            advanceGeometry!(config, v)
            advanceState!(config, :weno5, v)
            advanceState!(config, :cip, v)
            config.time[] = step * config.dt   # recomputed, not accumulated
            (step == half || step == config.steps) &&
                push!(snaps, snapshot(config, config.time[]))
        end
        nothing
    end
    return (; snapshots = snaps, seconds = stats.time, bytes = stats.bytes,
              gctime = stats.gctime)
end

# -------------------------------------------------------------------------- diagnostics
relativeErrors(got, exact) = (l1 = sum(abs, got .- exact) / sum(abs, exact),
                              l2 = sqrt(sum(abs2, got .- exact)) / sqrt(sum(abs2, exact)))

# Grid-based estimates: nodal sums, h³ per node. They integrate the *node* samples, not the
# cut cells, so they must not be read as cut-cell volumes.
voxelVolume(phi, h) = h^3 * count(<(0), phi)
"""Heaviside smeared over one cell; exact for a planar interface, O(h²) for a curved one."""
smearedVolume(phi, h) = h^3 * sum(clamp(0.5 - φ / h, 0.0, 1.0) for φ in phi)

"""
    cutCellVolume(config, field)

Cut-cell-integrated volume of {φ < 0} by the divergence theorem over the interface facets
of the package's own `cut` snapshot: `V = (1/3) ∮ x·n dA`. Valid as an enclosed volume only
when the negative phase does not touch ∂D (`wallContacts == 0`), otherwise the box walls
would contribute a surface the snapshot does not describe.
"""
function cutCellVolume(config, field)
    c = try
        cut(field, :init)
    catch err
        err isa ArgumentError || rethrow()      # only the documented "no artifact" path
        return (; volume = missing, wall_contacts = missing, note = sprint(showerror, err))
    end
    s = c.boundary
    volume = 0.0
    for f in eachindex(s.facet_to_bgcell)
        base = s.point_pointers[f]
        a = SVector(s.physical_points[s.point_ids[base]])
        b = SVector(s.physical_points[s.point_ids[base + 1]])
        e = SVector(s.physical_points[s.point_ids[base + 2]])
        # (1/3)·centroid·n·area with centroid = (a+b+e)/3 and area = |(b-a)×(e-a)|/2
        volume += dot(a + b + e, SVector(s.normals[f])) * norm(cross(b - a, e - a)) / 18
    end
    return (; volume, wall_contacts = wallContacts(c, config.grid), note = "")
end

"""Number of background cells on ∂D that are not entirely outside the negative phase."""
function wallContacts(c::Cut, grid)
    dims = ntuple(d -> length(grid.axes[d]) - 1, 3)
    touched = 0
    for k in 1:dims[3], j in 1:dims[2], i in 1:dims[1]
        (i == 1 || j == 1 || k == 1 || i == dims[1] || j == dims[2] || k == dims[3]) ||
            continue
        c.bg_states[i + (j - 1) * dims[1] + (k - 1) * dims[1] * dims[2]] == CUT_OUT ||
            (touched += 1)
    end
    return touched
end

"""    diagnostics(config, run) -> NamedTuple

Every number is computed from the saved arrays; nothing is hardcoded except the analytic
reference volume 4πR³/3 and the initial data the run started from.
"""
function diagnostics(config, run)
    snaps = run.snapshots
    rows = map(snaps) do s
        # `cut` snapshots a Field, so the saved array is re-wrapped (no copy: retained)
        cell = cutCellVolume(config, Field(config.discretization, Dict(:init => s.phi)))
        isempty(cell.note) || @warn "cut-cell volume unavailable" cell.note maxlog = 1
        (; t = s.t, voxel = voxelVolume(s.phi, config.h),
           smeared = smearedVolume(s.phi, config.h),
           cut = cell.volume, wall_contacts = cell.wall_contacts,
           phi_min = minimum(s.phi), phi_max = maximum(s.phi),
           weno5_min = minimum(s.c_weno5), weno5_max = maximum(s.c_weno5),
           cip_min = minimum(s.c_cip), cip_max = maximum(s.c_cip))
    end
    errors = (; phi = relativeErrors(snaps[end].phi, config.phi0),
                weno5 = relativeErrors(snaps[end].c_weno5, config.c0),
                cip = relativeErrors(snaps[end].c_cip, config.c0))
    finite = all(all(isfinite, a) for s in snaps for a in (s.phi, s.c_weno5, s.c_cip))
    return (; nodes = config.nodes, dt = config.dt, steps = config.steps, cfl = config.cfl,
              h = config.h, vmax = config.vmax, t_end = config.t_end,
              seconds = run.seconds, bytes = run.bytes, gctime = run.gctime,
              exact_volume = EXACT_VOLUME, finite, rows, errors)
end

# ------------------------------------------------------------------------ persistence
"""
    writeRun(dir, config, run, m)

`params.txt` (what the reader needs to reshape `fields.bin`), `metrics.csv` (one row per
snapshot) and `errors.csv`. `fields.bin` holds `Float64` records of `nodes³` in
column-major order (first axis fastest), in this order:

    1  phi0              2  c0
    3  phi       @ t0     4  phi       @ T/2    5  phi       @ T
    6  c_weno5   @ t0     7  c_weno5   @ T/2    8  c_weno5   @ T
    9  c_cip     @ t0    10  c_cip     @ T/2   11  c_cip     @ T
"""
function writeRun(dir, config, run, m)
    mkpath(dir)
    open(joinpath(dir, "params.txt"), "w") do io
        for (k, v) in (("nodes", config.nodes), ("steps", config.steps),
                       ("dt", config.dt), ("cfl", config.cfl), ("h", config.h),
                       ("vmax", config.vmax), ("t_end", config.t_end),
                       ("box", join((first(config.axes[1]), last(config.axes[1])), " ")))
            println(io, k, " ", v)
        end
    end
    open(joinpath(dir, "fields.bin"), "w") do io
        write(io, config.phi0)
        write(io, config.c0)
        for get in (s -> s.phi, s -> s.c_weno5, s -> s.c_cip), s in run.snapshots
            write(io, get(s))
        end
    end
    open(joinpath(dir, "metrics.csv"), "w") do io
        println(io, "time,voxel,smeared,cutcell,wall_contacts,phi_min,phi_max,",
                    "weno5_min,weno5_max,cip_min,cip_max")
        for r in m.rows
            println(io, join((r.t, r.voxel, r.smeared, r.cut, r.wall_contacts,
                              r.phi_min, r.phi_max, r.weno5_min, r.weno5_max,
                              r.cip_min, r.cip_max), ','))
        end
    end
    open(joinpath(dir, "errors.csv"), "w") do io
        println(io, "quantity,l1,l2")
        for (name, e) in (("phi", m.errors.phi), ("c_weno5", m.errors.weno5),
                          ("c_cip", m.errors.cip))
            println(io, join((name, e.l1, e.l2), ','))
        end
    end
    return dir
end

# ------------------------------------------------------------------------------ report
function report(io::IO, m)
    @printf(io, "\n=== VortexReversal3D — %d³ nodes, %d steps, dt = %.6g (CFL %.2f) ===\n",
            m.nodes, m.steps, m.dt, m.cfl)
    @printf(io, "wall clock %.2f s (GC %.2f s), %.3g GB allocated, all values finite: %s\n",
            m.seconds, m.gctime, m.bytes / 1e9, m.finite)
    @printf(io, "%-10s %7s %12s %12s %12s %12s\n", "t", "", "V_voxel", "V_smeared",
            "V_cutcell", "wall cells")
    for r in m.rows
        @printf(io, "%-10.3f %7s %12.6g %12.6g %12s %12s\n", r.t, "", r.voxel, r.smeared,
                r.cut === missing ? "n/a" : @sprintf("%.6g", r.cut),
                r.wall_contacts === missing ? "n/a" : string(r.wall_contacts))
    end
    @printf(io, "exact volume 4πR³/3 = %.9g\n", m.exact_volume)
    for r in m.rows
        @printf(io, "t=%.2f  grid volume drift %+.3e / %+.3e (voxel/smeared), cut-cell %s\n",
                r.t, (r.voxel - m.rows[1].voxel) / m.rows[1].voxel,
                (r.smeared - m.rows[1].smeared) / m.rows[1].smeared,
                r.cut === missing ? "n/a" :
                    @sprintf("%+.3e", (r.cut - m.rows[1].cut) / m.rows[1].cut))
    end
    for (name, e) in (("phi    ", m.errors.phi), ("c/weno5", m.errors.weno5),
                      ("c/cip  ", m.errors.cip))
        @printf(io, "return error %s  L1 %.4e   L2 %.4e\n", name, e.l1, e.l2)
    end
    for r in m.rows
        @printf(io, "t=%.2f  phi in [%.4f, %.4f]  c_weno5 in [%9.2e, %.4f]  c_cip in [%9.2e, %.4f]\n",
                r.t, r.phi_min, r.phi_max, r.weno5_min, r.weno5_max, r.cip_min, r.cip_max)
    end
    return io
end

# ------------------------------------------------------------------------------- entry
function main()
    resolutions = isempty(ARGS) ? (40, 64) : (parse(Int, ARGS[1]),)
    outdir = joinpath(@__DIR__, "output")
    mkpath(outdir)
    dirs = String[]

    walltime = @elapsed for nodes in resolutions
        config = configuration(nodes)
        verifyVelocityOrdering(config)
        run = simulate(config)
        m = diagnostics(config, run)
        report(stdout, m)
        dir = joinpath(outdir, "n$nodes")
        writeRun(dir, config, run, m)
        open(joinpath(dir, "metrics.txt"), "w") do io
            report(io, m)
        end
        push!(dirs, dir)
    end
    @printf("\ntotal simulation wall clock: %.2f s\n", walltime)

    visualizeAll(dirs, outdir)
    return dirs
end

include(joinpath(@__DIR__, "visualize.VortexReversal3D.jl"))
main()
