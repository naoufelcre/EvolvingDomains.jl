#!/usr/bin/env julia
# =============================================================================
# StrongDeformation3D — three-dimensional reversible transport with a strongly
# deforming initial geometry, plus an isolated allocation study of the kernels.
#
#   D = [0,1]^3,  0 ≤ t ≤ T = 1
#   v(x,t) = A cos(πt/T) V(x)        (divergence free, tangent to ∂D, reverses at T/2)
#   ∂_t φ + v·∇φ = 0                 φ_0 = a_min (sqrt(Σ((x-c_i)/a_i)²) - 1)
#   ∂_t c + v·∇c = 0                 c_0 = exp(-|x - (0.42,0.48,0.56)|² / 0.10²)
#
# Because v(·,T-t) = -v(·,t) the exact fields at t = T equal the initial ones, so every
# error below is a *return* error and needs no manufactured solution. No reinitialization
# and no analytic geometry ever overwrites the numerical state.
#
# EvolvingDomains is used as a data representation (Grid / Discretization / Field), two
# transport kernels (eulerianTransport! = WENO5 + SSP-RK3, lagrangianTransport! = CIP) and
# the `cut` snapshot for geometry. Everything else — velocity, time loop, diagnostics,
# allocation study, persistence — is plain array code here.
#
# Run:  julia --project=validation/env \
#         validation/VortexReversal3D/round2/run.StrongDeformation3D.jl [48] [64]
# =============================================================================

# The native CutCells artifact is not published in this checkout (no Artifacts.toml), so
# `cut` needs the documented environment override pointing at the library built in this
# repository. Never touch package sources.
let stem = "libevolvingdomains_cut_cells." * (Sys.isapple() ? "dylib" : Sys.islinux() ? "so" : "dll"),
    local_lib = normpath(joinpath(@__DIR__, "..", "..", "lib", "cutcells", "lib", stem))
    isempty(get(ENV, "EVOLVINGDOMAINS_CUTCELLS_LIBRARY", "")) && isfile(local_lib) &&
        (ENV["EVOLVINGDOMAINS_CUTCELLS_LIBRARY"] = local_lib)
end

using EvolvingDomains
using StaticArrays: SVector
using StaticArrays
using LinearAlgebra: norm, dot, cross
using Printf
using Profile

# ----------------------------------------------------------------------------- model
const T_END = 1.0
const AMP = 2.0                       # primary deformation strength
const CENTRE = SVector(0.41, 0.48, 0.56)
const SEMI = SVector(0.14, 0.11, 0.09)
const AMIN = 0.09
const BUMP = SVector(0.42, 0.48, 0.56)
const BUMP_WIDTH = 0.10
const EXACT_VOLUME = 4π / 3 * prod(SEMI)
const SNAP_TIMES = (0.0, 0.25, 0.5, 0.75, 1.0)

"""Spatial factor V of the velocity. Divergence free and tangent to ∂D by construction:
every component is a product of one squared and two full sines, so each ∂_i V_i contains
sin(2πx_i) and cancels against the two other divergence terms."""
prescribedVelocity(p) =
    SVector(2sinpi(p[1])^2 * sinpi(2p[2]) * sinpi(2p[3]),
            -sinpi(2p[1]) * sinpi(p[2])^2 * sinpi(2p[3]),
            -sinpi(2p[1]) * sinpi(2p[2]) * sinpi(p[3])^2)

velocity(p, t) = AMP * cos(π * t / T_END) * prescribedVelocity(p)

"""Zero level set is exactly the ellipsoid {Σ((x-c)/a)² = 1}, negative inside. It is *not* a
signed-distance function: |∇φ_0| ≠ 1 in general, so |φ_0(p)| is only a first-order proxy for
the distance from p to the interface. All distances below are normalised by `AMIN`."""
phi0At(p) = AMIN * (sqrt(sum(abs2, (p - CENTRE) ./ SEMI)) - 1)

c0At(p) = exp(-sum(abs2, p - BUMP) / BUMP_WIDTH^2)

# ------------------------------------------------------------------------ source check
"""Stage A: the loaded EvolvingDomains must be *this* checkout. Abort otherwise."""
function verifySource(io::IO)
    src = realpath(pathof(EvolvingDomains))
    repo = realpath(joinpath(@__DIR__, "..", "..", ".."))
    pkgdir(EvolvingDomains) == repo || error(
        "EvolvingDomains resolves to $(realpath(pkgdir(EvolvingDomains))); expected $repo. " *
        "Aborting.")
    startswith(src, joinpath(repo, "src")) || error("source path $(src) is not under $repo/src")
    rev = strip(read(`git -C $repo rev-parse HEAD`, String))
    dirty = let s = strip(read(`git -C $repo status --porcelain`, String))
        isempty(s) ? 0 : count(==('\n'), s) + 1
    end
    lib = get(ENV, "EVOLVINGDOMAINS_CUTCELLS_LIBRARY", "")
    @printf(io, "EvolvingDomains source : %s\n", src)
    @printf(io, "  package dir          : %s\n", pkgdir(EvolvingDomains))
    @printf(io, "  version (Project.toml): %s\n", pkgversion(EvolvingDomains))
    @printf(io, "  git revision         : %s (%d modified/untracked paths in the checkout)\n",
            rev, dirty)
    @printf(io, "  CutCells native lib  : %s%s\n", isempty(lib) ? "<UNSET>" : lib,
            isempty(lib) || isfile(lib) ? "" : "  [MISSING]")
    @printf(io, "  julia %s, %d thread(s), StaticArrays %s\n", VERSION, Threads.nthreads(),
            pkgversion(StaticArrays))
    return repo
end

# ---------------------------------------------------------------------- configuration
"""
    configuration(nodes; cfl = 0.4)

Grid, initial data, fields, step size. `dt` comes from a CFL bound on the nodal velocity
maximum, rounded down so that `steps` is a multiple of 4 and therefore t = 0.25, 0.5, 0.75
are all step boundaries.
"""
function configuration(nodes::Int;
                      cfl::Float64 = parse(Float64, get(ENV, "SD3D_CFL", "0.4")))
    axes = ntuple(_ -> collect(range(0.0, 1.0; length = nodes)), 3)
    grid = Grid(axes)
    discretization = Discretization(Dict(:init => grid))
    # node order, first axis fastest: `vec` of an (nx,ny,nz) array (README, "Velocity")
    points = vec([SVector(x, y, z) for x in axes[1], y in axes[2], z in axes[3]])
    phi0 = [phi0At(SVector(x, y, z)) for x in axes[1], y in axes[2], z in axes[3]]
    c0 = [c0At(SVector(x, y, z)) for x in axes[1], y in axes[2], z in axes[3]]
    field(values) = Field(discretization, Dict(:init => copy(values)))

    h = minimum(minimum(diff(axis)) for axis in axes)
    vmax = AMP * maximum(norm ∘ prescribedVelocity, points)
    steps = 4 * ceil(Int, T_END / (4 * (cfl * h / vmax)))
    dt = T_END / steps

    return (; nodes, axes, grid, discretization, points, h, cfl, vmax, dt, steps,
              amplitude = AMP, t_end = T_END, phi0, c0,
              geometry = field(phi0),
              scalarEuler = field(c0),        # passive scalar, WENO5 + SSP-RK3
              scalarCip = field(c0),          # passive scalar, CIP source buffer
              scalarCipTarget = field(c0),    # CIP destination; must not alias the source
              velocity = Vector{SVector{3,Float64}}(undef, length(points)),
              time = Ref(0.0),
              weno = WENO5Cache(),
              cip = CIPCache(vec(copy(c0)), grid))
end

"""Discrete sanity check of the prescribed velocity: ∇·V = 0 inside and V·n = 0 on ∂D."""
function verifyFlow(config, io::IO)
    h = config.h
    div = 0.0
    for p in (SVector(0.21, 0.37, 0.62), SVector(0.5, 0.5, 0.5), SVector(0.83, 0.24, 0.41))
        d = sum((prescribedVelocity(p + h .* SVector(δ...))[i] -
                 prescribedVelocity(p - h .* SVector(δ...))[i]) / (2h)
                for (i, δ) in enumerate(((1, 0, 0), (0, 1, 0), (0, 0, 1))))
        div = max(div, abs(d))
    end
    axis = config.axes[1]
    wall = maximum(abs(prescribedVelocity(SVector(x, y, z))[1])
                   for x in (0.0, 1.0), y in axis, z in axis)
    wall = max(wall, maximum(abs(prescribedVelocity(SVector(x, y, z))[2])
                             for y in (0.0, 1.0), x in axis, z in axis))
    wall = max(wall, maximum(abs(prescribedVelocity(SVector(x, y, z))[3])
                             for z in (0.0, 1.0), x in axis, y in axis))
    @printf(io, "flow self-check (%d³ nodes): max|∇·V| = %.2e, max|V·n| on ∂D = %.2e, max|V| = %.6f, vmax = A·max|V| = %.6f\n",
            config.nodes, div, wall, config.vmax / AMP, config.vmax)
    (div < 1e-10 && wall < 1e-12) ||
        error("prescribed velocity is not divergence free / not tangent to ∂D")
    return nothing
end

"""
    verifyVelocityOrdering(config)

Self-check of the untyped, load-bearing node ordering, for **both** kernels.

The test field is affine and the test velocity is `v = (a, 0, c)·g(y)` with `g(y) = 0.5 + y`.
Because `v_y = 0`, `y` is constant along every trajectory, so the exact solution is known in
closed form and both kernels transport it *exactly* in the interior: the field stays linear,
WENO5 reconstructs a linear function exactly, SSP-RK3 integrates a per-node linear ODE
exactly, and every CIP Hermite profile is a straight line. The velocity varies with `y` only,
so it is not invariant under any axis permutation: applying the same test to a transposed
node order must fail loudly, which is what makes the check meaningful.

`v` is zeroed on the two `x` shells and on `z = 0` so that no CIP departure point leaves the
box; the comparison is made a few nodes inside that discontinuity, where WENO5's clamped
stencil can no longer reach.
"""
function verifyVelocityOrdering(config, io::IO)
    φ(p, t) = 1.0 + 2.0p[1] - 0.5p[2] + 3.0p[3] - t * (0.5 + p[2]) * (2 * 0.01 + 3 * 0.006)
    shift(p) = SVector(0.01 * (0.5 + p[2]), 0.0, 0.006 * (0.5 + p[2]))
    dt, n = 0.2, config.nodes
    source = Field(config.discretization,
                   Dict(:init => reshape([φ(p, 0.0) for p in config.points], n, n, n)))
    exact = reshape([φ(p, dt) for p in config.points], n, n, n)
    v = zeros(SVector{3,Float64}, n, n, n)
    for k in 1:n, j in 1:n, i in 1:n
        (i >= 3 && i <= n - 2 && k >= 2) || continue
        v[i, j, k] = shift(SVector(config.axes[1][i], config.axes[2][j], config.axes[3][k]))
    end
    r = (7:n-6, 5:n-4, 5:n-4)

    function errors(vv)
        f = Field(config.discretization, Dict(:init => copy(source.data[:init])))
        eulerianTransport!(f, vec(vv), dt; level_id = :init, cache = WENO5Cache())
        weno = maximum(abs.(f.data[:init][r...] .- exact[r...]))
        t = Field(config.discretization, Dict(:init => zeros(n, n, n)))
        lagrangianTransport!(t, source, vec(vv), dt; level_id = :init,
                             cache = CIPCache(vec(copy(source.data[:init])), config.grid))
        cip = maximum(abs.(t.data[:init][r...] .- exact[r...]))
        return (; weno, cip)
    end

    good = errors(v)
    # same support, velocity attached to the wrong node: the CIP departure point would leave
    # the box, so the zero shells are restored after permuting
    bad = permutedims(v, (2, 1, 3))
    for k in 1:n, j in 1:n, i in 1:n
        (i >= 3 && i <= n - 2 && k >= 2) || (bad[i, j, k] = SVector(0.0, 0.0, 0.0))
    end
    bad = errors(bad)
    @printf(io, "velocity ordering self-check (%d³, affine field): WENO5 %.2e, CIP %.2e; with a transposed node order WENO5 %.2e, CIP %.2e\n",
            n, good.weno, good.cip, bad.weno, bad.cip)
    (max(good.weno, good.cip) < 1e-7 && min(bad.weno, bad.cip) > 1e-4) ||
        error("velocity node ordering does not match the kernels (correct order $(good), transposed $(bad))")
    return (; good, bad)
end

# ----------------------------------------------------------------------------- dynamics
"""
    dynamics(config, t)

Prescribed nodal velocity at time `t`, in the package ordering. The buffer is reused: a fresh
`Vector{SVector{3,Float64}}` per step would be 6 MB at 64³.
"""
function dynamics(config, t)
    v = config.velocity
    @inbounds for k in eachindex(config.points)
        v[k] = velocity(config.points[k], t)
    end
    return v
end

# ------------------------------------------------------------------------------- steps
"""Advance the level set with WENO5 + SSP-RK3 on the frozen velocity `v`, in place.

No reinitialization is used: the benchmark measures transport plus return error, and
`reinitialize!` would rewrite φ from a reconstructed distance function, destroying exactly
the quantity under test.
"""
function advanceGeometry!(config, v)
    eulerianTransport!(config.geometry, v, config.dt; level_id = :init, cache = config.weno)
    return config
end

"""
    advanceScalar!(config, scheme, v)

Advance the passive scalar with `:weno5` (eulerianTransport!) or `:cip`
(lagrangianTransport!). The CIP operator is non-aliasing, so it writes into a distinct target
buffer which is then copied back into the source: that copy is part of the CIP step's cost
and is measured separately.
"""
function advanceScalar!(config, scheme::Symbol, v)
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

The explicit time loop: velocity evaluated once per step at the midpoint `t + dt/2` and
shared by the level set and both scalars, then `t → t + dt`. Midpoint freezing is second
order in `dt` and exactly time-symmetric (cos(π(1-t)) = -cos(πt)), so it does not break the
reversibility of the benchmark.
"""
function simulate(config)
    config.time[] = 0.0
    snaps = [snapshot(config, 0.0)]
    q, half, three_q = config.steps ÷ 4, config.steps ÷ 2, 3 * config.steps ÷ 4

    # compilation, cache warm-up and first allocation, excluded from the reported wall clock
    warm = (level = Field(config.discretization, Dict(:init => copy(config.phi0))),
            scalar = Field(config.discretization, Dict(:init => copy(config.c0))),
            target = Field(config.discretization, Dict(:init => copy(config.c0))))
    v0 = dynamics(config, 0.0)
    eulerianTransport!(warm.level, v0, config.dt; level_id = :init, cache = config.weno)
    eulerianTransport!(warm.scalar, v0, config.dt; level_id = :init, cache = config.weno)
    lagrangianTransport!(warm.target, warm.scalar, v0, config.dt; level_id = :init,
                         cache = config.cip)

    stats = @timed begin
        for step in 1:config.steps
            v = dynamics(config, (step - 0.5) * config.dt)
            advanceGeometry!(config, v)
            advanceScalar!(config, :weno5, v)
            advanceScalar!(config, :cip, v)
            config.time[] = step * config.dt           # recomputed, not accumulated
            (step == q || step == half || step == three_q || step == config.steps) &&
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

"""Symmetric difference of two enclosed sets, both estimated on the same nodal grid, plus the
overlap measures that separate genuine deformation from spurious volume loss: a shrinking
`|Ω|` inflates the symmetric difference on its own, so `jaccard` (= 1 − D_sym/2 for equal
volumes) and `overlap` (share of the initial set still covered) are reported alongside it."""
function symmetricDifference(phiA, phiB, h)
    HA(k) = clamp(0.5 - k / h, 0.0, 1.0)                     # ≈ 1 inside (φ < 0)
    sa, sb = sum(HA, phiA), sum(HA, phiB)
    sym = sum(abs(HA(a) - HA(b)) for (a, b) in zip(phiA, phiB))
    both = sum(HA(a) * HA(b) for (a, b) in zip(phiA, phiB))
    either = sa + sb - both
    voxel = count(((a < 0) ⊻ (b < 0)) for (a, b) in zip(phiA, phiB))
    return (; smeared = sym / sa, voxel = voxel / count(<(0), phiA),
              jaccard = both / either, overlap = both / sa)
end

"""Trilinear sample of a nodal array on the uniform [0,1]³ node lattice."""
function trilinear(phi::Array{Float64,3}, p, nodes::Int)
    u = clamp.(p .* (nodes - 1) .+ 1, 1.0, Float64(nodes))
    i0 = min.(floor.(Int, u), nodes - 1)
    w = u .- i0
    (i, j, k) = i0
    @inbounds sum(phi[i + a, j + b, k + c] *
                  (a == 0 ? 1 - w[1] : w[1]) *
                  (b == 0 ? 1 - w[2] : w[2]) *
                  (c == 0 ? 1 - w[3] : w[3])
                  for a in 0:1, b in 0:1, c in 0:1)
end

"""Fibonnaci-sphere sample of the *exact* initial ellipsoid surface: the reference point set
for the backward interface distance."""
ellipsoidSurface(m::Int) = [let z = 1 - 2(i - 0.5) / m, r = sqrt(1 - z^2),
                                θ = i * π * (3 - sqrt(5))
        CENTRE + SVector(SEMI[1] * r * cos(θ), SEMI[2] * r * sin(θ), SEMI[3] * z)
    end for i in 1:m]

"""Facet-integrated volume of {φ < 0} from the package's own `cut` snapshot, by the
divergence theorem over the interface facets: V = (1/3) ∮ x·n dA. Valid as an enclosed
volume only when the negative phase does not touch ∂D (`wall_contacts == 0`)."""
function cutGeometry(config, phi)
    c = try
        cut(Field(config.discretization, Dict(:init => phi)), :init)
    catch err
        err isa ArgumentError || rethrow()          # only the documented "no artifact" path
        return (; volume = missing, wall_contacts = missing, points = NTuple{3,Float64}[],
                  note = sprint(showerror, err))
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
    return (; volume, wall_contacts = wallContacts(c, config), points = s.physical_points,
              note = "")
end

"""Number of background cells on ∂D that are not entirely outside the negative phase."""
function wallContacts(c::Cut, config)
    dims = ntuple(d -> length(config.axes[d]) - 1, 3)
    touched = 0
    for k in 1:dims[3], j in 1:dims[2], i in 1:dims[1]
        (i == 1 || j == 1 || k == 1 || i == dims[1] || j == dims[2] || k == dims[3]) ||
            continue
        c.bg_states[i + (j - 1) * dims[1] + (k - 1) * dims[1] * dims[2]] == CUT_OUT ||
            (touched += 1)
    end
    return touched
end

"""
    interfaceDistance(config, phi, reference_points)

Two-sided, sampled interface distance, normalised by `AMIN`.

* forward: every reconstructed interface vertex `p` of `phi`, evaluated in the **exact**
  algebraic φ₀ (`phi0At`), so the reference carries no interpolation error;
* backward: every point of the exact initial ellipsoid surface, where `phi` is sampled by
  trilinear interpolation (error O(h²)).

Both directions are sampled point sets, not a maximisation over the whole surface, and
|φ₀(p)| is only a first-order surrogate for the Euclidean distance because φ₀ is not a
signed-distance function. Hence this is a *sampled, two-sided Hausdorff-type* measure, not an
exact Hausdorff distance.
"""
function interfaceDistance(config, phi, cutpoints, reference_points)
    fwd = [abs(phi0At(SVector(p))) / AMIN for p in cutpoints]
    bwd = [abs(trilinear(phi, p, config.nodes)) / AMIN for p in reference_points]
    return (; mean_fwd = isempty(fwd) ? missing : sum(fwd) / length(fwd),
              max_fwd = isempty(fwd) ? missing : maximum(fwd),
              mean_bwd = sum(bwd) / length(bwd), max_bwd = maximum(bwd),
              n_fwd = length(fwd), n_bwd = length(bwd))
end

"""    diagnostics(config, run) -> NamedTuple

Every number is computed from the saved arrays; nothing is hardcoded except the analytic
reference volume 4πabc/3 and the initial data the run started from.
"""
function diagnostics(config, run)
    snaps = run.snapshots
    reference = ellipsoidSurface(4000)
    geom = map(snaps) do s
        g = cutGeometry(config, s.phi)
        isempty(g.note) || @warn "cut-cell geometry unavailable" g.note maxlog = 1
        sym = symmetricDifference(config.phi0, s.phi, config.h)
        iface = interfaceDistance(config, s.phi, g.points, reference)
        (; t = s.t, voxel = voxelVolume(s.phi, config.h),
           smeared = smearedVolume(s.phi, config.h), cut = g.volume,
           wall_contacts = g.wall_contacts, dsym_smeared = sym.smeared,
           dsym_voxel = sym.voxel, jaccard = sym.jaccard, overlap = sym.overlap,
           iface, points = g.points,
           phi_min = minimum(s.phi), phi_max = maximum(s.phi),
           weno5_min = minimum(s.c_weno5), weno5_max = maximum(s.c_weno5),
           cip_min = minimum(s.c_cip), cip_max = maximum(s.c_cip))
    end
    errors = (; phi = relativeErrors(snaps[end].phi, config.phi0),
                c_weno5 = relativeErrors(snaps[end].c_weno5, config.c0),
                c_cip = relativeErrors(snaps[end].c_cip, config.c0))
    finite = all(all(isfinite, a) for s in snaps for a in (s.phi, s.c_weno5, s.c_cip))
    return (; nodes = config.nodes, dt = config.dt, steps = config.steps, cfl = config.cfl,
              h = config.h, vmax = config.vmax, amplitude = config.amplitude,
              t_end = config.t_end, seconds = run.seconds, bytes = run.bytes,
              gctime = run.gctime, exact_volume = EXACT_VOLUME, finite, geom, errors)
end

# ------------------------------------------------------------------- allocation study
"""
    measure(op, reps)

Steady-state cost of `op`: one warm-up call, then the minimum `@allocated` and the median
wall time over `reps` further calls. BenchmarkTools is not in the pinned validation
environment, so repeated `@allocated`/`@timed` is used; for these kernels the allocation
counter is deterministic and identical across repetitions.
"""
function measure(op::F, reps::Int = 5) where {F}
    op()
    bytes = minimum(Base.@allocated(op()) for _ in 1:reps)
    times = sort!([Base.@timed(op()).time for _ in 1:reps])
    return (; bytes, seconds = times[cld(reps, 2)])
end

"""
    allocationStudy(config)

Isolated steady-state cost of every operation of the loop, measured after compilation and
cache initialisation, on dedicated work fields so the run state is untouched. `cip_scalar`
is measured with the source buffer refreshed exactly as in the loop, because CIP rebuilds its
interpolation profile whenever the source changed since the previous call — measuring the
bare call twice would silently skip `_reconstructProfile!`.
"""
function allocationStudy(config)
    n = length(config.points)
    dt, weno = config.dt, WENO5Cache()
    geom = Field(config.discretization, Dict(:init => copy(config.phi0)))
    src = Field(config.discretization, Dict(:init => copy(config.c0)))
    ref = Field(config.discretization, Dict(:init => copy(config.c0)))
    wscal = Field(config.discretization, Dict(:init => copy(config.c0)))
    dst = Field(config.discretization, Dict(:init => copy(config.c0)))
    v = dynamics(config, 0.5)
    cip = CIPCache(vec(copy(config.c0)), config.grid)
    # compile the CIP path and bring its cache into the loop's steady state
    lagrangianTransport!(dst, src, v, dt; level_id = :init, cache = cip)
    copyto!(src.data[:init], dst.data[:init])

    specs = (
        ("velocity",        () -> (dynamics(config, 0.5); nothing)),
        ("weno5_geometry",  () -> (eulerianTransport!(geom, v, dt; level_id = :init,
                                                     cache = weno); nothing)),
        ("weno5_scalar",    () -> (eulerianTransport!(wscal, v, dt; level_id = :init,
                                                     cache = weno); nothing)),
        ("cip_scalar",      () -> (copyto!(src.data[:init], ref.data[:init]);
                                   lagrangianTransport!(dst, src, v, dt; level_id = :init,
                                                        cache = cip); nothing)),
        # the loop's steady state: the previous call left cache.snapshot == dst, and the
        # caller copied dst back into src, so the interpolation profile is still valid and
        # `_reconstructProfile!` is skipped. This is the honest per-step cost.
        ("cip_scalar_loop", () -> (lagrangianTransport!(dst, src, v, dt; level_id = :init,
                                                        cache = cip);
                                   copyto!(src.data[:init], dst.data[:init]); nothing)),
        ("cip_copyback",    () -> (copyto!(src.data[:init], dst.data[:init]); nothing)),
        ("snapshot_copy",   () -> (copy(geom.data[:init]); nothing)),
        ("cut",             () -> (cut(geom, :init); nothing)),
        # the whole loop body, with the CIP source untouched by the WENO5 scalar exactly as in
        # `advanceScalar!` (they are separate Fields), so this row must equal the sum of parts
        ("full_step",       () -> (v2 = dynamics(config, 0.5);
                                   eulerianTransport!(geom, v2, dt; level_id = :init,
                                                      cache = weno);
                                   eulerianTransport!(wscal, v2, dt; level_id = :init,
                                                      cache = weno);
                                   lagrangianTransport!(dst, src, v2, dt; level_id = :init,
                                                        cache = cip);
                                   copyto!(src.data[:init], dst.data[:init]); nothing)),
    )
    rows = map(specs) do (name, op)
        m = measure(op)
        (; operation = name, bytes = m.bytes, per_node = m.bytes / n,
           seconds = m.seconds, nodes = n)
    end
    @printf("\n%-16s %14s %12s %12s\n", "operation", "bytes", "bytes/node", "seconds")
    for r in rows
        @printf("%-16s %14d %12.1f %12.4g\n", r.operation, r.bytes, r.per_node, r.seconds)
    end
    return rows
end

# ------------------------------------------------------------------- allocation profile
"""
    allocationProfile(config)

Attribute the measured bytes to source locations with `Profile.Allocs` (sample_rate = 1), with
no modification of the package. Each allocation is labelled by the innermost and the outermost
frame that belongs to this checkout, which separates "which function called `similar`" from
"which function is responsible for the whole temporary".
"""
function allocationProfile(config)
    dt = config.dt
    geom = Field(config.discretization, Dict(:init => copy(config.phi0)))
    src = Field(config.discretization, Dict(:init => copy(config.c0)))
    dst = Field(config.discretization, Dict(:init => copy(config.c0)))
    v = dynamics(config, 0.5)
    weno = WENO5Cache()
    cip = CIPCache(vec(copy(config.c0)), config.grid)
    lagrangianTransport!(dst, src, v, dt; level_id = :init, cache = cip)
    copyto!(src.data[:init], dst.data[:init])

    specs = (
        ("weno5_geometry", () -> (eulerianTransport!(geom, v, dt; level_id = :init,
                                                    cache = weno); nothing)),
        ("cip_scalar_loop", () -> (lagrangianTransport!(dst, src, v, dt; level_id = :init,
                                                        cache = cip);
                                   copyto!(src.data[:init], dst.data[:init]); nothing)),
        ("cut",            () -> (cut(geom, :init); nothing)),
    )
    rows = NamedTuple[]
    for (name, op) in specs
        op()
        Profile.Allocs.clear()
        Profile.Allocs.@profile sample_rate = 1.0 op()
        res = Profile.Allocs.fetch()
        tally = Dict{NTuple{2,String},NTuple{2,Int}}()
        for a in res.allocs
            frames = [f for f in a.stacktrace
                      if occursin("EvolvingDomains-Redesign/src", string(f.file))]
            key = isempty(frames) ? ("<outside>", "<outside>") :
                  (string(first(frames).func), string(last(frames).func))
            c, b = get(tally, key, (0, 0))
            tally[key] = (c + 1, b + a.size)
        end
        for (key, (count, bytes)) in sort!(collect(tally); by = kv -> -kv[2][2])
            push!(rows, (; operation = name, inner = key[1], outer = key[2], count, bytes))
        end
        @printf("\nallocation profile, %s: %d allocations, %d bytes\n", name,
                sum(first, values(tally)), sum(last, values(tally)))
        for (key, (count, bytes)) in first(sort!(collect(tally); by = kv -> -kv[2][2]),
                                           min(8, length(tally)))
            @printf("  %8d B  x%-6d %-28s <- %s\n", bytes, count, key[1], key[2])
        end
    end
    return rows
end

# ------------------------------------------------------------------------ persistence
"""
    writeRun(dir, config, m, allocs)

`params.txt` (what a reader needs to reshape everything), `metrics.csv`, `errors.csv`,
`interface.csv`, `allocs.csv`, one `interface_points_t<S>.csv` per snapshot, and `fields.bin`
holding `Float64` records of `nodes³` in column-major order (first axis fastest, i.e.
`reshape(v, nodes, nodes, nodes)` is indexed (x, y, z)):

    1  phi0             2  c0
    3  phi       @ t0   4  c_weno5   @ t0   5  c_cip  @ t0
    ... repeated at t = 0.25, 0.5, 0.75, 1 (records 6..17)
"""
function writeRun(dir, config, m, run, allocs, profile)
    mkpath(dir)
    open(joinpath(dir, "params.txt"), "w") do io
        for (k, val) in (("nodes", config.nodes), ("steps", config.steps),
                         ("dt", config.dt), ("cfl", config.cfl), ("h", config.h),
                         ("vmax", config.vmax), ("amplitude", config.amplitude),
                         ("t_end", config.t_end), ("records", 2 + 3length(SNAP_TIMES)),
                         ("times", join(SNAP_TIMES, " ")))
            println(io, k, " ", val)
        end
    end
    open(joinpath(dir, "fields.bin"), "w") do io
        write(io, config.phi0)
        write(io, config.c0)
        for s in run.snapshots
            write(io, s.phi)
            write(io, s.c_weno5)
            write(io, s.c_cip)
        end
    end
    open(joinpath(dir, "metrics.csv"), "w") do io
        println(io, "time,voxel,smeared,cutcell,wall_contacts,phi_min,phi_max,",
                    "weno5_min,weno5_max,cip_min,cip_max,dsym_smeared,dsym_voxel,jaccard,overlap")
        for r in m.geom
            println(io, join((r.t, r.voxel, r.smeared,
                              r.cut === missing ? "" : r.cut,
                              r.wall_contacts === missing ? "" : r.wall_contacts,
                              r.phi_min, r.phi_max, r.weno5_min, r.weno5_max,
                              r.cip_min, r.cip_max, r.dsym_smeared, r.dsym_voxel,
                              r.jaccard, r.overlap), ','))
        end
    end
    open(joinpath(dir, "errors.csv"), "w") do io
        println(io, "quantity,l1,l2")
        for (name, e) in (("phi", m.errors.phi), ("c_weno5", m.errors.c_weno5),
                          ("c_cip", m.errors.c_cip))
            println(io, join((name, e.l1, e.l2), ','))
        end
    end
    open(joinpath(dir, "interface.csv"), "w") do io
        println(io, "time,mean_fwd,max_fwd,mean_bwd,max_bwd,n_fwd,n_bwd")
        for r in m.geom
            i = r.iface
            println(io, join((r.t, i.mean_fwd, i.max_fwd, i.mean_bwd, i.max_bwd,
                              i.n_fwd, i.n_bwd), ','))
        end
    end
    open(joinpath(dir, "allocs.csv"), "w") do io
        println(io, "operation,bytes,bytes_per_node,seconds,nodes")
        for r in allocs
            println(io, join((r.operation, r.bytes, r.per_node, r.seconds, r.nodes), ','))
        end
    end
    open(joinpath(dir, "allocs_profile.csv"), "w") do io
        println(io, "operation,outer_frame,inner_frame,count,bytes")
        for r in profile
            println(io, join((r.operation, r.outer, r.inner, r.count, r.bytes), ','))
        end
    end
    for (g, t) in zip(m.geom, SNAP_TIMES)
        suffix = replace(string(t), r"\.0$" => "")
        open(joinpath(dir, "interface_points_t$(suffix).csv"), "w") do io
            println(io, "x,y,z")
            for p in g.points
                println(io, join((p[1], p[2], p[3]), ','))
            end
        end
    end
    return dir
end

# ------------------------------------------------------------------------------ report
function report(io::IO, m, allocs, profile)
    @printf(io, "\n=== StrongDeformation3D — %d³ nodes, %d steps, dt = %.6g (CFL %.2f), A = %.2f ===\n",
            m.nodes, m.steps, m.dt, m.cfl, m.amplitude)
    @printf(io, "wall clock %.2f s (GC %.2f s), %.3g GB allocated, all values finite: %s\n",
            m.seconds, m.gctime, m.bytes / 1e9, m.finite)
    @printf(io, "peak RSS after the time loop: %.3f GB (Sys.maxrss high-water mark, taken before the allocation study)\n",
            m.maxrss / 2^30)
    @printf(io, "exact initial volume 4πabc/3 = %.9g; wall contacts: %s\n",
            m.exact_volume,
            any(r -> r.wall_contacts === missing ? false : r.wall_contacts != 0, m.geom) ?
                "YES (cut volume invalid)" : "none")
    @printf(io, "\n%-6s %11s %11s %11s %8s %8s %8s %8s\n", "t", "V_voxel", "V_smeared",
            "V_cutcell", "D_sym_s", "D_sym_v", "Jaccard", "overlap")
    for r in m.geom
        @printf(io, "%-6.2f %11.6g %11.6g %11s %8.3f %8.3f %8.3f %8.3f\n", r.t, r.voxel,
                r.smeared, r.cut === missing ? "n/a" : @sprintf("%.6g", r.cut),
                r.dsym_smeared, r.dsym_voxel, r.jaccard, r.overlap)
    end
    @printf(io, "phi range: %s\n",
            join((@sprintf("t=%.2f [%.4f, %.4f]", r.t, r.phi_min, r.phi_max)
                  for r in m.geom), "  "))
    for r in m.geom
        @printf(io, "t=%.2f  volume drift %+.3e (voxel) / %+.3e (smeared) / %s (cut)\n",
                r.t, (r.voxel - m.geom[1].voxel) / m.geom[1].voxel,
                (r.smeared - m.geom[1].smeared) / m.geom[1].smeared,
                r.cut === missing ? "n/a" :
                    @sprintf("%+.3e", (r.cut - m.geom[1].cut) / m.geom[1].cut))
    end
    @printf(io, "\ninterface distance (normalised by a_min = %.3g):\n", AMIN)
    @printf(io, "%-6s %10s %10s %10s %10s %8s\n", "t", "mean_fwd", "max_fwd",
            "mean_bwd", "max_bwd", "n_pts")
    for r in m.geom
        i = r.iface
        @printf(io, "%-6.2f %10.3e %10.3e %10.3e %10.3e %8d\n", r.t, i.mean_fwd, i.max_fwd,
                i.mean_bwd, i.max_bwd, i.n_fwd)
    end
    for (name, e) in (("phi    ", m.errors.phi), ("c/weno5", m.errors.c_weno5),
                      ("c/cip  ", m.errors.c_cip))
        @printf(io, "return error %s  L1 %.4e   L2 %.4e\n", name, e.l1, e.l2)
    end
    @printf(io, "\nallocation, steady state (nodes = %d, bytes/node/step):\n", m.nodes)
    for r in allocs
        @printf(io, "  %-16s %12d B  %9.1f B/node  %10.4g s\n", r.operation, r.bytes,
                r.per_node, r.seconds)
    end
    total = sum(r.bytes for r in allocs if r.operation ∉ ("full_step", "cut"))
    @printf(io, "  %-16s %12d B  (sum of every per-step operation)\n", "sum_parts", total)
    @printf(io, "  loop total       %12d B over %d steps\n", m.bytes, m.steps)
    @printf(io, "  loop per step    %.4g B/step = %.1f B/node/step = %.4g s\n",
            m.bytes / m.steps, m.bytes / m.nodes^3 / m.steps, m.seconds / m.steps)
    @printf(io, "  cut per call     %.4g B = %.1f B/node = %.4g s (diagnostic, not in the loop)\n",
            allocs[findfirst(r -> r.operation == "cut", allocs)].bytes,
            allocs[findfirst(r -> r.operation == "cut", allocs)].per_node,
            allocs[findfirst(r -> r.operation == "cut", allocs)].seconds)
    @printf(io, "\nallocation profile (Profile.Allocs, one call each):\n")
    for r in first(profile, min(12, length(profile)))
        @printf(io, "  %-14s %10d B x%-5d %-26s <- %s\n", r.operation, r.bytes, r.count,
                r.inner, r.outer)
    end
    return io
end

# ------------------------------------------------------------------------------- entry
function main()
    # finest first: Sys.maxrss() is a process high-water mark, so the fine case's reading is
    # only free of the profiler's own footprint if nothing bigger ran before it
    resolutions = sort(isempty(ARGS) ? [48, 64] : [parse(Int, a) for a in ARGS]; rev = true)
    # SD3D_OUT lets a refinement run (e.g. a CFL-halving study) write next to the main
    # results without touching them; the visualization module only scans "output".
    outdir = joinpath(@__DIR__, "output" * get(ENV, "SD3D_OUT", ""))
    mkpath(outdir)
    dirs = String[]
    conv = Tuple{Int,String,Float64,Float64}[]
    checks = String[]

    let buf = IOBuffer()                       # Stage A evidence, written to output/
        verifySource(buf)
        text = take!(buf) |> String
        print(stdout, text)
        push!(checks, text)
    end

    walltime = @elapsed for nodes in resolutions
        config = configuration(nodes)
        let buf = IOBuffer()                   # Stage B evidence
            verifyFlow(config, buf)
            verifyVelocityOrdering(config, buf)
            text = take!(buf) |> String
            print(stdout, text)
            push!(checks, text)
        end
        run = simulate(config)
        rss = Sys.maxrss()                      # before diagnostics/allocation study allocate
        m = merge(diagnostics(config, run), (; maxrss = rss))
        allocs = allocationStudy(config)
        profile = allocationProfile(config)
        report(stdout, m, allocs, profile)
        dir = joinpath(outdir, "n$nodes")
        writeRun(dir, config, m, run, allocs, profile)
        open(joinpath(dir, "metrics.txt"), "w") do io
            report(io, m, allocs, profile)
        end
        push!(dirs, dir)
        append!(conv, ((nodes, "phi", m.errors.phi.l1, m.errors.phi.l2),
                       (nodes, "c_weno5", m.errors.c_weno5.l1, m.errors.c_weno5.l2),
                       (nodes, "c_cip", m.errors.c_cip.l1, m.errors.c_cip.l2)))
    end
    @printf("\ntotal simulation wall clock: %.2f s\n", walltime)

    write(joinpath(outdir, "verification.txt"), join(checks, "\n"))

    # return-error-vs-resolution table for the visualization module
    open(joinpath(outdir, "convergence.csv"), "w") do io
        println(io, "nodes,quantity,l1,l2")
        for row in conv
            println(io, join(row, ','))
        end
    end
    return dirs
end

main()
