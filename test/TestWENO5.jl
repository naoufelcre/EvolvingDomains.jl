# Focused ED1 regressions for WENO5 sign selection.
#
# ED1 changes `weno5_rhs!` so that it reads the frozen velocity first and
# evaluates only the upwind reconstruction that the sign selects, instead of
# eagerly evaluating all four reconstructions and discarding two.  The operator,
# the `v > 0` branch, the final arithmetic, the RK3 coefficients, the full grid
# (no narrow band), the cache and the Float64 precision must all be unchanged.
#
# The reference below is a *frozen eager four-reconstruction* implementation,
# written out independently of `src/Geometric/Stencils.jl` and of the optimized
# production kernel.  It reproduces the pre-ED1 behaviour exactly (evaluate
# `weno5⁻` and `weno5⁺` in both directions, then select) so that a change to the
# production arithmetic or to the sign rule is caught by a real float-bit
# comparison, rather than by a comparison of the optimized code with itself.
#
# Runnable directly:
#     julia --project=. --check-bounds=yes -t 1 test/TestWENO5.jl
#     julia --project=. --check-bounds=yes -t 4 test/TestWENO5.jl
# It is also wired into test/runtests.jl.

using Test
using EvolvingDomains
using EvolvingDomains.Geometric
using EvolvingDomains.Kinematic
using Gridap
using Gridap.TensorValues

const W5 = EvolvingDomains.Kinematic.WENO5

# ---------------------------------------------------------------------------
# Frozen independent eager four-reconstruction reference (pre-ED1 behaviour).
# The formulas are copied from the mathematics of the operator, not called from
# the package, so a change to Stencils.jl cannot silently redefine the reference.
# ---------------------------------------------------------------------------

@inline _ref_offset(dim::Int) = CartesianIndex(ntuple(k -> k == dim ? 1 : 0, 2))

@inline function _ref_get(data::Vector{Float64}, nx::Int, ny::Int, I::CartesianIndex{2})
    i = clamp(I[1], 1, nx)   # same clamped boundary handling as CartesianMeshField
    j = clamp(I[2], 1, ny)
    @inbounds data[i + (j - 1) * nx]
end

@inline _ref_Dminus(data, nx, ny, I, dim, h) =
    (_ref_get(data, nx, ny, I) - _ref_get(data, nx, ny, I - _ref_offset(dim))) / h

@inline _ref_Dplus(data, nx, ny, I, dim, h) =
    (_ref_get(data, nx, ny, I + _ref_offset(dim)) - _ref_get(data, nx, ny, I)) / h

@inline function _ref_weno5_weights(v1, v2, v3, v4, v5)
    ϵ = 1e-6
    S1 = 13/12 * (v1 - 2*v2 + v3)^2 + 1/4 * (v1 - 4*v2 + 3*v3)^2
    S2 = 13/12 * (v2 - 2*v3 + v4)^2 + 1/4 * (v2 - v4)^2
    S3 = 13/12 * (v3 - 2*v4 + v5)^2 + 1/4 * (3*v3 - 4*v4 + v5)^2

    α1 = 0.1 / (S1 + ϵ)^2
    α2 = 0.6 / (S2 + ϵ)^2
    α3 = 0.3 / (S3 + ϵ)^2
    sum_α = α1 + α2 + α3

    return α1/sum_α, α2/sum_α, α3/sum_α
end

function _ref_weno5_minus(data, nx, ny, I, dim, h)
    offset = _ref_offset(dim)
    v1 = _ref_Dminus(data, nx, ny, I - 2*offset, dim, h)
    v2 = _ref_Dminus(data, nx, ny, I - offset,   dim, h)
    v3 = _ref_Dminus(data, nx, ny, I,            dim, h)
    v4 = _ref_Dminus(data, nx, ny, I + offset,   dim, h)
    v5 = _ref_Dminus(data, nx, ny, I + 2*offset, dim, h)

    w1, w2, w3 = _ref_weno5_weights(v1, v2, v3, v4, v5)

    q1 = v1/3 - 7*v2/6 + 11*v3/6
    q2 = -v2/6 + 5*v3/6 + v4/3
    q3 = v3/3 + 5*v4/6 - v5/6

    return w1*q1 + w2*q2 + w3*q3
end

function _ref_weno5_plus(data, nx, ny, I, dim, h)
    offset = _ref_offset(dim)
    v1 = _ref_Dplus(data, nx, ny, I + 2*offset, dim, h)
    v2 = _ref_Dplus(data, nx, ny, I + offset,   dim, h)
    v3 = _ref_Dplus(data, nx, ny, I,            dim, h)
    v4 = _ref_Dplus(data, nx, ny, I - offset,   dim, h)
    v5 = _ref_Dplus(data, nx, ny, I - 2*offset, dim, h)

    w1, w2, w3 = _ref_weno5_weights(v1, v2, v3, v4, v5)

    q1 = v1/3 - 7*v2/6 + 11*v3/6
    q2 = -v2/6 + 5*v3/6 + v4/3
    q3 = v3/3 + 5*v4/6 - v5/6

    return w1*q1 + w2*q2 + w3*q3
end

# Eager four-reconstruction RHS: exactly the pre-ED1 production kernel.
function ref_rhs_eager!(rhs::Vector{Float64}, data::Vector{Float64},
                        info::CartesianGridInfo, velocity::AbstractVector)
    nx, ny = info.dims
    dx, dy = info.spacing
    for j in 1:ny
        @inbounds for i in 1:nx
            I = CartesianIndex(i, j)
            idx = i + (j - 1) * nx

            dx_L = _ref_weno5_minus(data, nx, ny, I, 1, dx)
            dx_R = _ref_weno5_plus(data, nx, ny, I, 1, dx)
            dy_L = _ref_weno5_minus(data, nx, ny, I, 2, dy)
            dy_R = _ref_weno5_plus(data, nx, ny, I, 2, dy)

            v = velocity[idx]
            grad_x = (v[1] > 0) ? dx_L : dx_R
            grad_y = (v[2] > 0) ? dy_L : dy_R

            rhs[idx] = -(v[1] * grad_x + v[2] * grad_y)
        end
    end
    return rhs
end

# Reference SSP-RK3, arithmetic copied verbatim from weno5_step!.
mutable struct RefCache
    rhs::Vector{Float64}
    stage::Vector{Float64}
    phi0::Vector{Float64}
end
RefCache() = RefCache(Float64[], Float64[], Float64[])

function ref_step_eager!(phi::Vector{Float64}, info::CartesianGridInfo,
                         velocity::AbstractVector, dt::Float64, c::RefCache)
    n = length(phi)
    if length(c.rhs) != n
        c.rhs = zeros(Float64, n)
        c.stage = zeros(Float64, n)
        c.phi0 = zeros(Float64, n)
    end
    rhs, stage, phi0 = c.rhs, c.stage, c.phi0

    @. phi0 = phi
    ref_rhs_eager!(rhs, phi, info, velocity)
    @. stage = phi + dt * rhs
    @. phi = stage
    ref_rhs_eager!(rhs, phi, info, velocity)
    @. stage = 0.75 * phi0 + 0.25 * stage + 0.25 * dt * rhs
    @. phi = stage
    ref_rhs_eager!(rhs, phi, info, velocity)
    @. phi = (1.0/3.0) * phi0 + (2.0/3.0) * stage + (2.0/3.0) * dt * rhs
    return phi
end

# ---------------------------------------------------------------------------
# Fixtures and helpers
# ---------------------------------------------------------------------------

bits(x::Vector{Float64}) = reinterpret(UInt64, x)

weno_grid_meta(origin, spacing, dims) =
    CartesianGridInfo(origin, spacing, dims, (dims[1] - 1, dims[2] - 1))

# Tiny, shifted, anisotropic and non-square grids.  Every fixture includes the
# domain boundary nodes, so the clamped-boundary path is always exercised.
const FIXTURES = [
    ("tiny_2x2",       (0.0, 0.0),    (1.0, 1.0),   (2, 2)),
    ("tiny_3x3",       (0.0, 0.0),    (1.0, 1.0),   (3, 3)),
    ("small_2x4",      (0.0, 0.0),    (1.0, 0.5),   (2, 4)),
    ("shifted_aniso",  (1.5, -3.0),   (0.7, 0.3),   (8, 4)),
    ("tall_aniso",     (-2.0, 0.5),   (2.0, 0.1),   (5, 9)),
    ("wide_aniso",     (0.0, 0.0),    (0.25, 0.8),  (13, 3)),
]

fixture_points(info) = [(info.origin[1] + (i - 1) * info.spacing[1],
                         info.origin[2] + (j - 1) * info.spacing[2])
                        for j in 1:info.dims[2] for i in 1:info.dims[1]]

const FIELDS = [
    ("smooth", pts -> [sin(2π * p[1]) * cos(2π * p[2]) + 0.3 * p[1]^2 - 0.2 * p[2]
                       for p in pts]),
    ("kink",   pts -> [abs(p[1] - 0.4) - 0.7 * abs(p[2] - 0.6) + 0.05 * p[1] * p[2]
                       for p in pts]),
    ("step",   pts -> [(p[1] < 0.5 ? -1.0 : 1.0) + 0.1 * p[2] for p in pts]),
]

const VEL_MODES = [
    ("mixed",        pts -> [VectorValue(sin(4π * p[1]), cos(3π * p[2])) for p in pts]),
    ("zeros",        pts -> [VectorValue(0.0, 0.0) for p in pts]),
    ("neg_zeros",    pts -> [VectorValue(-0.0, -0.0) for p in pts]),
    ("onehot_x+",    pts -> [VectorValue(1.0, 0.0) for p in pts]),
    ("onehot_x-",    pts -> [VectorValue(-1.0, 0.0) for p in pts]),
    ("onehot_y+",    pts -> [VectorValue(0.0, 1.0) for p in pts]),
    ("onehot_y-",    pts -> [VectorValue(0.0, -1.0) for p in pts]),
    ("sign_pattern", pts -> [VectorValue(sign(sin(4π * p[1])) * 0.3,
                                         sign(cos(3π * p[2])) * 0.2) for p in pts]),
    ("tiny_vel",     pts -> [VectorValue(nextfloat(0.0), -nextfloat(0.0)) for p in pts]),
]

# ---------------------------------------------------------------------------
# 1. RHS bit-identity: optimized production kernel vs frozen eager reference
# ---------------------------------------------------------------------------
@testset "ED1 RHS bit-identity (package vs frozen eager reference)" begin
    n_checks = 0
    for (fxname, origin, spacing, dims) in FIXTURES
        info = weno_grid_meta(origin, spacing, dims)
        pts = fixture_points(info)
        nx, ny = info.dims
        @assert length(pts) == nx * ny
        for (fname, fbuild) in FIELDS
            data = fbuild(pts)
            for (vname, vbuild) in VEL_MODES
                vel = vbuild(pts)
                r_pkg = fill(NaN, nx * ny)
                r_ref = fill(NaN, nx * ny)
                W5.weno5_rhs!(r_pkg, CartesianMeshField(copy(data), info), vel)
                ref_rhs_eager!(r_ref, copy(data), info, vel)
                @test bits(r_pkg) == bits(r_ref)
                n_checks += 1
            end
        end
    end
    @test n_checks == length(FIXTURES) * length(FIELDS) * length(VEL_MODES)
end

# ---------------------------------------------------------------------------
# 2. Boundary nodes participate explicitly (not just the interior)
# ---------------------------------------------------------------------------
@testset "ED1 boundary-node RHS bit-identity" begin
    for (fxname, origin, spacing, dims) in FIXTURES
        info = weno_grid_meta(origin, spacing, dims)
        pts = fixture_points(info)
        nx, ny = info.dims
        data = FIELDS[2][2](pts)                       # kinked field
        vel = VEL_MODES[8][2](pts)                     # sign pattern
        r_pkg = fill(NaN, nx * ny)
        r_ref = fill(NaN, nx * ny)
        W5.weno5_rhs!(r_pkg, CartesianMeshField(copy(data), info), vel)
        ref_rhs_eager!(r_ref, copy(data), info, vel)

        bnd = Int[]
        for j in 1:ny, i in 1:nx
            if i == 1 || i == nx || j == 1 || j == ny
                push!(bnd, i + (j - 1) * nx)
            end
        end
        @test !isempty(bnd)
        @test bits(r_pkg[bnd]) == bits(r_ref[bnd])     # real bit compare on boundary
    end
end

# ---------------------------------------------------------------------------
# 3. Zero and signed-zero velocity: exact output bits, including sign of zero
# ---------------------------------------------------------------------------
@testset "ED1 zero and signed-zero velocity bits" begin
    for (fxname, origin, spacing, dims) in FIXTURES
        info = weno_grid_meta(origin, spacing, dims)
        pts = fixture_points(info)
        nx, ny = info.dims
        data = FIELDS[2][2](pts)                       # kinked, opposite one-sided slopes
        for vpair in ((0.0, 0.0), (-0.0, -0.0), (0.0, -0.0), (-0.0, 0.0))
            vel = [VectorValue(vpair[1], vpair[2]) for _ in 1:(nx * ny)]
            r_pkg = fill(NaN, nx * ny)
            r_ref = fill(NaN, nx * ny)
            W5.weno5_rhs!(r_pkg, CartesianMeshField(copy(data), info), vel)
            ref_rhs_eager!(r_ref, copy(data), info, vel)
            @test bits(r_pkg) == bits(r_ref)
            @test all(x -> x == 0.0, r_pkg)           # zero times finite gradient
        end
    end
end

# ---------------------------------------------------------------------------
# 4. Repeated full SSP-RK3 steps, separate caches
# ---------------------------------------------------------------------------
@testset "ED1 repeated SSP-RK3 bit-identity" begin
    step_counts = 0
    for (fxname, origin, spacing, dims) in FIXTURES
        info = weno_grid_meta(origin, spacing, dims)
        pts = fixture_points(info)
        for (fname, fbuild) in FIELDS
            data = fbuild(pts)
            for (vname, vbuild) in VEL_MODES
                vel = vbuild(pts)
                dt = 0.02
                steps = 5
                a = copy(data)
                b = copy(data)
                ca = WENO5Cache()
                cb = RefCache()
                for _ in 1:steps
                    weno5_step!(a, info, vel, dt, ca)
                    ref_step_eager!(b, info, vel, dt, cb)
                end
                @test bits(a) == bits(b)
                step_counts += 1
            end
        end
    end
    @test step_counts == length(FIXTURES) * length(FIELDS) * length(VEL_MODES)
end

# dt sweeps including dt = 0 and a negative dt; coefficient order must survive.
@testset "ED1 SSP-RK3 dt edge cases" begin
    info = weno_grid_meta((0.3, -1.1), (0.4, 0.9), (6, 5))
    pts = fixture_points(info)
    vel = VEL_MODES[1][2](pts)
    data = FIELDS[1][2](pts)
    for dt in (0.0, -0.02, 1.0e-12, 0.5)
        a = copy(data)
        b = copy(data)
        weno5_step!(a, info, vel, dt, WENO5Cache())
        ref_step_eager!(b, info, vel, dt, RefCache())
        @test bits(a) == bits(b)
    end
end

# ---------------------------------------------------------------------------
# 5. Cache reuse across dimension changes (resizing) and cache-less path
# ---------------------------------------------------------------------------
@testset "ED1 cache resize across dimension changes" begin
    sequence = [weno_grid_meta((0.0, 0.0), (1.0, 1.0), (2, 2)),
                weno_grid_meta((0.0, 0.0), (0.5, 0.25), (8, 4)),
                weno_grid_meta((-1.0, 2.0), (0.3, 0.7), (3, 9)),
                weno_grid_meta((0.0, 0.0), (0.5, 0.25), (8, 4))]  # back to a previous size
    shared = WENO5Cache()
    for info in sequence
        n = prod(info.dims)
        pts = fixture_points(info)
        data = FIELDS[1][2](pts)
        vel = VEL_MODES[1][2](pts)
        a = copy(data)
        b = copy(data)
        weno5_step!(a, info, vel, 0.01, shared)
        ref_step_eager!(b, info, vel, 0.01, RefCache())
        @test length(shared.rhs) == n
        @test length(shared.stage) == n
        @test length(shared.phi0) == n
        @test bits(a) == bits(b)
    end
end

@testset "ED1 cache-less singleton path" begin
    # Uses WENO_CACHE.  Run serially (the module singleton is pre-existing and
    # shared); only one weno5_step! call is in flight at a time.
    for info in (weno_grid_meta((0.0, 0.0), (1.0, 1.0), (4, 4)),
                 weno_grid_meta((0.2, -0.4), (0.6, 0.3), (7, 5)))
        pts = fixture_points(info)
        data = FIELDS[3][2](pts)
        vel = VEL_MODES[8][2](pts)
        a = copy(data)
        b = copy(data)
        weno5_step!(a, info, vel, 0.015)                 # no cache argument
        ref_step_eager!(b, info, vel, 0.015, RefCache())
        @test bits(a) == bits(b)
    end
end

# ---------------------------------------------------------------------------
# 6. Geometry caller: advance! (uses the geometry-owned cache and scratch)
# ---------------------------------------------------------------------------
@testset "ED1 geometry advance! caller bit-identity" begin
    model = CartesianDiscreteModel((0.0, 1.0, 0.0, 1.0), (6, 5))
    info = grid_info(model)
    pts = vec(collect(Gridap.Geometry.get_node_coordinates(model)))
    phi0 = [sin(2π * p[1]) * cos(2π * p[2]) + 0.3 * p[1]^2 - 0.2 * p[2] for p in pts]

    vel_src = StaticFunctionVelocity(x -> VectorValue(sin(4π * x[1]), -cos(3π * x[2])))
    vel = sample_velocity(vel_src, info, 0.0)

    geom = EvolvingDiscreteGeometry(vec(copy(phi0)), model)
    dt = 0.02
    nsteps = 3
    for _ in 1:nsteps
        advance!(geom, vel, dt)
    end

    ref = copy(phi0)
    rc = RefCache()
    for _ in 1:nsteps
        ref_step_eager!(ref, info, vel, dt, rc)
    end

    @test geom.cache.weno_cache.rhs !== nothing
    @test bits(current_levelset(geom)) == bits(ref)
    @test bits(geom.levelset) == bits(ref)
end

# ---------------------------------------------------------------------------
# 7. Bounded physics checks on the PRODUCTION kernels (not the reference).
#    Boundary equivalence is checked above; these are smooth-interior checks and
#    are deliberately kept separate.  The interior margin is 9 nodes, the exact
#    SSP-RK3 dependency reach (3 RHS evaluations x 3 stencil nodes), so clamping
#    at the domain wall cannot leak into the measured region.
# ---------------------------------------------------------------------------
@testset "ED1 production smooth-field spatial RHS order on interior" begin
    function rhs_err(n)
        info = weno_grid_meta((0.0, 0.0), (1.0 / (n - 1), 1.0 / (n - 1)), (n, n))
        pts = fixture_points(info)
        f(x, y) = sin(2π * x) * cos(2π * y)
        phi = CartesianMeshField([f(p...) for p in pts], info)
        vel = fill(VectorValue(1.0, 1.0), length(pts))
        rhs = zeros(length(pts))
        W5.weno5_rhs!(rhs, phi, vel)
        nx, ny = info.dims
        e = 0.0
        for j in 10:ny-9, i in 10:nx-9            # +-9 margin
            k = i + (j - 1) * nx
            x, y = pts[k]
            exact = -(2π * cos(2π * x) * cos(2π * y) - 2π * sin(2π * x) * sin(2π * y))
            e = max(e, abs(rhs[k] - exact))
        end
        return e
    end
    errs = [rhs_err(n) for n in (32, 64, 128)]
    rates = [log2(errs[i] / errs[i + 1]) for i in 1:2]
    @test errs[end] < errs[1] / 100
    @test rates[end] > 4.5
end

@testset "ED1 production one-step RK3 advection dt sweep" begin
    # One full SSP-RK3 step of phi_t + 0.3 phi_x = 0 on a smooth field.  As
    # dt -> 0 the one-step error is O(dt^4) + O(dt*h^5): the early sweep is
    # time-dominated (rate ~4), the late sweep approaches the spatial floor and
    # the rate drops.  This is a local one-step trend; it is NOT a global-order
    # claim (global fixed-T error would be O(dt^3) + O(h^5)).
    n = 128
    info = weno_grid_meta((0.0, 0.0), (1.0 / (n - 1), 1.0 / (n - 1)), (n, n))
    pts = fixture_points(info)
    init(x, y) = sin(2π * x) * cos(2π * y)
    phi0 = [init(p...) for p in pts]
    vel = fill(VectorValue(0.3, 0.0), length(pts))
    nx, ny = info.dims
    function step_err(dt)
        phi = copy(phi0)
        weno5_step!(phi, info, vel, dt, WENO5Cache())
        e = 0.0
        for j in 10:ny-9, i in 10:nx-9
            k = i + (j - 1) * nx
            x, y = pts[k]
            e = max(e, abs(phi[k] - init(x - dt * 0.3, y)))
        end
        return e
    end
    dts = [0.08, 0.04, 0.02, 0.01, 0.005]
    errs = [step_err(dt) for dt in dts]
    rates = [log2(errs[i] / errs[i + 1]) for i in 1:4]
    @test rates[1] > 3.5 && rates[2] > 3.5   # time-dominated regime, close to local order 4
    @test errs[end] < errs[1] / 1.0e3
end

@testset "ED1 production linear-profile exactness" begin
    # WENO5 reproduces a linear field exactly: phi = x + y, v = (0.3, -0.2) gives
    # the nonzero constant RHS -(0.3 - 0.2), so a one-step update must equal
    # data - dt*(0.3 - 0.2).  A cancelling velocity would hide RK stage errors.
    n = 40
    info = weno_grid_meta((0.0, 0.0), (1.0 / (n - 1), 1.0 / (n - 1)), (n, n))
    pts = fixture_points(info)
    data = [p[1] + p[2] for p in pts]
    c = 0.3 - 0.2                    # nonzero net advection coefficient
    vel = fill(VectorValue(0.3, -0.2), length(pts))
    nx, ny = info.dims
    rhs = zeros(length(pts))
    W5.weno5_rhs!(rhs, CartesianMeshField(copy(data), info), vel)
    e_rhs = 0.0
    for j in 10:ny-9, i in 10:nx-9
        e_rhs = max(e_rhs, abs(rhs[i + (j - 1) * nx] - (-c)))
    end
    dt = 0.05
    phi = copy(data)
    weno5_step!(phi, info, vel, dt, WENO5Cache())
    e_step = 0.0
    for j in 10:ny-9, i in 10:nx-9
        k = i + (j - 1) * nx
        e_step = max(e_step, abs(phi[k] - (data[k] - dt * c)))
    end
    @test e_rhs < 1.0e-12
    @test e_step < 1.0e-12
end
