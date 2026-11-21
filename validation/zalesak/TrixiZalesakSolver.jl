# Optional Trixi DGSEM solver and minimal scalar Field exchange for Zalesak.
# Used by runZalesak2DTrixi.jl; run that script rather than this file.
module TrixiZalesakSolver

using Trixi
using StaticArrays: SVector

const OMEGA = 20π

"""Three-component state: transported scalar and frozen nodal velocity components."""
struct FrozenVelocityAdvection2D <: Trixi.AbstractEquations{2,3} end

@inline Trixi.varnames(::typeof(Trixi.cons2cons), ::FrozenVelocityAdvection2D) = ("phi", "vx", "vy")
@inline Trixi.varnames(::typeof(Trixi.cons2prim), ::FrozenVelocityAdvection2D) = ("phi", "vx", "vy")
@inline Trixi.have_constant_speed(::FrozenVelocityAdvection2D) = Trixi.False()
@inline Trixi.max_abs_speeds(u, ::FrozenVelocityAdvection2D) = (abs(u[2]), abs(u[3]))
@inline Trixi.max_abs_speed_naive(ul, ur, d::Integer, ::FrozenVelocityAdvection2D) =
    max(abs(ul[d+1]), abs(ur[d+1]))
@inline Trixi.max_abs_speed_naive(ul, ur, normal::AbstractVector, eq::FrozenVelocityAdvection2D) =
    max(abs(ul[2]*normal[1] + ul[3]*normal[2]), abs(ur[2]*normal[1] + ur[3]*normal[2]))
@inline Trixi.max_abs_speed(ul, ur, d::Integer, eq::FrozenVelocityAdvection2D) =
    Trixi.max_abs_speed_naive(ul, ur, d, eq)
@inline Trixi.max_abs_speed(ul, ur, normal::AbstractVector, eq::FrozenVelocityAdvection2D) =
    Trixi.max_abs_speed_naive(ul, ur, normal, eq)
@inline Trixi.flux(u, d::Integer, ::FrozenVelocityAdvection2D) =
    d == 1 ? SVector(u[2] * u[1], 0.0, 0.0) : SVector(u[3] * u[1], 0.0, 0.0)
@inline function Trixi.flux(u, normal::AbstractVector, ::FrozenVelocityAdvection2D)
    return SVector((u[2]*normal[1] + u[3]*normal[2]) * u[1], 0.0, 0.0)
end
@inline Trixi.cons2prim(u, ::FrozenVelocityAdvection2D) = u
@inline Trixi.cons2entropy(u, ::FrozenVelocityAdvection2D) = u

# The velocities are prescribed coefficients, so numerical dissipation acts only on phi.
@inline function _scalar_lax_friedrichs(ul, ur, direction, eq::FrozenVelocityAdvection2D)
    fl, fr = Trixi.flux(ul, direction, eq), Trixi.flux(ur, direction, eq)
    λ = Trixi.max_abs_speed(ul, ur, direction, eq)
    return SVector(0.5*(fl[1]+fr[1]) - 0.5*λ*(ur[1]-ul[1]), 0.0, 0.0)
end

_velocity(x, y) = (-OMEGA * (y - 0.5), OMEGA * (x - 0.5))

function _initial_state(x, t, equation::FrozenVelocityAdvection2D, initial_function)
    vx, vy = _velocity(x[1], x[2])
    return SVector(initial_function(x[1], x[2]), vx, vy)
end

"""Make a private scalar bilinear interpolant from a named ED `Field` level."""
function field_to_dg(field, points; level_id::Symbol=:init)
    hasproperty(field, :discretization) && hasproperty(field, :data) ||
        throw(ArgumentError("expected an EvolvingDomains Field"))
    parent = parentmodule(typeof(field))
    g = if isdefined(parent, :grid) && applicable(getfield(parent, :grid), field.discretization, level_id)
        try getfield(parent, :grid)(field.discretization, level_id) catch err
            err isa KeyError || rethrow()
            throw(ArgumentError("unknown Field level $level_id"))
        end
    else
        try field.discretization.levels[level_id] catch err
            err isa KeyError || rethrow()
            throw(ArgumentError("unknown Field level $level_id"))
        end
    end
    length(g.axes) == 2 || throw(ArgumentError("Field-to-DG exchange requires a 2D grid"))
    data = if isdefined(parent, :values) && applicable(getfield(parent, :values), field, level_id)
        try getfield(parent, :values)(field, level_id) catch err
            err isa KeyError || rethrow()
            throw(ArgumentError("unknown Field level $level_id"))
        end
    else
        try field.data[level_id] catch err
            err isa KeyError || rethrow()
            throw(ArgumentError("unknown Field level $level_id"))
        end
    end
    size(data) == map(length, g.axes) || throw(ArgumentError("Field values must have shape (length(x), length(y))"))
    all(v -> v isa Real && isfinite(v), data) || throw(ArgumentError("Field values must be finite real scalars"))
    # Snapshot both axes and data; exchange never borrows mutable Field storage.
    x, y, u = copy(g.axes[1]), copy(g.axes[2]), copy(data)
    sample(a, b) = _bilinear(u, x, y, a, b)
    if points isa Tuple && length(points) == 2 && all(p -> p isa Real, points)
        return sample(points...)
    end
    return map(p -> sample(p[1], p[2]), points)
end

function _bilinear(u, x, y, a, b)
    isfinite(a) && isfinite(b) || throw(ArgumentError("sample coordinates must be finite"))
    first(x) <= a <= last(x) && first(y) <= b <= last(y) ||
        throw(ArgumentError("Field-to-DG sample lies outside the source Field domain"))
    i = clamp(searchsortedlast(x, a), 1, length(x)-1)
    j = clamp(searchsortedlast(y, b), 1, length(y)-1)
    r = clamp((a-x[i])/(x[i+1]-x[i]), 0.0, 1.0)
    s = clamp((b-y[j])/(y[j+1]-y[j]), 0.0, 1.0)
    return (1-r)*(1-s)*u[i,j] + r*(1-s)*u[i+1,j] + (1-r)*s*u[i,j+1] + r*s*u[i+1,j+1]
end

function _lagrange_weights(nodes, ξ)
    n = length(nodes)
    return [prod((ξ - nodes[j]) / (nodes[i] - nodes[j]) for j in 1:n if j != i) for i in 1:n]
end

# At element faces use the element with the lowest Cartesian index (left/bottom trace).
function _dg_value(u, nodes, cells, x, y, a, b)
    (first(x) <= a <= last(x) && first(y) <= b <= last(y)) ||
        throw(ArgumentError("DG-to-grid sample lies outside the Trixi domain"))
    hx, hy = (last(x)-first(x))/cells, (last(y)-first(y))/cells
    ix = clamp(ceil(Int, (a-first(x))/hx), 1, cells)
    iy = clamp(ceil(Int, (b-first(y))/hy), 1, cells)
    e = ix + (iy-1)*cells
    ξ = 2 * (a - (first(x)+(ix-1)*hx))/hx - 1
    η = 2 * (b - (first(y)+(iy-1)*hy))/hy - 1
    wx, wy = _lagrange_weights(nodes, ξ), _lagrange_weights(nodes, η)
    val = zero(eltype(u))
    for j in eachindex(wy), i in eachindex(wx)
        val += wx[i] * wy[j] * u[1, i, j, e]
    end
    return val
end

"""
    run_trixi(grid, initial_function, times; cells=64, polydeg=3, cfl=0.25)

Evolve scalar `phi` under the prescribed divergence-free rotation with Trixi's
structured DGSEM. Velocities are frozen components in the state; only `phi` has
nonzero physical flux. `boundary_condition_do_nothing` supplies the interior trace
at incoming as well as outgoing boundaries (not periodic). Returns x-first nodal
snapshots, solver timing and native DG integral diagnostics.
"""
function run_trixi(output_grid, initial_function, times;
                   cells::Integer=64, polydeg::Integer=3, cfl::Real=0.25)
    length(output_grid.axes) == 2 || throw(ArgumentError("output grid must be 2D"))
    cells > 0 && polydeg >= 1 || throw(ArgumentError("cells must be positive and polydeg ≥ 1"))
    isfinite(cfl) && 0 < cfl <= 0.4 || throw(ArgumentError("cfl must lie in (0, 0.4]"))
    ts = Float64.(collect(times))
    !isempty(ts) && all(isfinite, ts) && issorted(ts) && first(ts) == 0 && last(ts) > 0 ||
        throw(ArgumentError("times must be finite, strictly increasing from 0 to a positive final time"))
    length(unique(ts)) == length(ts) || throw(ArgumentError("times must be strictly increasing"))
    x, y = output_grid.axes
    first(x) >= 0 && last(x) <= 1 && first(y) >= 0 && last(y) <= 1 ||
        throw(ArgumentError("the Trixi rotation domain is [0,1]^2"))

    equations = FrozenVelocityAdvection2D()
    initial_condition(xy, t, eq) = _initial_state(xy, t, eq, initial_function)
    mesh = Trixi.StructuredMesh((cells, cells), (0.0, 0.0), (1.0, 1.0), periodicity=false)
    solver = Trixi.DGSEM(polydeg=polydeg, surface_flux=_scalar_lax_friedrichs)
    semi = Trixi.SemidiscretizationHyperbolic(mesh, equations, initial_condition, solver;
        boundary_conditions=Trixi.boundary_condition_do_nothing)
    ode = Trixi.semidiscretize(semi, (0.0, last(ts)))
    dt = cfl / ((2polydeg + 1) * 20π * cells)
    integrator = Trixi.init(ode, Trixi.SimpleSSPRK33(); dt)
    nodes = solver.basis.nodes
    states = Array{Float64,4}[reshape(copy(integrator.u), 3, polydeg+1, polydeg+1, cells^2)]
    actual_times, steps = [0.0], 0
    elapsed = @elapsed for target in ts[2:end]
        while integrator.t < target - 16eps(target)
            integrator.dtcache = min(dt, target - integrator.t)
            integrator.dt = integrator.dtcache
            Trixi.step!(integrator)
            steps += 1
        end
        isapprox(integrator.t, target; atol=32eps(target), rtol=0) ||
            error("Trixi integration missed requested snapshot time $target")
        push!(actual_times, target)
        push!(states, reshape(copy(integrator.u), 3, polydeg+1, polydeg+1, cells^2))
    end
    snapshots = [ [ _dg_value(u, nodes, cells, (0.0,1.0), (0.0,1.0), a, b)
                    for a in x, b in y ] for u in states ]
    integrals = [Trixi.integrate((v, eq) -> SVector(v[1]), vec(u), semi; normalize=false)[1]
                 for u in states]
    return (; x=collect(x), y=collect(y), times=actual_times, snapshots, elapsed,
            dt, steps, integral=integrals, cells, polydeg, cfl,
            metadata=(equation=:divergence_free_rotation, boundary=:interior_extrapolation,
                      state_variables=(:phi,:vx,:vy), scalar_dofs=(polydeg+1)^2*cells^2,
                      integrator=:Trixi_SimpleSSPRK33, dt_is_maximum=true,
                      timing_scope=:stepping_and_snapshot_copies,
                      compilation_warmup=:caller_dependent,
                      native_integral=:quadrature_integral_of_phi_not_indicator_area))
end

"""
    self_check()

Minimal exchange contract: a constant survives the DG solve, the sampled initial state
reproduces an affine level set, and out-of-domain samples are rejected. This is what makes
a cross-solver difference interpretable; the physics gates live with the case, not here.
"""
function self_check()
    x = collect(range(0.0, 1.0; length=5))
    y = collect(range(0.0, 1.0; length=4))
    out = (axes=(x, y),)
    constant = run_trixi(out, (a, b) -> 2.0, [0.0, 1e-4]; cells=8, polydeg=3, cfl=0.2)
    maximum(abs.(constant.snapshots[end] .- 2.0)) < 1e-10 ||
        error("DG did not preserve a constant field")
    affine = run_trixi(out, (a, b) -> 2a - 3b + 1, [0.0, 0.0 + 1e-5]; cells=8, polydeg=3,
                       cfl=0.2)
    affine.snapshots[1] ≈ [2a - 3b + 1 for a in x, b in y] ||
        error("DG initial sampling is not exact on an affine level set")
    _throws_argument_error(() -> run_trixi((axes=(x, y),), (a, b) -> 1.0, [0.0, 1.0];
                                           cells=0, polydeg=3, cfl=0.2))
    _throws_argument_error(() -> run_trixi(out, (a, b) -> 1.0, [0.0, 1.0]; cells=4,
                                           polydeg=3, cfl=0.5))
    _throws_argument_error(() -> run_trixi((axes=(x .+ 2, y),), (a, b) -> 1.0, [0.0, 1.0];
                                           cells=4, polydeg=3, cfl=0.2))
    true
end

function _throws_argument_error(f)
    try
        f()
    catch err
        err isa ArgumentError && return true
        rethrow()
    end
    error("expected ArgumentError")
end

end # module
