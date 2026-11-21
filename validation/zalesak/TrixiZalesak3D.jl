# Optional independent Trixi 3D DGSEM solver for the canonical Zalesak sphere.
module TrixiZalesak3D

using Trixi
using StaticArrays: SVector

"""Four-component state: transported scalar and frozen nodal velocity components."""
struct FrozenVelocityAdvection3D <: Trixi.AbstractEquations{3,4} end

@inline Trixi.varnames(::typeof(Trixi.cons2cons), ::FrozenVelocityAdvection3D) =
    ("phi", "vx", "vy", "vz")
@inline Trixi.varnames(::typeof(Trixi.cons2prim), ::FrozenVelocityAdvection3D) =
    ("phi", "vx", "vy", "vz")
@inline Trixi.have_constant_speed(::FrozenVelocityAdvection3D) = Trixi.False()
@inline Trixi.max_abs_speeds(u, ::FrozenVelocityAdvection3D) = (abs(u[2]), abs(u[3]), abs(u[4]))
@inline Trixi.max_abs_speed_naive(ul, ur, d::Integer, ::FrozenVelocityAdvection3D) =
    max(abs(ul[d+1]), abs(ur[d+1]))
@inline Trixi.max_abs_speed_naive(ul, ur, normal::AbstractVector, ::FrozenVelocityAdvection3D) =
    max(abs(ul[2]*normal[1] + ul[3]*normal[2] + ul[4]*normal[3]),
        abs(ur[2]*normal[1] + ur[3]*normal[2] + ur[4]*normal[3]))
@inline Trixi.max_abs_speed(ul, ur, d::Integer, eq::FrozenVelocityAdvection3D) =
    Trixi.max_abs_speed_naive(ul, ur, d, eq)
@inline Trixi.max_abs_speed(ul, ur, normal::AbstractVector, eq::FrozenVelocityAdvection3D) =
    Trixi.max_abs_speed_naive(ul, ur, normal, eq)
@inline function Trixi.flux(u, d::Integer, ::FrozenVelocityAdvection3D)
    d == 1 && return SVector(u[2]*u[1], 0.0, 0.0, 0.0)
    d == 2 && return SVector(u[3]*u[1], 0.0, 0.0, 0.0)
    return SVector(u[4]*u[1], 0.0, 0.0, 0.0)
end
@inline function Trixi.flux(u, normal::AbstractVector, ::FrozenVelocityAdvection3D)
    a = u[2]*normal[1] + u[3]*normal[2] + u[4]*normal[3]
    return SVector(a*u[1], 0.0, 0.0, 0.0)
end
@inline Trixi.cons2prim(u, ::FrozenVelocityAdvection3D) = u
@inline Trixi.cons2entropy(u, ::FrozenVelocityAdvection3D) = u

# Only phi receives numerical dissipation; frozen velocity components stay exact.
@inline function _scalar_lax_friedrichs(ul, ur, direction, eq::FrozenVelocityAdvection3D)
    fl, fr = Trixi.flux(ul, direction, eq), Trixi.flux(ur, direction, eq)
    λ = Trixi.max_abs_speed(ul, ur, direction, eq)
    return SVector(0.5*(fl[1] + fr[1]) - 0.5λ*(ur[1] - ul[1]), 0.0, 0.0, 0.0)
end

_velocity(x, y, z) = (π*y, -π*x, 0.0)

function _initial_state(x, t, ::FrozenVelocityAdvection3D, initial_function)
    vx, vy, vz = _velocity(x...)
    return SVector(initial_function(x...), vx, vy, vz)
end

function _lagrange_weights(nodes, ξ)
    n = length(nodes)
    return [prod((ξ - nodes[j])/(nodes[i] - nodes[j]) for j in 1:n if j != i) for i in 1:n]
end

# At element faces use the lowest-index Cartesian element trace.
function _dg_value(u, nodes, cells, limits, a, b, c)
    lo, hi = limits
    (lo <= a <= hi && lo <= b <= hi && lo <= c <= hi) ||
        throw(ArgumentError("DG-to-grid sample lies outside the Trixi domain"))
    h = (hi - lo)/cells
    ix = clamp(ceil(Int, (a - lo)/h), 1, cells)
    iy = clamp(ceil(Int, (b - lo)/h), 1, cells)
    iz = clamp(ceil(Int, (c - lo)/h), 1, cells)
    e = ix + (iy - 1)*cells + (iz - 1)*cells^2
    ξ = 2*(a - (lo + (ix - 1)*h))/h - 1
    η = 2*(b - (lo + (iy - 1)*h))/h - 1
    ζ = 2*(c - (lo + (iz - 1)*h))/h - 1
    wx, wy, wz = _lagrange_weights(nodes, ξ), _lagrange_weights(nodes, η),
                 _lagrange_weights(nodes, ζ)
    val = zero(eltype(u))
    for k in eachindex(wz), j in eachindex(wy), i in eachindex(wx)
        val += wx[i]*wy[j]*wz[k]*u[1, i, j, k, e]
    end
    return val
end

"""
    run_trixi3d(grid, initial_function, times; cells=16, polydeg=3, cfl=0.25)

Evolve the canonical 3D Zalesak sphere with Trixi's structured DGSEM on
`[-1,1]^3`, velocity `π(y,-x,0)`, frozen velocity components and scalar-only
Lax–Friedrichs dissipation. `boundary_condition_do_nothing` uses interior-trace
extrapolation, not a periodic boundary. Returns x-first snapshots on `grid`.
"""
function run_trixi3d(output_grid, initial_function, times;
                     cells::Integer=16, polydeg::Integer=3, cfl::Real=0.25)
    length(output_grid.axes) == 3 || throw(ArgumentError("output grid must be 3D"))
    cells > 0 && polydeg >= 1 || throw(ArgumentError("cells must be positive and polydeg ≥ 1"))
    isfinite(cfl) && 0 < cfl <= 0.4 || throw(ArgumentError("cfl must lie in (0,0.4]"))
    ts = Float64.(collect(times))
    !isempty(ts) && all(isfinite, ts) && issorted(ts) && first(ts) == 0 && last(ts) > 0 ||
        throw(ArgumentError("times must be finite, increasing from 0 to a positive final time"))
    length(unique(ts)) == length(ts) || throw(ArgumentError("times must be strictly increasing"))
    for axis in output_grid.axes
        first(axis) >= -1 && last(axis) <= 1 ||
            throw(ArgumentError("the Trixi Zalesak domain is [-1,1]^3"))
    end

    lims = (-1.0, 1.0)
    equations = FrozenVelocityAdvection3D()
    initial_condition(x, t, eq) = _initial_state(x, t, eq, initial_function)
    mesh = Trixi.StructuredMesh((cells, cells, cells), (lims[1], lims[1], lims[1]),
        (lims[2], lims[2], lims[2]), periodicity=false)
    solver = Trixi.DGSEM(polydeg=polydeg, surface_flux=_scalar_lax_friedrichs)
    semi = Trixi.SemidiscretizationHyperbolic(mesh, equations, initial_condition, solver;
        boundary_conditions=Trixi.boundary_condition_do_nothing)
    ode = Trixi.semidiscretize(semi, (0.0, last(ts)))
    dt = cfl*(lims[2] - lims[1])/cells/((2polydeg + 1)*(2π))
    integrator = Trixi.init(ode, Trixi.SimpleSSPRK33(); dt)
    nodes = solver.basis.nodes
    states = [reshape(copy(integrator.u), 4, polydeg+1, polydeg+1, polydeg+1, cells^3)]
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
        push!(states, reshape(copy(integrator.u), 4, polydeg+1, polydeg+1,
            polydeg+1, cells^3))
    end
    x, y, z = output_grid.axes
    snapshots = [[_dg_value(u, nodes, cells, lims, a, b, c)
        for a in x, b in y, c in z] for u in states]
    integrals = [Trixi.integrate((v, _) -> SVector(v[1]), vec(u), semi; normalize=false)[1]
        for u in states]
    return (; x=collect(x), y=collect(y), z=collect(z), times=actual_times,
        snapshots, elapsed, dt, steps, cells, polydeg, cfl, integral=integrals,
        metadata=(equation=:canonical_trixi_zalesak_3d,
            boundary=:interior_extrapolation, state_variables=(:phi, :vx, :vy, :vz),
            scalar_dofs=(polydeg+1)^3*cells^3, integrator=:Trixi_SimpleSSPRK33,
            dt_is_maximum=true, native_integral=:quadrature_integral_of_phi_not_phase_volume))
end

"""Minimal exchange contract: a constant survives the 3D DG solve."""
function self_check()
    out = (axes=(collect(-1.0:0.5:1.0), collect(-1.0:0.5:1.0), collect(-1.0:0.5:1.0)),)
    constant = run_trixi3d(out, (x, y, z) -> 2.0, [0.0, 1e-5]; cells=4, polydeg=3, cfl=0.2)
    maximum(abs.(constant.snapshots[end] .- 2.0)) < 1e-10 ||
        error("3D DG did not preserve a constant field")
    true
end

end # module
