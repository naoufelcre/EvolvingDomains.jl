# Prescribed-motion smoke example: an explicit simulation loop that moves a level-set
# geometry and a transported nodal state together, using only the exported API.
#
# This exercises the workflow, not its accuracy: a coarse grid, a few steps, and the only
# balance checked is the unweighted nodal sum of the transport map's support. No physical
# mass or area conservation is claimed.
#
# From the root of an instantiated EvolvingDomains checkout:
# julia --startup-file=no --project=. skills/evolvingdomains/examples/prescribed_motion.jl

using EvolvingDomains: Grid, Discretization, Field, cut, eulerianTransport!, WENO5Cache,
    TransportMap, lagrangianTransport!, cutCellIds
using StaticArrays: SVector

function prescribed_motion()
    # Dimensionless toy: 33×33 nodes and five steps keep the check small; motion is 0.025.
    nsteps = 5
    # A small step for this short smoke case, not a general stability prescription.
    dt = 0.01
    axis = collect(range(-1.0, 1.0; length = 33))
    grid = Grid((axis, copy(axis)))
    d = Discretization(Dict(:init => grid))

    # Application-defined prescribed dynamics: uniform translation to the right. Nothing
    # here derives velocity from an equation; the example supplies it.
    velocity_of(grid::Grid) =
        fill(SVector(0.5, 0.0), length(grid.axes[1]) * length(grid.axes[2]))

    # Separate fields on the same grid: `levelset` carries the moving geometry, `state` is
    # the transported quantity, and `scratch` is its distinct target buffer.
    radius = 0.35                     # compact support of the state, well inside r = 0.5
    function state_bump(x, y)
        r2 = x^2 + y^2
        return r2 < radius^2 ? (1 - r2 / radius^2)^2 : 0.0
    end
    levelset = Field(d, Dict(:init => [hypot(x, y) - 0.5 for x in axis, y in axis]))
    state = Field(d, Dict(:init => [state_bump(x, y) for x in axis, y in axis]))
    scratch = Field(d, Dict(:init => zeros(length(axis), length(axis))))
    # The state's nonzeros all lie in the in-phase region, so the masked sum below covers
    # the whole support instead of hiding the part outside it.
    inside = findall(!iszero, vec(state.data[:init]))
    @assert all(φ -> φ < 0, vec(levelset.data[:init])[inside]) "State support must lie in the negative phase (φ < 0)"

    data = (state.data[:init], scratch.data[:init], levelset.data[:init])
    @assert !Base.mightalias(data[1], data[2])
    @assert !Base.mightalias(data[1], data[3])
    @assert !Base.mightalias(data[2], data[3])
    cache = WENO5Cache()
    previous = cut(levelset, :init)
    t = 0.0
    for step in 1:nsteps
        velocity = velocity_of(grid)
        old_points = copy(previous.boundary.physical_points)

        # Eulerian advection updates the level set in place; the cut must then be rebuilt
        # so the transport map knows the geometry at t and at t + dt.
        eulerianTransport!(levelset, velocity, dt; level_id = :init, cache)
        current = cut(levelset, :init)
        transport = TransportMap(previous, current, grid, velocity, dt)

        # This short interior-support run retains contributions on the target support.
        # Check the supported unweighted nodal sum, not physical mass.
        expected = sum(vec(state.data[:init])[transport.source_mask])
        lagrangianTransport!(scratch, state, transport)
        state, scratch = scratch, state     # state holds t + dt; scratch is free again
        t += dt

        @assert all(isfinite, levelset.data[:init])
        @assert all(isfinite, state.data[:init])
        @assert !isempty(cutCellIds(current))     # the geometry is genuinely cut each step
        @assert previous.boundary.physical_points == old_points
        @assert current.boundary.physical_points != old_points
        # These tolerances check Float64 bookkeeping, not spatial or temporal accuracy.
        @assert isapprox(t, step * dt; rtol = 1e-12)
        @assert isapprox(sum(state.data[:init]), expected; rtol = 1e-10) "Step $step: expected supported nodal sum $expected, got $(sum(state.data[:init]))"

        previous = current
    end
    println("Prescribed-motion smoke checks passed: steps=$nsteps, time=$t, unweighted nodal sum=$(sum(state.data[:init]))")
    return t
end

prescribed_motion()
