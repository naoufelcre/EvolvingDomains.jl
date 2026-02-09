module Kinematic

using StaticArrays: SVector
using Gridap.TensorValues: VectorValue
using ..Geometric: EvolvingDiscreteGeometry, grid_info, set_levelset!

# Include the implementation
include("VelocitySource.jl")
include("TransportOperators/WENO5.jl")
include("TransportOperators/SemiLagrangian.jl")

using .WENO5
using .SemiLagrangian

export AbstractVelocitySource
export StaticFunctionVelocity, TimeDependentVelocity
export sample_velocity, is_time_dependent
export advance!, weno5_step!

# Export SemiLagrangian tools
export compute_raymaps, RayMaps

"""
    advance!(geom::EvolvingDiscreteGeometry, Δt) -> geom

Evolve the geometry by time Δt using WENO5 transport.
Updates the geometry via `set_levelset!` to ensure history is preserved.
"""
function advance!(geom::EvolvingDiscreteGeometry, velocity_field, Δt::Real)
    # Get components
    phi_old = geom.levelset
    grid = grid_info(geom.grid)

    # Create a working copy for the new state
    # weno5_step! modifies its first argument in place
    phi_new = copy(phi_old)

    # Perform transport step (Euler + WENO5)
    weno5_step!(phi_new, grid, velocity_field, Float64(Δt))

    # Update geometry (triggers cache management)
    set_levelset!(geom, phi_new)

    return geom
end

end # module
