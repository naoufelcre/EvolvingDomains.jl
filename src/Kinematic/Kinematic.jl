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
export TransportMap, advect!

# NOTE on naming: `advance!` evolves the geometry's level set (moves the interface).
# `advect!` transports a scalar/vector field on a fixed TransportMap.
# These are distinct operations; the naming distinction is intentional.

"""
    advance!(geom::EvolvingDiscreteGeometry, velocity_field, Δt) -> geom

Evolve the geometry by time Δt using WENO5 + SSP-RK3 transport.
Updates the geometry via `set_levelset!` to ensure cache history is preserved.

A persistent scratch buffer (`geom.cache.weno_scratch`) is allocated on the first call
and reused on subsequent calls, eliminating one `Vector{Float64}` allocation per step.
The WENO5Cache (RK3 stage buffers) is stored in `geom.cache.weno_cache` and is
owned by the geometry object, ensuring thread-safety for multi-geometry workflows.
"""
function advance!(geom::EvolvingDiscreteGeometry, velocity_field, Δt::Real)
    grid = grid_info(geom.grid)
    n = length(geom.levelset)

    # Lazy-allocate a persistent scratch buffer for the new-state levelset.
    # This avoids copy(phi_old) on every call.
    if isnothing(geom.cache.weno_scratch) || length(geom.cache.weno_scratch) != n
        geom.cache.weno_scratch = Vector{Float64}(undef, n)
    end
    phi_new = geom.cache.weno_scratch

    # Copy current levelset into scratch buffer (weno5_step! modifies in-place)
    copyto!(phi_new, geom.levelset)

    # Perform WENO5 + SSP-RK3 transport step
    # Uses the geometry-owned weno_cache for RK3 stage buffers.
    weno5_step!(phi_new, grid, velocity_field, Float64(Δt), geom.cache.weno_cache)

    # Update geometry (triggers cache invalidation and history shift)
    set_levelset!(geom, phi_new)

    return geom
end

end # module
