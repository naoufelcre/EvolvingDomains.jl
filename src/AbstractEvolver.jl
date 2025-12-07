# =============================================================================
# Abstract Interface for Level Set Evolvers
# =============================================================================
# This defines the contract that any backend must implement.
# Replace LevelSetMethods.jl by implementing these methods for your custom solver.

"""
    AbstractLevelSetEvolver

Abstract type for level set evolution backends.

Any backend must implement:
- `evolve!(evolver, Δt)`: Advance the level set by time Δt
- `current_values(evolver)`: Return the current nodal level set values
- `current_time(evolver)`: Return the current time
- `reinitialize!(evolver)`: Restore signed distance property
- `grid_coords(evolver)`: Return the grid node coordinates
"""
abstract type AbstractLevelSetEvolver end

# =============================================================================
# Required Interface Methods
# =============================================================================

"""
    evolve!(evolver::AbstractLevelSetEvolver, Δt::Real)

Advance the level set by time Δt. Modifies `evolver` in-place.
"""
function evolve!(evolver::AbstractLevelSetEvolver, Δt::Real)
    error("evolve! not implemented for $(typeof(evolver))")
end

"""
    current_values(evolver::AbstractLevelSetEvolver) -> Vector{Float64}

Return the current nodal values of the level set function.
"""
function current_values(evolver::AbstractLevelSetEvolver)
    error("current_values not implemented for $(typeof(evolver))")
end

"""
    current_time(evolver::AbstractLevelSetEvolver) -> Float64

Return the current simulation time.
"""
function current_time(evolver::AbstractLevelSetEvolver)
    error("current_time not implemented for $(typeof(evolver))")
end

"""
    reinitialize!(evolver::AbstractLevelSetEvolver)

Restore the signed distance property of the level set.
"""
function reinitialize!(evolver::AbstractLevelSetEvolver)
    error("reinitialize! not implemented for $(typeof(evolver))")
end

"""
    grid_coords(evolver::AbstractLevelSetEvolver) -> Vector{VectorValue}

Return the grid node coordinates as a vector of Gridap VectorValues.
"""
function grid_coords(evolver::AbstractLevelSetEvolver)
    error("grid_coords not implemented for $(typeof(evolver))")
end

# =============================================================================
# Optional Interface Methods (with defaults)
# =============================================================================

"""
    update_velocity!(evolver::AbstractLevelSetEvolver, velocity_source, t)

Update the evolver's internal velocity representation from an external source.
Called before each time step when using time-dependent or FE-coupled velocities.

Default implementation does nothing (for static velocity cases).
Override this for evolvers that support velocity updates.
"""
function update_velocity!(evolver::AbstractLevelSetEvolver, velocity_source, t)
    # Default: no-op for static velocity evolvers
    return evolver
end

"""
    supports_velocity_update(evolver::AbstractLevelSetEvolver) -> Bool

Return whether this evolver supports runtime velocity updates.
"""
supports_velocity_update(::AbstractLevelSetEvolver) = false
