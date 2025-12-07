# =============================================================================
# EvolvingDiscreteGeometry — Backend-Agnostic Evolving Domain
# =============================================================================

"""
    EvolvingDiscreteGeometry{E<:AbstractLevelSetEvolver}

A wrapper that couples a level set evolver with a GridapEmbedded DiscreteGeometry.

# Fields
- `evolver`: The backend that evolves the level set (implements AbstractLevelSetEvolver)
- `bg_model`: The background CartesianDiscreteModel
- `geometry`: The current DiscreteGeometry (regenerated after each advance)

# Example
```julia
eg = EvolvingDiscreteGeometry(evolver, model)
advance!(eg, 0.01)
geo = current_geometry(eg)
cut_geo = cut(model, geo)
```
"""
mutable struct EvolvingDiscreteGeometry{E<:AbstractLevelSetEvolver}
    evolver::E
    bg_model::CartesianDiscreteModel
    geometry::GridapEmbedded.LevelSetCutters.DiscreteGeometry
end

"""
    EvolvingDiscreteGeometry(evolver, bg_model)

Construct an EvolvingDiscreteGeometry from an evolver and background model.
"""
function EvolvingDiscreteGeometry(evolver::AbstractLevelSetEvolver, bg_model::CartesianDiscreteModel)
    geo = _rebuild_geometry(evolver, bg_model)
    return EvolvingDiscreteGeometry(evolver, bg_model, geo)
end

# =============================================================================
# Main API
# =============================================================================

"""
    advance!(eg::EvolvingDiscreteGeometry, Δt::Real; velocity=nothing)

Advance the geometry by time Δt. This evolves the level set and regenerates
the DiscreteGeometry.

# Arguments
- `Δt`: Time step size
- `velocity`: Optional velocity source to use for this advance. Can be:
  - `nothing`: Use the velocity configured at evolver construction
  - `AbstractVelocitySource`: Update evolver velocity before advancing
  - `CellField/FEFunction`: Automatically wrapped in FEVelocitySource

# Example
```julia
# Using configured velocity (from evolver construction)
advance!(eg, 0.01)

# With FE-coupled velocity
velocity_fh = solve_stokes(current_geometry(eg))
advance!(eg, 0.01; velocity=FEVelocitySource(velocity_fh, model))

# Or directly (auto-wrapped if evolver has FEVelocitySource)
advance!(eg, 0.01; velocity=velocity_fh)
```
"""
function advance!(eg::EvolvingDiscreteGeometry, Δt::Real; velocity=nothing)
    t = current_time(eg.evolver)
    
    # Update velocity if provided
    if !isnothing(velocity)
        update_velocity!(eg.evolver, velocity, t)
    end
    
    evolve!(eg.evolver, Δt)
    eg.geometry = _rebuild_geometry(eg.evolver, eg.bg_model)
    return eg
end

"""
    current_geometry(eg::EvolvingDiscreteGeometry) -> DiscreteGeometry

Return the current DiscreteGeometry (for use with `cut`).
"""
current_geometry(eg::EvolvingDiscreteGeometry) = eg.geometry

"""
    current_time(eg::EvolvingDiscreteGeometry) -> Float64

Return the current simulation time.
"""
current_time(eg::EvolvingDiscreteGeometry) = current_time(eg.evolver)

"""
    current_levelset(eg::EvolvingDiscreteGeometry) -> Vector{Float64}

Return the current nodal level set values.
"""
current_levelset(eg::EvolvingDiscreteGeometry) = current_values(eg.evolver)

"""
    reinitialize!(eg::EvolvingDiscreteGeometry)

Restore the signed distance property and regenerate the geometry.
"""
function reinitialize!(eg::EvolvingDiscreteGeometry)
    reinitialize!(eg.evolver)
    eg.geometry = _rebuild_geometry(eg.evolver, eg.bg_model)
    return eg
end

# =============================================================================
# Internal Helpers
# =============================================================================

"""
    _rebuild_geometry(evolver, bg_model) -> DiscreteGeometry

Construct a DiscreteGeometry from the current evolver state.
"""
function _rebuild_geometry(evolver::AbstractLevelSetEvolver, bg_model::CartesianDiscreteModel)
    vals = current_values(evolver)
    coords = grid_coords(evolver)
    return GridapEmbedded.LevelSetCutters.DiscreteGeometry(vals, coords, name="evolving_ls")
end
