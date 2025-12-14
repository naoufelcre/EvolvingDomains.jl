# =============================================================================
# EvolvingDiscreteGeometry — Backend-Agnostic Evolving Domain
# =============================================================================

"""
    EvolvingDiscreteGeometry{E<:AbstractLevelSetEvolver}

A wrapper that couples a level set evolver with a GridapEmbedded DiscreteGeometry.

# Fields
- `evolver`: The backend that evolves the level set (implements AbstractLevelSetEvolver)
- `bg_model`: The background CartesianDiscreteModel
- `geometry`: The current DiscreteGeometry (lazily regenerated when dirty)
- `geometry_dirty`: Whether geometry needs rebuild (set after advance!)
- `cached_cut`: Cached EmbeddedDiscretization (lazily computed)

# Performance Features
- Lazy geometry rebuild: only rebuilds when accessed AND dirty
- Cached cut geometry: avoids redundant cut() calls within same time step

# Example
```julia
eg = EvolvingDiscreteGeometry(evolver, model)
advance!(eg, 0.01)
geo = current_geometry(eg)   # Rebuilds geometry (lazy)
cut_geo = current_cut(eg)    # Computes cut, caches result
cut_geo2 = current_cut(eg)   # Returns cached, no recomputation
```
"""
mutable struct EvolvingDiscreteGeometry{E<:AbstractLevelSetEvolver}
    evolver::E
    bg_model::CartesianDiscreteModel
    geometry::GridapEmbedded.LevelSetCutters.DiscreteGeometry
    geometry_dirty::Bool
    cached_cut::Union{Nothing, GridapEmbedded.Interfaces.EmbeddedDiscretization}
    # Curvature caching
    cached_curvature::Union{Nothing, Gridap.FESpaces.SingleFieldFEFunction}
    curvature_dirty::Bool
end

"""
    EvolvingDiscreteGeometry(evolver, bg_model)

Construct an EvolvingDiscreteGeometry from an evolver and background model.
"""
function EvolvingDiscreteGeometry(evolver::AbstractLevelSetEvolver, bg_model::CartesianDiscreteModel)
    geo = _rebuild_geometry(evolver, bg_model)
    return EvolvingDiscreteGeometry(evolver, bg_model, geo, false, nothing, nothing, true)
end

# =============================================================================
# Main API
# =============================================================================

"""
    advance!(eg::EvolvingDiscreteGeometry, Δt::Real; velocity=nothing, lazy=true)

Advance the geometry by time Δt. This evolves the level set.

# Arguments
- `Δt`: Time step size
- `velocity`: Optional velocity source to use for this advance
- `lazy`: If true (default), defer geometry rebuild until accessed

# Performance
With `lazy=true`, geometry and cut are only rebuilt when accessed via
`current_geometry()` or `current_cut()`. This avoids redundant rebuilds
if you access the same geometry multiple times.

# Example
```julia
# Using configured velocity (from evolver construction)
advance!(eg, 0.01)

# With FE-coupled velocity
velocity_fh = solve_stokes(current_geometry(eg))
advance!(eg, 0.01; velocity=FEVelocitySource(velocity_fh, model))
```
"""
function advance!(eg::EvolvingDiscreteGeometry, Δt::Real; velocity=nothing, lazy=true)
    t = current_time(eg.evolver)
    
    # Update velocity if provided
    if !isnothing(velocity)
        update_velocity!(eg.evolver, velocity, t)
    end
    
    evolve!(eg.evolver, Δt)
    
    # Mark geometry as dirty (needs rebuild)
    eg.geometry_dirty = true
    eg.cached_cut = nothing  # Invalidate cut cache
    eg.curvature_dirty = true  # Invalidate curvature cache
    
    # Optionally rebuild immediately (non-lazy mode)
    if !lazy
        _ensure_geometry_fresh!(eg)
    end
    
    return eg
end

"""
    current_geometry(eg::EvolvingDiscreteGeometry) -> DiscreteGeometry

Return the current DiscreteGeometry (for use with `cut`).
Rebuilds geometry lazily if dirty.
"""
function current_geometry(eg::EvolvingDiscreteGeometry)
    _ensure_geometry_fresh!(eg)
    return eg.geometry
end

"""
    current_cut(eg::EvolvingDiscreteGeometry) -> EmbeddedDiscretization

Return the current cut geometry, with caching.
Computes `cut(bg_model, geometry)` only if not already cached.

This is more efficient than calling `cut(model, current_geometry(eg))`
multiple times within the same time step.
"""
function current_cut(eg::EvolvingDiscreteGeometry)
    _ensure_geometry_fresh!(eg)
    if isnothing(eg.cached_cut)
        eg.cached_cut = cut(eg.bg_model, eg.geometry)
    end
    return eg.cached_cut
end

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

Restore the signed distance property. Marks geometry as dirty.
"""
function reinitialize!(eg::EvolvingDiscreteGeometry)
    reinitialize!(eg.evolver)
    eg.geometry_dirty = true
    eg.cached_cut = nothing
    eg.curvature_dirty = true
    return eg
end

"""
    invalidate_cache!(eg::EvolvingDiscreteGeometry)

Force invalidation of cached geometry and cut. Useful when external
changes affect the geometry (e.g., manual level set modification).
"""
function invalidate_cache!(eg::EvolvingDiscreteGeometry)
    eg.geometry_dirty = true
    eg.cached_cut = nothing
    eg.cached_curvature = nothing
    eg.curvature_dirty = true
    return eg
end

# =============================================================================
# External Solver Interface
# =============================================================================

"""
    grid_info(eg::EvolvingDiscreteGeometry) -> CartesianGridInfo

Get structured grid metadata for external solvers.

# Example
```julia
info = grid_info(eg)
Δx, Δy = info.spacing
nx, ny = info.dims
```
"""
grid_info(eg::EvolvingDiscreteGeometry) = CartesianGridInfo(eg.bg_model)

"""
    set_levelset!(eg::EvolvingDiscreteGeometry, ϕ_new::Vector{Float64})

Update level set from external solver. Marks geometry as dirty.

This is the primary entry point for external Cartesian-grid hyperbolic
solvers to inject their computed level set values.

# Arguments
- `eg`: The evolving geometry
- `ϕ_new`: New level set values (must match node count)

# Example
```julia
# External solver evolves the level set
ϕ = current_levelset(eg)
ϕ_new = external_hyperbolic_step(grid_info(eg), ϕ, u_data, Δt)

# Inject back and optionally reinitialize
set_levelset!(eg, ϕ_new)
reinitialize!(eg)
```
"""
function set_levelset!(eg::EvolvingDiscreteGeometry, ϕ_new::Vector{Float64})
    set_values!(eg.evolver, ϕ_new)
    eg.geometry_dirty = true
    eg.cached_cut = nothing
    eg.curvature_dirty = true
    return eg
end

# =============================================================================
# Internal Helpers
# =============================================================================

"""
    _ensure_geometry_fresh!(eg::EvolvingDiscreteGeometry)

Rebuild geometry if dirty. Internal helper for lazy evaluation.
"""
function _ensure_geometry_fresh!(eg::EvolvingDiscreteGeometry)
    if eg.geometry_dirty
        eg.geometry = _rebuild_geometry(eg.evolver, eg.bg_model)
        eg.geometry_dirty = false
        # Note: don't clear cached_cut here, it's handled separately
    end
    return nothing
end

"""
    _rebuild_geometry(evolver, bg_model) -> DiscreteGeometry

Construct a DiscreteGeometry from the current evolver state.
"""
function _rebuild_geometry(evolver::AbstractLevelSetEvolver, bg_model::CartesianDiscreteModel)
    vals = current_values(evolver)
    coords = grid_coords(evolver)
    return GridapEmbedded.LevelSetCutters.DiscreteGeometry(vals, coords, name="evolving_ls")
end

