# =============================================================================
# EvolvingGeometry.jl
# =============================================================================
#
# This module provides the main EvolvingDiscreteGeometry type which combines:
# - Level set evolution via LevelSetMethods.jl
# - Integration with GridapEmbedded for CutFEM
# - Simple API: advance!, current_cut, set_velocity!
#
# =============================================================================

# Note: All imports come from the parent module (EvolvingDomains.jl)
# LevelSetMethods is now a direct dependency, not accessed via Main

using StaticArrays: SVector
using GridapEmbedded
import LevelSetMethods as LSM

"""
    EvolvingDiscreteGeometry

Main container for evolving domain simulations.

# Example
```julia
model = CartesianDiscreteModel((-2,2,-2,2), (50,50))
geom = EvolvingDiscreteGeometry(model, x -> norm(x) - 1.0)
set_velocity!(geom, TimeDependentVelocity((x,t) -> (-x[2], x[1])))

for step in 1:100
    advance!(geom, 0.01)
    cut = current_cut(geom)  # For CutFEM
end
```
"""
mutable struct EvolvingDiscreteGeometry{Eq}
    equation::Eq                           # LevelSetEquation from LevelSetMethods
    bg_model::DiscreteModel                # Background Gridap model (Unstructured/Quadtree/Cartesian)
    lsm_model::CartesianDiscreteModel      # Cartesian Gridap model for Level Set Evolution
    velocity::Union{Nothing, AbstractVelocitySource}
    velocity_buffer::Vector{SVector{2,Float64}}  # Cached velocity for update_func
    coords_tuples::Vector{NTuple{2,Float64}}     # Node coords for velocity sampling (on lsm_model)

    # Cached derived quantities (lazy)
    cached_cut::Union{Nothing, GridapEmbedded.Interfaces.EmbeddedDiscretization}
    dirty::Bool

    # Time tracking
    t::Float64
    step::Int

    # Reinitialization config
    reinit_freq::Int
    last_reinit_step::Int
end

"""
    EvolvingDiscreteGeometry(model, initial_sdf; velocity=nothing, reinit_freq=10)

Construct an evolving geometry using the same Cartesian model for both evolution and CutFEM.
"""
function EvolvingDiscreteGeometry(model::CartesianDiscreteModel, initial_sdf::Function;
                                   velocity::Union{Nothing, AbstractVelocitySource}=nothing,
                                   reinit_freq::Int=10)
    return EvolvingDiscreteGeometry(model, model, initial_sdf; velocity=velocity, reinit_freq=reinit_freq)
end

"""
    EvolvingDiscreteGeometry(bg_model, lsm_model, initial_sdf; velocity=nothing, reinit_freq=10)

Construct an evolving geometry with decoupled models:
- `bg_model`: The background model for CutFEM (can be unstructured, Quadtree, etc.)
- `lsm_model`: The Cartesian model for Level Set evolution (MUST be Cartesian)
"""
function EvolvingDiscreteGeometry(bg_model::DiscreteModel, lsm_model::CartesianDiscreteModel, initial_sdf::Function;
                                   velocity::Union{Nothing, AbstractVelocitySource}=nothing,
                                   reinit_freq::Int=10)
    info = grid_info(lsm_model)

    # === Input Validation ===
    length(info.dims) == 2 || error("Only 2D grids are supported (got $(length(info.dims))D)")
    
    # Test SDF at origin to catch NaN/Inf early
    test_point = info.origin
    test_val = initial_sdf(test_point)
    isfinite(test_val) || error("initial_sdf must return finite values (got $test_val at $test_point)")

    # Create LevelSetMethods grid
    origin = info.origin
    corner = info.origin .+ info.spacing .* info.cells
    partition = info.dims  # Node counts

    grid = LSM.CartesianGrid(origin, corner, partition)

    # Create initial level set
    ϕ = LSM.LevelSet(initial_sdf, grid)

    # Get node coordinates for velocity sampling
    # LevelSetMethods grid coords are iterable
    coords_tuples = [NTuple{2,Float64}(Tuple(c)) for c in vec(collect(grid))]

    # Sample initial velocity
    n_nodes = prod(partition)
    vel_buffer = if velocity !== nothing
        [SVector{2,Float64}(sample_velocity(velocity, c, 0.0)...) for c in coords_tuples]
    else
        fill(SVector{2,Float64}(0.0, 0.0), n_nodes)
    end

    # Create velocity MeshField
    vel_array = reshape(vel_buffer, Tuple(partition))
    u_mesh = LSM.MeshField(vel_array, grid, nothing)

    # Create update function for dynamic velocities
    update_fn = if velocity !== nothing && is_time_dependent(velocity)
        (u, ϕ_current, t) -> begin
            for (i, c) in enumerate(coords_tuples)
                v = sample_velocity(velocity, c, t)
                vel_buffer[i] = SVector{2,Float64}(v[1], v[2])
            end
            nothing
        end
    else
        (u, ϕ_current, t) -> nothing
    end

    # Create advection term
    terms = (LSM.AdvectionTerm(u_mesh, LSM.WENO5(), update_fn),)

    # Create level set equation
    eq = LSM.LevelSetEquation(;
        terms = terms,
        levelset = ϕ,
        integrator = LSM.RK3(),
        bc = LSM.PeriodicBC(),
        t = 0.0
    )

    return EvolvingDiscreteGeometry(
        eq, bg_model, lsm_model, velocity, vel_buffer, coords_tuples,
        nothing, true,  # cached_cut, dirty
        0.0, 0,         # t, step
        reinit_freq, 0  # reinit config
    )
end

"""
    advance!(geom::EvolvingDiscreteGeometry, Δt; lazy=true) -> geom

Evolve the geometry by time Δt.
"""
function advance!(geom::EvolvingDiscreteGeometry, Δt::Real; lazy::Bool=true)
    # Integrate to new time
    tf = LSM.current_time(geom.equation) + Δt
    LSM.integrate!(geom.equation, tf)

    # Update state
    geom.t = tf
    geom.step += 1
    geom.dirty = true
    geom.cached_cut = nothing

    # Check reinitialization
    if geom.reinit_freq > 0 && (geom.step - geom.last_reinit_step) >= geom.reinit_freq
        reinitialize!(geom)
    end

    if !lazy
        _ensure_fresh!(geom)
    end

    return geom
end

"""
    current_geometry(geom::EvolvingDiscreteGeometry) -> DiscreteGeometry
"""
function current_geometry(geom::EvolvingDiscreteGeometry)
    ϕ = current_levelset(geom)
    return _build_discrete_geometry(ϕ, geom.bg_model, geom.lsm_model)
end

"""
    current_cut(geom::EvolvingDiscreteGeometry) -> EmbeddedDiscretization

Return the current cut discretization (lazy rebuild).
"""
function current_cut(geom::EvolvingDiscreteGeometry)
    _ensure_fresh!(geom)
    return geom.cached_cut
end

"""
    current_levelset(geom::EvolvingDiscreteGeometry) -> Vector{Float64}
"""
function current_levelset(geom::EvolvingDiscreteGeometry)
    state = LSM.current_state(geom.equation)
    return vec(LSM.values(state))
end

"""
    current_time(geom::EvolvingDiscreteGeometry) -> Float64
"""
current_time(geom::EvolvingDiscreteGeometry) = geom.t

"""
    set_levelset!(geom::EvolvingDiscreteGeometry, ϕ_new::Vector{Float64})
"""
function set_levelset!(geom::EvolvingDiscreteGeometry, ϕ_new::AbstractVector{Float64})
    state = LSM.current_state(geom.equation)
    vals = LSM.values(state)
    copyto!(vec(vals), ϕ_new)
    geom.dirty = true
    geom.cached_cut = nothing
    return geom
end

"""
    reinitialize!(geom::EvolvingDiscreteGeometry)
"""
function reinitialize!(geom::EvolvingDiscreteGeometry)
    # Fix: reinitialize! is defined for LevelSet, stored in .state
    LSM.reinitialize!(geom.equation.state) 
    geom.last_reinit_step = geom.step
    geom.dirty = true
    geom.cached_cut = nothing
    return geom
end


"""
    set_velocity!(geom::EvolvingDiscreteGeometry, vel::AbstractVelocitySource)
"""
function set_velocity!(geom::EvolvingDiscreteGeometry, vel::AbstractVelocitySource)
    geom.velocity = vel
    # Eagerly update velocity buffer for next advance
    for (i, c) in enumerate(geom.coords_tuples)
        v = sample_velocity(vel, c, geom.t)
        geom.velocity_buffer[i] = SVector{2,Float64}(v[1], v[2])
    end
    return geom
end

"""
    grid_info(geom::EvolvingDiscreteGeometry) -> CartesianGridInfo
"""
grid_info(geom::EvolvingDiscreteGeometry) = grid_info(geom.lsm_model)

# =============================================================================
# Internal Helpers
# =============================================================================

function _ensure_fresh!(geom::EvolvingDiscreteGeometry)
    if geom.dirty || geom.cached_cut === nothing
        geo = current_geometry(geom)
        geom.cached_cut = cut(geom.bg_model, geo)
        geom.dirty = false
    end
    return nothing
end

function _build_discrete_geometry(ϕ::Vector{Float64}, bg_model::DiscreteModel, lsm_model::CartesianDiscreteModel)
    
    # 1. Create FE Function on the LSM grid
    # If models are identical, we can skip interpolation for efficiency
    if bg_model === lsm_model
        reffe_lsm = ReferenceFE(lagrangian, Float64, 1)
        V_lsm = FESpace(lsm_model, reffe_lsm)
        ϕ_lsm_fe = FEFunction(V_lsm, ϕ)
        return GridapEmbedded.LevelSetCutters.DiscreteGeometry(ϕ_lsm_fe, bg_model)
    end

    # 2. If background model is different, we must interpolate.
    # Gridap.interpolate works best with a generic function for cross-mesh operations.
    # We construct a manual bilinear interpolator for robustness.
    
    interpolator = _create_bilinear_interpolator(ϕ, lsm_model)
    
    reffe_bg = ReferenceFE(lagrangian, Float64, 1)
    V_bg = FESpace(bg_model, reffe_bg)
    ϕ_bg_fe = interpolate(interpolator, V_bg)
    
    return GridapEmbedded.LevelSetCutters.DiscreteGeometry(ϕ_bg_fe, bg_model)
end

function _create_bilinear_interpolator(ϕ::Vector{Float64}, model::CartesianDiscreteModel)
    desc = get_cartesian_descriptor(model)
    origin = Tuple(desc.origin)
    spacing = Tuple(desc.sizes) # Cell sizes
    partition = Tuple(desc.partition) 
    nx, ny = partition .+ 1
    dims = (nx, ny)
    
    return x -> begin
        # InterpolationUtils expects a tuple/vector
        return bilinear_interpolate_scalar(ϕ, origin, spacing, dims, (x[1], x[2]))
    end
end
