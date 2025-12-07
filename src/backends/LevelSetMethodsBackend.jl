# =============================================================================
# LevelSetMethods.jl Backend
# =============================================================================
# This backend uses LevelSetMethods.jl for level set evolution.
# It supports both static and time-dependent/FE-coupled velocities via
# a buffered MeshField approach with update_func callbacks.

# Note: LevelSetMethods.jl is loaded conditionally to avoid hard dependency
# Users must explicitly `using LevelSetMethods` before using this backend.

using StaticArrays

"""
    LevelSetMethodsEvolver <: AbstractLevelSetEvolver

A level set evolver backend using LevelSetMethods.jl.

# Fields
- `equation`: The LevelSetEquation from LevelSetMethods.jl
- `coords`: Node coordinates as Gridap VectorValues  
- `velocity_source`: The velocity source (AbstractVelocitySource or nothing)
- `velocity_buffer`: Mutable buffer for velocity values (used with update_func)

# Velocity Modes
- **Static**: Velocity sampled once at construction, never updated
- **Dynamic**: Velocity re-sampled before each RK stage via update_func callback

# Example
```julia
# Static velocity
evolver = LevelSetMethodsEvolver(;
    bg_model = model,
    initial_ls = ϕ0,
    velocity = StaticFunctionVelocity(x -> (1.0, 0.0))
)

# FE-coupled velocity
vel_source = FEVelocitySource(velocity_fh, model)
evolver = LevelSetMethodsEvolver(;
    bg_model = model,
    initial_ls = ϕ0,
    velocity = vel_source
)
```
"""
mutable struct LevelSetMethodsEvolver{Eq, C, V, B} <: AbstractLevelSetEvolver
    equation::Eq
    coords::C
    velocity_source::V      # Nothing or AbstractVelocitySource
    velocity_buffer::B      # Vector{SVector{2,Float64}} or nothing
end

# =============================================================================
# Constructor
# =============================================================================

"""
    LevelSetMethodsEvolver(;
        bg_model,
        initial_ls,
        velocity,
        integrator = :RK3,
        spatial_scheme = :WENO5,
        reinit_freq = nothing,
        bc = :Neumann
    )

Construct a LevelSetMethodsEvolver from a Gridap model and parameters.

# Arguments
- `bg_model`: CartesianDiscreteModel defining the computational domain
- `initial_ls`: Function `x -> Float64` defining the initial level set
- `velocity`: Either:
  - A function `x -> (vx, vy)` (static velocity, for backwards compatibility)
  - An `AbstractVelocitySource` (static, time-dependent, or FE-coupled)
- `integrator`: Time integrator `:ForwardEuler`, `:RK2`, or `:RK3` (default)
- `spatial_scheme`: Spatial scheme `:Upwind` or `:WENO5` (default)
- `bc`: Boundary conditions `:Periodic`, `:Neumann` (default), or `:Dirichlet`

Requires `using LevelSetMethods` before calling.
"""
function LevelSetMethodsEvolver(;
    bg_model::CartesianDiscreteModel,
    initial_ls::Function,
    velocity,
    integrator = :RK3,
    spatial_scheme = :WENO5,
    reinit_freq = nothing,
    bc = :Neumann
)
    # Check that LevelSetMethods is loaded
    if !isdefined(Main, :LevelSetMethods)
        error("LevelSetMethods.jl must be loaded. Run `using LevelSetMethods` first.")
    end
    LSM = Main.LevelSetMethods
    
    # Get grid parameters from Gridap model
    (origin, corner, partition) = cartesian_descriptor(bg_model)
    
    # LevelSetMethods uses node counts, Gridap uses cell counts
    node_partition = partition .+ 1
    
    # Create LevelSetMethods grid
    grid = LSM.CartesianGrid(origin, corner, node_partition)
    
    # Create initial level set
    ϕ = LSM.LevelSet(initial_ls, grid)
    
    # Get Gridap-compatible coordinates
    coords = get_node_coords_as_vector(bg_model)
    
    # Convert grid coords to tuple format for velocity sampling
    grid_coords_tuples = [Tuple(c) for c in coords]
    
    # Wrap velocity in VelocitySource if needed (backwards compatibility)
    vel_source = if velocity isa AbstractVelocitySource
        velocity
    elseif velocity isa Function
        # Infer if time-dependent by checking method signature
        if _has_time_arg(velocity)
            TimeDependentVelocity(velocity)
        else
            StaticFunctionVelocity(velocity)
        end
    else
        error("velocity must be a Function or AbstractVelocitySource")
    end
    
    # Select integrator
    int = if integrator == :RK2
        LSM.RK2()
    elseif integrator == :RK3
        LSM.RK3()
    elseif integrator == :ForwardEuler
        LSM.ForwardEuler()
    else
        LSM.RK3()
    end
    
    # Select boundary conditions
    boundary = if bc == :Periodic
        LSM.PeriodicBC()
    elseif bc == :Neumann
        LSM.NeumannBC()
    elseif bc == :Dirichlet
        LSM.DirichletBC()
    else
        LSM.NeumannBC()
    end
    
    # Select spatial scheme
    scheme = spatial_scheme == :Upwind ? LSM.Upwind() : LSM.WENO5()
    
    # Determine if we need dynamic velocity updates
    needs_update = is_time_dependent(vel_source)
    
    # Create velocity representation for LevelSetMethods
    # For static velocity: use MeshField (sampled once, most efficient)
    # For dynamic velocity: use MeshField + update_func (buffered, ~10x faster than function)
    grid_coords_tuples = [Tuple(c) for c in coords]
    
    # Sample initial velocity
    initial_vel = sample_velocity(vel_source, grid_coords_tuples, 0.0)
    
    # Create velocity buffer (will be mutated by update_func for dynamic velocities)
    vel_buffer = copy(initial_vel)
    
    # Create MeshField from buffer
    # The array is reshaped to match grid dimensions; buffer shares underlying data
    vel_array = reshape(vel_buffer, Tuple(node_partition))
    u_mesh = LSM.MeshField(vel_array, grid, nothing)
    
    # Create update function for dynamic velocities
    # This is called before each RK stage, allowing us to sample FE velocity once per stage
    update_fn = if needs_update
        # Closure that captures vel_source, grid_coords_tuples, and vel_buffer
        (u, ϕ_current, t) -> begin
            new_vel = sample_velocity(vel_source, grid_coords_tuples, t)
            # Update buffer in-place (MeshField.vals points to vel_array which aliases vel_buffer)
            for i in eachindex(new_vel)
                vel_buffer[i] = new_vel[i]
            end
            nothing
        end
    else
        (u, ϕ_current, t) -> nothing  # No-op for static velocity
    end
    
    # Create AdvectionTerm with update_func (requires LevelSetMethods from main branch)
    terms = (LSM.AdvectionTerm(u_mesh, scheme, update_fn),)
    
    # Create equation
    eq = LSM.LevelSetEquation(;
        terms = terms,
        levelset = ϕ,
        integrator = int,
        bc = boundary,
        t = 0.0
    )
    
    return LevelSetMethodsEvolver(eq, coords, vel_source, vel_buffer)
end

"""Check if a function accepts a time argument (x, t) vs just (x)."""
function _has_time_arg(f::Function)
    # Try calling with 2 args - if it works, it's time-dependent
    try
        # Use a test point
        f((0.0, 0.0), 0.0)
        return true
    catch e
        if e isa MethodError
            return false
        end
        rethrow()
    end
end

# =============================================================================
# Interface Implementation
# =============================================================================

function evolve!(e::LevelSetMethodsEvolver, Δt::Real)
    LSM = Main.LevelSetMethods
    tf = LSM.current_time(e.equation) + Δt
    LSM.integrate!(e.equation, tf)
    return e
end

function current_values(e::LevelSetMethodsEvolver)
    LSM = Main.LevelSetMethods
    state = LSM.current_state(e.equation)
    # LevelSetMethods stores values as a matrix, iterate in grid order
    vals = LSM.values(state)
    # Flatten in column-major order (Julia default) which matches grid iteration
    return vec(vals)
end

function current_time(e::LevelSetMethodsEvolver)
    LSM = Main.LevelSetMethods
    return LSM.current_time(e.equation)
end

function reinitialize!(e::LevelSetMethodsEvolver)
    LSM = Main.LevelSetMethods
    LSM.reinitialize!(LSM.current_state(e.equation))
    return e
end

function grid_coords(e::LevelSetMethodsEvolver)
    # Return coordinates in the same order as level set values
    LSM = Main.LevelSetMethods
    state = LSM.current_state(e.equation)
    grid = LSM.mesh(state)
    # Extract grid points as VectorValues
    # LevelSetMethods uses column-major ordering
    coords = Vector{VectorValue{2,Float64}}()
    for pt in grid
        push!(coords, VectorValue(pt...))
    end
    return coords
end

# =============================================================================
# Velocity Update Support
# =============================================================================

supports_velocity_update(e::LevelSetMethodsEvolver) = !isnothing(e.velocity_source)

"""
    update_velocity!(evolver::LevelSetMethodsEvolver, new_source::AbstractVelocitySource, t)

Update the evolver's velocity source. The new source will be sampled at the
next RK stage via the update_func callback.
"""
function update_velocity!(e::LevelSetMethodsEvolver, new_source::AbstractVelocitySource, t)
    e.velocity_source = new_source
    return e
end

"""
    update_velocity!(evolver::LevelSetMethodsEvolver, fh, t)

Convenience method: if an FEVelocitySource is already set, update its wrapped function.
Otherwise, create a new FEVelocitySource.
"""
function update_velocity!(e::LevelSetMethodsEvolver, fh, t)
    if e.velocity_source isa FEVelocitySource
        update_velocity!(e.velocity_source, fh)
    else
        # Cannot update - source type mismatch
        @warn "Cannot update velocity: evolver velocity source is not FEVelocitySource"
    end
    return e
end
