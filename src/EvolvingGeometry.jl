# =============================================================================
# EvolvingGeometry.jl
# =============================================================================

using StaticArrays: SVector
using GridapEmbedded
import LevelSetMethods as LSM

"""
    GeometryCache{T_Trans, T_Ext}

Typed cache for memoized Transfer/Extension operators.
Uses granular versioning to track validity against `geom.step`.
"""
mutable struct GeometryCache{T_Trans, T_Ext}
    # Trackers (Initialize to -1 so Step 0 triggers a build)
    transfer_step::Int
    extension_step::Int

    # Artifacts (Allow them to be empty initially via Union{Nothing, T})
    transfer_op::Union{Nothing, T_Trans}
    extension_op::Union{Nothing, T_Ext}

    # Inner constructor
    function GeometryCache{T_Trans, T_Ext}() where {T_Trans, T_Ext}
        new{T_Trans, T_Ext}(-1, -1, nothing, nothing)
    end
end

"""
    EvolvingDiscreteGeometry{Eq, T_Trans, T_Ext}

Main container for evolving domain simulations.
"""
mutable struct EvolvingDiscreteGeometry{Eq, T_Trans, T_Ext}
    equation::Eq                           # LevelSetEquation from LevelSetMethods
    model::CartesianDiscreteModel          # Cartesian Gridap model
    velocity::Union{Nothing, AbstractVelocitySource}
    velocity_buffer::Vector{SVector{2,Float64}}  # Cached velocity for LSM update
    coords_tuples::Vector{NTuple{2,Float64}}     # Node coords for velocity sampling

    # Memoization Cache
    cache::GeometryCache{T_Trans, T_Ext}

    # Cached derived quantities
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
    EvolvingDiscreteGeometry(model, initial_sdf; velocity=nothing, reinit_freq=10, bc=NeumannBC())

Construct an evolving geometry on a Cartesian model.
"""
function EvolvingDiscreteGeometry(model::CartesianDiscreteModel, initial_sdf::Function;
                                   velocity::Union{Nothing, AbstractVelocitySource}=nothing,
                                   reinit_freq::Int=10,
                                   bc=LSM.NeumannBC())
    info = grid_info(model)

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
    coords_tuples = [NTuple{2,Float64}(Tuple(c)) for c in vec(collect(grid))]

    # Create velocity MeshField (initial velocity is zero)
    n_nodes = prod(partition)
    vel_buffer = fill(SVector{2,Float64}(0.0, 0.0), n_nodes)
    vel_array = reshape(vel_buffer, Tuple(partition))
    u_mesh = LSM.MeshField(vel_array, grid, nothing)

    update_fn = (u, ϕ_current, t) -> nothing
    terms = (LSM.AdvectionTerm(u_mesh, LSM.WENO5(), update_fn),)

    eq = LSM.LevelSetEquation(;
        terms = terms,
        levelset = ϕ,
        integrator = LSM.RK3(0.1),
        t = 0.0,
        bc = bc
    )

    # Initialize Cache with default types
    T_Trans = GridMeshTransfer{typeof(model)}
    T_Ext = ClosestPointExtension
    cache = GeometryCache{T_Trans, T_Ext}()

    return EvolvingDiscreteGeometry(
        eq, model, velocity, vel_buffer, coords_tuples,
        cache, nothing, true,
        0.0, 0,
        reinit_freq, 0
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
    geom.step += 1 # This invalidates the cache trackers!
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

Reinitialize the level set to a signed distance function.
Requires the ReinitializationExt extension from LevelSetMethods.jl
(auto-loaded when Interpolations and NearestNeighbors are present).
"""
function reinitialize!(geom::EvolvingDiscreteGeometry)
    # Check if extension is loaded
    if isempty(methods(LSM.reinitialize!))
        error("reinitialize! requires LevelSetMethods ReinitializationExt. " *
              "Ensure both `Interpolations` and `NearestNeighbors` are loaded.")
    end
    
    # Use invokelatest for world-age compatibility with extensions
    Base.invokelatest(LSM.reinitialize!, geom.equation)
    
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
    for (i, c) in enumerate(geom.coords_tuples)
        v = sample_velocity(vel, c, geom.t)
        geom.velocity_buffer[i] = SVector{2,Float64}(v[1], v[2])
    end
    return geom
end




"""
    get_transfer_op(geom::EvolvingDiscreteGeometry{Eq, T_Trans, T_Ext}) -> T_Trans

Retrieve the transfer operator from cache, or rebuild it if the cache is stale.
"""
function get_transfer_op(geom::EvolvingDiscreteGeometry{Eq, T_Trans, T_Ext}) where {Eq, T_Trans, T_Ext}
    # 1. Sync Check
    if geom.cache.transfer_step != geom.step
        # 2. Lazy Build
        new_op::T_Trans = setup_transfer(geom.model, geom.model)

        # 3. Update Cache
        geom.cache.transfer_op = new_op
        geom.cache.transfer_step = geom.step
    end
    # 4. Safe Return (Type Assertion)
    return geom.cache.transfer_op::T_Trans
end

"""
    prolong(geom::EvolvingDiscreteGeometry, u_mesh)

Transfer a field from the FE mesh to the fine level set grid.
Automatically manages cache synchronization.
"""
function TransferOperator.prolong(geom::EvolvingDiscreteGeometry, u_mesh)
    op = get_transfer_op(geom)
    return TransferOperator.prolong(op, u_mesh)
end

"""
    get_extension_op(geom::EvolvingDiscreteGeometry{Eq, T_Trans, T_Ext}) -> T_Ext

Retrieve the extension operator from cache, or rebuild it using Gridap gradients.
"""
function get_extension_op(geom::EvolvingDiscreteGeometry{Eq, T_Trans, T_Ext}) where {Eq, T_Trans, T_Ext}
    # 1. Sync Check
    if geom.cache.extension_step != geom.step
        # 2. Lazily Build Operator (Internal FD gradients)
        ϕ_vals = current_levelset(geom)
        info = grid_info(geom.model)
        new_op = ClosestPointExtension(info, ϕ_vals)

        # 3. Update Cache
        geom.cache.extension_op = new_op
        geom.cache.extension_step = geom.step
    end
    return geom.cache.extension_op::T_Ext
end

"""
    extend(geom::EvolvingDiscreteGeometry, u_grid)

Extend a grid-based field (e.g., velocity prolonged from mesh) to the void domain.
Uses closest-point mapping via cached extension operator.
"""
function extend(geom::EvolvingDiscreteGeometry, u_grid)
    op = get_extension_op(geom)
    return extend_field(op, u_grid)
end

export get_extension_op, extend

"""
    set_velocity!(geom::EvolvingDiscreteGeometry, u_ext::AbstractArray)

Write an extended velocity field into the geometry's internal buffer for advection.
Supports both Vector{SVector} (flat) and Matrix{SVector} (grid).
"""
function set_velocity!(geom::EvolvingDiscreteGeometry, u_ext::AbstractArray)
    # Explicitly convert to SVector to handle VectorValue or other types
    for (i, v) in enumerate(vec(u_ext))
        geom.velocity_buffer[i] = SVector{2,Float64}(v[1], v[2])
    end
    return geom
end

"""
    update_transfer_cache!(geom::EvolvingDiscreteGeometry, op)

Manually push a transfer operator (e.g., from CoarseAgFEM) into the cache.
Automatically synchronizes the transfer step.
"""
function update_transfer_cache!(geom::EvolvingDiscreteGeometry, op)
    geom.cache.transfer_op = op
    geom.cache.transfer_step = geom.step
    return geom
end

"""
    update_extension_cache!(geom::EvolvingDiscreteGeometry, op)

Manually push an extension operator into the cache.
Automatically synchronizes the extension step.
"""
function update_extension_cache!(geom::EvolvingDiscreteGeometry, op)
    geom.cache.extension_op = op
    geom.cache.extension_step = geom.step
    return geom
end

export set_velocity!, update_transfer_cache!, update_extension_cache!

# =============================================================================
# Internal Helpers
# =============================================================================

function _ensure_fresh!(geom::EvolvingDiscreteGeometry)
    if geom.dirty || geom.cached_cut === nothing
        # Build DiscreteGeometry inline (no separate function needed)
        ϕ = current_levelset(geom)
        reffe = ReferenceFE(lagrangian, Float64, 1)
        V = FESpace(geom.model, reffe)
        ϕ_fe = FEFunction(V, ϕ)
        geo = GridapEmbedded.LevelSetCutters.DiscreteGeometry(ϕ_fe, geom.model)
        geom.cached_cut = cut(geom.model, geo)
        geom.dirty = false
    end
    return nothing
end


