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
    bg_model::DiscreteModel                # Background Gridap model
    lsm_model::CartesianDiscreteModel      # Cartesian Gridap model for Level Set Evolution
    velocity::Union{Nothing, AbstractVelocitySource}
    velocity_buffer::Vector{SVector{2,Float64}}  # Cached velocity for LSM update
    coords_tuples::Vector{NTuple{2,Float64}}     # Node coords for velocity sampling

    # Memoization Cache
    cache::GeometryCache{T_Trans, T_Ext}

    # Cached derived quantities # Legacy ?
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
                                   reinit_freq::Int=10,
                                   bc=LSM.NeumannBC())
    return EvolvingDiscreteGeometry(model, model, initial_sdf; velocity=velocity, reinit_freq=reinit_freq, bc=bc)
end

"""
    EvolvingDiscreteGeometry(bg_model, lsm_model, initial_sdf; velocity=nothing, reinit_freq=10)

Construct an evolving geometry with decoupled models:
- `bg_model`: The background model for CutFEM (can be unstructured, Quadtree, etc.)
- `lsm_model`: The Cartesian model for Level Set evolution (MUST be Cartesian)
"""
function EvolvingDiscreteGeometry(bg_model::DiscreteModel, lsm_model::CartesianDiscreteModel, initial_sdf::Function;
                                   velocity::Union{Nothing, AbstractVelocitySource}=nothing,
                                   reinit_freq::Int=10,
                                   bc=LSM.NeumannBC())
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

    # Create velocity MeshField
    # (Initial velocity is zero; update explicitly via Bridge API)
    n_nodes = prod(partition)
    vel_buffer = fill(SVector{2,Float64}(0.0, 0.0), n_nodes)
    vel_array = reshape(vel_buffer, Tuple(partition))
    u_mesh = LSM.MeshField(vel_array, grid, nothing)

    # We no longer provide a complex update_fn inside the constructor.
    update_fn = (u, ϕ_current, t) -> nothing

    # Create advection term (with update_func)
    terms = (LSM.AdvectionTerm(u_mesh, LSM.WENO5(), update_fn),)

    eq = LSM.LevelSetEquation(;
        terms = terms,
        levelset = ϕ,
        integrator = LSM.RK3(0.1),
        t = 0.0,
        bc = bc
    )

    # Initialize Cache with default types
    T_Trans = GridMeshTransfer{typeof(bg_model)}
    T_Ext = ClosestPointExtension
    cache = GeometryCache{T_Trans, T_Ext}()

    return EvolvingDiscreteGeometry(
        eq, bg_model, lsm_model, velocity, vel_buffer, coords_tuples,
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
    # Reinitialize the level set to a signed distance function
    # Newton-based reinitialization provided by ReinitializationExt
    try
        phi = LSM.current_state(geom.equation)
        # Use invokelatest to ensure world-age compatibility with extensions
        Base.invokelatest(LSM.reinitialize!, phi)
    catch e
        @warn "Reinitialization dispatch failed. Check extension loading. Error: $e"
    end
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
    grid_info(geom::EvolvingDiscreteGeometry) -> CartesianGridInfo
"""
grid_info(geom::EvolvingDiscreteGeometry) = grid_info(geom.lsm_model)

"""
    get_transfer_op(geom::EvolvingDiscreteGeometry{Eq, T_Trans, T_Ext}) -> T_Trans

Retrieve the transfer operator from cache, or rebuild it if the cache is stale.
"""
function get_transfer_op(geom::EvolvingDiscreteGeometry{Eq, T_Trans, T_Ext}) where {Eq, T_Trans, T_Ext}
    # 1. Sync Check
    if geom.cache.transfer_step != geom.step
        # 2. Lazy Build
        new_op::T_Trans = setup_transfer(geom.lsm_model, geom.bg_model)

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
        # 2. Lazy Build via Gridap Gradients
        ϕ_vals = current_levelset(geom)
        info = grid_info(geom)

        # Create FE Function on LSM model
        reffe = ReferenceFE(lagrangian, Float64, 1)
        V = FESpace(geom.lsm_model, reffe)
        ϕ_fe = FEFunction(V, ϕ_vals)

        # Compute exact gradient
        # Compute exact gradient and interpolate to continuous space for nodal evaluation
        grad_ϕ_field = gradient(ϕ_fe)
        reffe_grad = ReferenceFE(lagrangian, VectorValue{2,Float64}, 1)
        V_grad = FESpace(geom.lsm_model, reffe_grad)
        grad_ϕ_fe = interpolate(grad_ϕ_field, V_grad)

        # Sample gradient at nodes
        # (Gridap nodes match LSM nodes for CartesianDiscreteModel)
        node_coords = get_node_coordinates(geom.lsm_model)
        grad_vals = [grad_ϕ_fe(x) for x in node_coords]

        new_op = ClosestPointExtension(info, ϕ_vals, vec(grad_vals))

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
