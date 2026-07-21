module Geometric

using StaticArrays: SVector
using GridapEmbedded
using Gridap: CartesianDiscreteModel
using Gridap.Geometry: get_node_coordinates, get_grid_topology, get_faces, num_nodes
using GridapEmbedded: cut
using Gridap.TensorValues: VectorValue

include("GridInfo.jl")
include("CartesianField.jl")
include("Stencils.jl")
include("GeometryDesign.jl")
include("TopologyFilter.jl")

using .CartesianField
using .Stencils

export CartesianMeshField, meshsize, get_interpolator
export D⁺, D⁻, D⁰, D²
export weno5⁻, weno5⁺

export GeometryCache, EvolvingDiscreteGeometry
export set_levelset!, current_levelset, current_cut, get_active_indices
export CartesianGridInfo, grid_info
export AbstractGeometry, Circle, Rectangle, Translate, signed_distance

export filter_small_phase_islands!

mutable struct WENO5Cache
    rhs::Union{Nothing, Vector{Float64}}
    stage::Union{Nothing, Vector{Float64}}
    phi0::Union{Nothing, Vector{Float64}}

    WENO5Cache() = new(nothing, nothing, nothing)
end
export WENO5Cache


mutable struct GeometryCache
    # Cross-module operator caches (typed as Any to avoid circular module deps).
    # transfer_op is a GridMeshTransfer; extension_op is a ClosestPointExtension.
    transfer_op::Any
    extension_op::Any

    # Current and previous EmbeddedDiscretization (computed lazily on first access)
    cut::Union{Nothing,GridapEmbedded.Interfaces.EmbeddedDiscretization}
    prev_cut::Union{Nothing,GridapEmbedded.Interfaces.EmbeddedDiscretization}

    # Cached active node indices — recomputed when `cut` is invalidated.
    # Avoids repeated topology walks (get_faces + unique!) across multiple callers.
    active_nodes::Union{Nothing,Vector{Int}}

    # Persistent scratch buffer for the WENO5 RK3 stages.
    # Allocated lazily on first advance! call; reused every subsequent step.
    # Eliminates one Vector{Float64} allocation per advance! call.
    weno_scratch::Union{Nothing,Vector{Float64}}

    # Per-geometry WENO5 cache for RK3 stage buffers.
    # Allocated once at construction, eliminates the need for a global singleton.
    weno_cache::WENO5Cache

    # Interface position/normal samples for curvature. Typed as Any to avoid a
    # circular dependency on the Curvature submodule, which is included below.
    interface_samples::Any

    function GeometryCache()
        new(nothing, nothing, nothing, nothing, nothing, nothing, WENO5Cache(), nothing)
    end
end

mutable struct EvolvingDiscreteGeometry
    levelset::Vector{Float64}
    grid::CartesianDiscreteModel
    cache::GeometryCache
end

"""
    EvolvingDiscreteGeometry(levelset::Vector{Float64}, grid::CartesianDiscreteModel)

Create a geometry object from a levelset vector and a grid.
The cache is initialized automatically.
"""
function EvolvingDiscreteGeometry(levelset::Vector{Float64}, grid::CartesianDiscreteModel)
    return EvolvingDiscreteGeometry(levelset, grid, GeometryCache())
end

"""
    EvolvingDiscreteGeometry(grid::CartesianDiscreteModel)

Initialize an empty geometry (all zeros) for a given grid.
Useful for cases where set_levelset! will be called immediately after.
"""
function EvolvingDiscreteGeometry(grid::CartesianDiscreteModel)
    n_nodes = num_nodes(grid)
    phi = zeros(n_nodes)
    return EvolvingDiscreteGeometry(phi, grid)
end

"""
    current_levelset(geom::EvolvingDiscreteGeometry) -> Vector{Float64}

Return the current level set field.
"""
current_levelset(geom::EvolvingDiscreteGeometry) = geom.levelset

"""
    current_cut(geom::EvolvingDiscreteGeometry)

Return the current GridapEmbedded cut geometry.
"""
current_cut(geom::EvolvingDiscreteGeometry) = geom.cache.cut

"""
    set_levelset!(geom::EvolvingDiscreteGeometry, ϕ::Vector{Float64})

Update the level set function and invalidate the geometry cache.
Previous geometry cut is moved to `prev_cut` for history access.
"""
function set_levelset!(geom::EvolvingDiscreteGeometry, ϕ::Vector{Float64})
    if length(ϕ) != length(geom.levelset)
        error("set_levelset!: dimension mismatch — expected $(length(geom.levelset)) values, got $(length(ϕ)).")
    end
    copyto!(geom.levelset, ϕ)

    # Shift current cut to history
    geom.cache.prev_cut = geom.cache.cut

    invalidate!(geom.cache)
    return geom
end

"""
    ensure_cut!(geom::EvolvingDiscreteGeometry)

Return the current cut, computing and caching it if the cache was invalidated.
"""
function ensure_cut!(geom::EvolvingDiscreteGeometry)
    if isnothing(geom.cache.cut)
        flat_coords = vec(collect(get_node_coordinates(geom.grid)))
        geo = GridapEmbedded.LevelSetCutters.DiscreteGeometry(geom.levelset, flat_coords)
        geom.cache.cut = cut(geom.grid, geo)
    end
    return geom.cache.cut
end

"""
    invalidate!(cache::GeometryCache)

Drop every cached quantity derived from the level set.

Kept in one place because the reset happens from several call sites; adding a cached
field and updating only some of them would leave stale geometry in circulation.
"""
function invalidate!(cache::GeometryCache)
    cache.cut = nothing
    cache.active_nodes = nothing
    cache.transfer_op = nothing
    cache.extension_op = nothing
    cache.interface_samples = nothing
    return cache
end

"""
    get_active_indices(geom::EvolvingDiscreteGeometry, state::Symbol=:current) -> Vector{Int}

Return the linear indices of the active **nodes** (nodes belonging to IN or CUT cells) for the specified state.
State can be `:current` or `:prev`.

The result for `:current` is cached in `geom.cache.active_nodes` and reused until the
next `set_levelset!` or `reinitialize!` call, avoiding repeated topology walks.
"""
function get_active_indices(geom::EvolvingDiscreteGeometry, state::Symbol=:current)
    cut_geo = nothing
    if state == :current
        if isnothing(geom.cache.cut)
            # We need to wrap the levelset vector into a geometry object that 'cut' accepts.
            flat_coords = vec(collect(get_node_coordinates(geom.grid)))
            geo = GridapEmbedded.LevelSetCutters.DiscreteGeometry(geom.levelset, flat_coords)
            geom.cache.cut = cut(geom.grid, geo)
        end

        # Return cached active nodes if available (avoids topology walk on repeated calls)
        if !isnothing(geom.cache.active_nodes)
            return geom.cache.active_nodes
        end

        cut_geo = geom.cache.cut
    elseif state == :prev
        cut_geo = geom.cache.prev_cut
        if isnothing(cut_geo)
            return Int[]
        end
    else
        error("Unknown state: $state. Use :current or :prev.")
    end

    # 1. Access Cell Classification
    const_IN = -1
    const_CUT = 0

    raw_status = cut_geo.ls_to_bgcell_to_inoutcut
    status = if eltype(raw_status) <: AbstractVector
        raw_status[1]
    else
        raw_status
    end

    active_cell_indices = findall(x -> x == const_IN || x == const_CUT, status)

    if isempty(active_cell_indices)
        if state == :current
            geom.cache.active_nodes = Int[]
        end
        return Int[]
    end

    # 2. Map Cells -> Nodes (Topology)
    grid = geom.grid
    topo = get_grid_topology(grid)
    cell_to_nodes = get_faces(topo, 2, 0)

    active_nodes = Int[]
    sizehint!(active_nodes, length(active_cell_indices) * 4)
    for cell_idx in active_cell_indices
        nodes = cell_to_nodes[cell_idx]
        append!(active_nodes, nodes)
    end
    sort!(unique!(active_nodes))

    # Cache the result so subsequent calls skip the topology walk
    if state == :current
        geom.cache.active_nodes = active_nodes
    end

    return active_nodes
end

# Curvature first: it owns the explicit reference curve (position, normal, curvature per
# subfacet) that reinitialisation also seeds from, so both read one cached object.
include("Curvature.jl")
using .Curvature

include("Reinitialization.jl")
using .Reinitialization
export reinitialize!
export InterfaceSamples, interface_samples, get_curvature, curvature_at
export reference_curve

"""
    interface_curvature(geom::EvolvingDiscreteGeometry; radius=nothing) -> Vector{Float64}

Curvature at each subfacet of the current embedded interface, in subfacet order.

Matches `num_cells(EmbeddedBoundary(cut))`, so it can be passed straight to
`CellField(κ, Γ)` and integrated against a `Measure` on Γ. `radius` defaults to four
cells; see `get_curvature` for what it trades off. Samples are cached and dropped
whenever the level set changes.
"""
function interface_curvature(geom::EvolvingDiscreteGeometry; radius=nothing)
    curve = reference_curve(geom)
    # The reference curve already carries curvature, fitted at 8h — 8h and not 4h because
    # the fit is degree 4, and order and window must rise together (degree 4 at 4h is
    # starved: κ error ~5-7% on a circle; at 8h it is ~1%). Only an explicit override refits.
    isnothing(radius) && return curve.curvatures
    return get_curvature(curve; radius=radius)
end
export interface_curvature, ensure_cut!, invalidate!

"""
    tangential_smooth!(geom::EvolvingDiscreteGeometry; strength=0.05, band=3) -> geom

Diffuse the level set *along* the interface, damping sub-grid wiggles that surface tension
cannot see.

Applies one explicit step of the tangential heat equation φ ← φ + ν ∂ττφ, where ∂ττ is the
second derivative along the tangent t = ∇φ⊥ / |∇φ|, over the band `|φ| < band·h`. It is
**tangential by construction**: for a clean signed-distance circle φ is constant along the
tangent, so ∂ττφ = 0 and the interface does not move — only variation *along* Γ is damped.

`strength` is dimensionless; the diffusivity is ν = strength·h² (h = min grid spacing), so
the operator is O(h²) and **vanishes under refinement** — a consistent regulariser, not a
shape change. Stability needs strength ≲ 0.1.

**Why this exists.** The curvature is fitted over an ~8h window, so it low-passes the
interface and, near its cutoff, its gain sign-flips — turning the physical decay s_m ∝ −m³
of short-wavelength ripples into growth. This restores the dissipation the window discards.
The principled alternative is a semi-implicit (Laplace–Beltrami) treatment of surface
tension, which this package does not yet provide accurately on cut cells; until it does,
this is the working regulariser for explicit surface-tension flows.

**Caveat — it is not volume-conserving.** Removing a wiggle of amplitude a from
r = R + a cos(mθ) sheds enclosed area πa²/2 (the incompressible physics would grow R to
compensate; a geometric filter does not). If area matters, follow this with a global normal
shift φ += (A − A₀)/P to restore it.
"""
function tangential_smooth!(geom::EvolvingDiscreteGeometry; strength::Real=0.05, band::Real=3)
    phi = geom.levelset
    info = grid_info(geom.grid)
    nx, ny = info.dims
    hx, hy = info.spacing
    hmin = min(hx, hy)
    ν = strength * hmin^2
    bw = band * hmin

    out = copy(phi)
    @inbounds for j in 2:ny-1, i in 2:nx-1
        k = i + (j - 1) * nx
        abs(phi[k]) < bw || continue
        gx = (phi[k+1] - phi[k-1]) / (2hx)
        gy = (phi[k+nx] - phi[k-nx]) / (2hy)
        g = hypot(gx, gy)
        g > 1e-12 || continue
        tx, ty = -gy / g, gx / g                     # unit tangent = ∇φ rotated 90°
        fxx = (phi[k+1]  - 2phi[k] + phi[k-1])  / hx^2
        fyy = (phi[k+nx] - 2phi[k] + phi[k-nx]) / hy^2
        fxy = (phi[k+nx+1] - phi[k+nx-1] - phi[k-nx+1] + phi[k-nx-1]) / (4hx * hy)
        out[k] = phi[k] + ν * (tx*tx*fxx + 2tx*ty*fxy + ty*ty*fyy)   # φ + ν ∂ττφ
    end
    copyto!(phi, out)

    invalidate!(geom.cache)
    return geom
end
export tangential_smooth!

end # module
