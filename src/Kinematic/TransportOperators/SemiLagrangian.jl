module SemiLagrangian

using Gridap
using Gridap.Geometry: get_node_coordinates, num_nodes
using Gridap.TensorValues
using StaticArrays

using ...Geometric: CartesianMeshField, CartesianGridInfo, EvolvingDiscreteGeometry,
    get_active_indices, get_interpolator, quadratic_interpolation_weights, grid_info

using ..Kinematic: AbstractVelocitySource, get_velocity, is_time_dependent

export TransportMap, advect!

#This module follows the method described in the paper

#An unconditionally stable fully conservative semi-Lagrangian method
#Lentine, Grétarsson, Fedkiw 2011
# https://doi.org/10.1016/j.jcp.2010.12.036
#
# The important details is that we integrate over the whole neigborhood of the backward ray destination using a supersampling.
# This avoids making too much mass dissapear in the first phase due to diverging characteristics and being obligated to compensate in the phase 𝟚
# Which could lead to unphysical results.

"""
    TransportMap(geom, velocity, dt; cache_weights::Bool=false)

A discretized representation of the flow between two time steps.
Contains all geometric and kinematic information required for conservative advection.
This object is field-independent and should be reused for all fields advecting
with the same velocity (e.g., components of a strain tensor).

The canonical velocity input is a frozen nodal field over the step. It is interpolated
bilinearly for off-grid characteristic tracing, so map construction does not require
physical time. A time-independent `AbstractVelocitySource` is also accepted directly.

The geometry must carry both time levels. Materialize its current cut before updating
the level set so that the update can preserve it as `prev_cut`.

Pass `cache_weights=true` to additionally store the conservative pull stencils,
which makes repeated `advect!` calls on the same map cheaper at the cost of extra
map memory. The default (`false`) leaves the map unextended and recomputes weights
per call, exactly as before.

"""
struct TransportMap
    # All fields below are construction snapshots: rays, masks, offsets and the
    # packed stencils are read during advect! but are not designed to be mutated.
    # (source_mask is the exception: advect! detects a changed mask and falls back.)
    #
    # --- Backward (Pull Phase) ---
    # target_active_idx[k] -> bundle of 4 departure_points
    active_indices::Vector{Int}
    backward_rays::Vector{SVector{4,Point{2,Float64}}}

    # --- Conservation ---
    # source_idx -> total_weight_pulled (Used to scale weights for conservation)
    demand_map::Vector{Float64}

    # --- Optional packed pull-weight cache ---
    # compute_conservative_weights(x_dep, meta, source_mask) is step-invariant:
    # it depends only on the cached backward rays, the grid and the source
    # support, never on the field or dt.  With `cache_weights=true` its nonzero
    # outputs are stored flat, in the original (ray 1..4, stencil 1..16) order:
    #
    #   wq = 0.25 * w_m   (stencil weight w_m with the supersample share folded
    #                      in; 0.25 is a power of two, so this is exact)
    #   weight_offsets[k] <= p < weight_offsets[k+1]  selects node k's entries.
    #
    # Dropping wq == 0 (they add nothing) and keeping the order makes the packed
    # sum reproduce the uncached loop bit for bit.  A 4x4 quadratic stencil has
    # at most 8 positive entries per ray, i.e. at most 32 pairs per node, instead
    # of the 64+64 of a dense cache.
    weight_offsets::Vector{Int}      # empty unless cache_valid
    weight_indices::Vector{Int32}    # linear source-node index per entry
    weight_values::Vector{Float64}   # wq = 0.25 * w_m per entry
    # source_mask snapshot; advect! ignores the cache when the live mask differs,
    # so editing the support never silently reuses stale weights.
    cached_source_mask::BitVector
    # False when cache_weights was not requested (default) or n_nodes exceeds
    # typemax(Int32); advect! then recomputes, so indices are never truncated.
    cache_valid::Bool

    # --- Forward (Push Phase / Leakage Correction) ---
    # source_idx -> arrival_point (Only for nodes not fully resolved by Pull phase)
    leakage_indices::Vector{Int}
    leakage_rays::Vector{Point{2,Float64}}

    # --- Metadata ---
    grid_meta::CartesianGridInfo
    # Source and target supports are the discrete geometry at t⁻ and t⁺.
    # The level set itself is deliberately not retained in the map.
    source_mask::BitVector
    target_mask::BitVector
    is_identity::Bool
end

function TransportMap(
    geom::EvolvingDiscreteGeometry,
    velocity::AbstractVector{<:VectorValue{2}},
    dt::Real;
    cache_weights::Bool=false,
)
    Base.require_one_based_indexing(velocity)
    meta = grid_info(geom.grid)
    expected = prod(meta.dims)
    length(velocity) == expected || throw(DimensionMismatch(
        "TransportMap velocity has $(length(velocity)) values; grid requires $expected."))
    all(v -> isfinite(v[1]) && isfinite(v[2]), velocity) ||
        throw(ArgumentError("TransportMap velocity contains non-finite values."))

    sampled = velocity isa Vector{VectorValue{2,Float64}} ? velocity :
        [VectorValue(Float64(v[1]), Float64(v[2])) for v in velocity]
    interpolator = get_interpolator(CartesianMeshField(sampled, meta))
    return _build_transport_map(geom, x -> interpolator(x[1], x[2]), dt; cache_weights)
end

function TransportMap(
    geom::EvolvingDiscreteGeometry,
    velocity::CartesianMeshField{T},
    dt::Real;
    cache_weights::Bool=false,
) where {T<:VectorValue{2}}
    meta = grid_info(geom.grid)
    velocity.grid == meta || throw(ArgumentError(
        "TransportMap velocity belongs to a different Cartesian grid."))
    return TransportMap(geom, velocity.data, dt; cache_weights)
end

function TransportMap(
    geom::EvolvingDiscreteGeometry,
    velocity::AbstractVelocitySource,
    dt::Real;
    cache_weights::Bool=false,
)
    is_time_dependent(velocity) && throw(ArgumentError(
        "TransportMap requires a frozen velocity; call sample_velocity at the desired time."))
    return _build_transport_map(geom, x -> get_velocity(velocity, x, 0.0), dt; cache_weights)
end

# Backwards-compatible raw positional constructor: returns an uncached map
# (empty packed arrays, cache_valid=false) so existing 9-field callers work.
function TransportMap(
    active_indices::Vector{Int},
    backward_rays::Vector{SVector{4,Point{2,Float64}}},
    demand_map::Vector{Float64},
    leakage_indices::Vector{Int},
    leakage_rays::Vector{Point{2,Float64}},
    grid_meta::CartesianGridInfo,
    source_mask::BitVector,
    target_mask::BitVector,
    is_identity::Bool,
)
    return TransportMap(active_indices, backward_rays, demand_map,
        Int[], Int32[], Float64[], BitVector(), false,
        leakage_indices, leakage_rays, grid_meta, source_mask, target_mask, is_identity)
end

function _build_transport_map(geom::EvolvingDiscreteGeometry, velocity_at, dt::Real;
    cache_weights::Bool=false)
    isfinite(dt) || throw(ArgumentError("TransportMap time step must be finite."))
    dt >= 0 || throw(ArgumentError("TransportMap time step must be non-negative."))

    grid = geom.grid
    meta = grid_info(grid)
    coords = get_node_coordinates(grid)
    n_nodes = num_nodes(grid)

    isnothing(geom.cache.prev_cut) && throw(ArgumentError(
        "TransportMap requires the previous geometry cut; call ensure_cut!(geom) " *
        "before updating its level set."))

    # Get nodes belonging to the previous and current geometries.
    active_current = get_active_indices(geom, :current)
    active_previous = get_active_indices(geom, :prev)

    source_mask = falses(n_nodes)
    target_mask = falses(n_nodes)
    source_mask[active_previous] .= true
    target_mask[active_current] .= true

    # 1. Backward Flow (Where active nodes come from)
    dx, dy = meta.spacing
    offsets = SVector(
        VectorValue(-0.25 * dx, -0.25 * dy), VectorValue(0.25 * dx, -0.25 * dy),
        VectorValue(-0.25 * dx, 0.25 * dy), VectorValue(0.25 * dx, 0.25 * dy)
    )

    stationary = true
    backward_rays = Vector{SVector{4,Point{2,Float64}}}(undef, length(active_current))
    for (k, i) in enumerate(active_current)
        rays_buffer = MVector{4,Point{2,Float64}}(undef)
        for j in 1:4
            x_departure = coords[i] + offsets[j]
            rays_buffer[j] = trace_ray(x_departure, velocity_at, -dt)
            stationary &= rays_buffer[j] == x_departure
        end
        backward_rays[k] = SVector(rays_buffer)
    end

    if stationary
        for i in eachindex(source_mask)
            source_mask[i] && (stationary &= trace_ray(coords[i], velocity_at, dt) == coords[i])
        end
    end
    stationary &= source_mask == target_mask

    # 2. Conservation Demand, plus the optional packed pull weights.  The weights
    # are only materialized when requested (default off) and when the linear node
    # index fits in Int32, so the default map allocates no large cache arrays.
    cache_valid = cache_weights && n_nodes <= typemax(Int32)
    n_active = length(backward_rays)
    demand = zeros(Float64, n_nodes)
    if cache_valid
        weight_offsets = Vector{Int}(undef, n_active + 1)
        weight_offsets[1] = 1
        weight_indices = Int32[]
        weight_values = Float64[]
        sizehint!(weight_indices, 32 * n_active)
        sizehint!(weight_values, 32 * n_active)
    else
        weight_offsets = Int[]
        weight_indices = Int32[]
        weight_values = Float64[]
    end
    for k in 1:n_active
        for x_dep in backward_rays[k]
            indices, weights = compute_conservative_weights(x_dep, meta, source_mask)
            for m in 1:16
                s_idx = indices[m]
                wq = 0.25 * weights[m]
                source_mask[s_idx] && (demand[s_idx] += wq)
                if cache_valid && wq > 0.0
                    push!(weight_indices, Int32(s_idx))
                    push!(weight_values, wq)
                end
            end
        end
        cache_valid && (weight_offsets[k + 1] = length(weight_indices) + 1)
    end

    # 3. Leakage Map (Forward rays for mass not 'pulled' by Pass 1)
    leak_idx = Int[]
    leak_rays = Point{2,Float64}[]
    for i in 1:n_nodes
        source_mask[i] || continue
        # If demand < 1.0, some mass at this source node might be left behind
        if demand[i] < 1.0
            push!(leak_idx, i)
            push!(leak_rays, trace_ray(coords[i], velocity_at, dt))
        end
    end

    return TransportMap(active_current, backward_rays, demand,
        weight_offsets, weight_indices, weight_values,
        cache_valid ? copy(source_mask) : BitVector(), cache_valid,
        leak_idx, leak_rays, meta, source_mask, target_mask, stationary)
end


"""
    advect!(target_data, source_data, map; type=:conservative)

Apply conservative CCISL redistribution defined by `map`. For a valid map, this
preserves the sum over its source support; use the velocity-based `advect!` overload
for intensive CIP.
"""
function advect!(
    target_data::Vector{Float64},
    source_data::Vector{Float64},
    map::TransportMap;
    type::Symbol=:conservative,
)
    type === :conservative || throw(ArgumentError(
        "TransportMap advect! supports type=:conservative; " *
        "pass a frozen velocity and dt for type=:intensive."))
    expected = length(map.source_mask)
    length(target_data) == expected || throw(DimensionMismatch(
        "advect!: target has $(length(target_data)) values; map requires $expected."))
    length(source_data) == expected || throw(DimensionMismatch(
        "advect!: source has $(length(source_data)) values; map requires $expected."))
    Base.mightalias(target_data, source_data) && throw(ArgumentError(
        "Conservative advect! source and target must not alias."))
    (Base.mightalias(target_data, map.demand_map) ||
     Base.mightalias(source_data, map.demand_map) ||
     Base.mightalias(target_data, map.weight_values) ||
     Base.mightalias(source_data, map.weight_values)) && throw(ArgumentError(
        "Conservative advect! fields must not alias TransportMap storage."))
    for i in eachindex(source_data)
        map.source_mask[i] && !isfinite(source_data[i]) && throw(ArgumentError(
            "Conservative advect! source contains non-finite values on its support."))
    end

    if map.is_identity
        fill!(target_data, 0.0)
        target_data[map.active_indices] .= source_data[map.active_indices]
        return target_data
    end

    fill!(target_data, 0.0)

    # --- Pass 1: Backward Pull (with 2x2 Supersampling) ---
    # The packed stencils are used only when the map carries a valid cache and
    # the live support still matches its construction snapshot; otherwise the
    # stencils are recomputed exactly as before.
    use_cache = map.cache_valid && map.source_mask == map.cached_source_mask
    if use_cache
        Base.Threads.@threads for k in eachindex(map.active_indices)
            target_idx = map.active_indices[k]
            val_accum = 0.0
            @inbounds for p in map.weight_offsets[k]:(map.weight_offsets[k + 1] - 1)
                s_idx = Int(map.weight_indices[p])
                wq = map.weight_values[p]
                req = map.demand_map[s_idx]
                # Scale by 1/demand on the source node if over-requested (mass
                # conservation).  wq already carries the 0.25 supersample share,
                # which is exact in binary.
                scale = req > 1.0 ? (1.0 / req) : 1.0
                val_accum += wq * scale * source_data[s_idx]
            end
            target_data[target_idx] = val_accum
        end
    else
        Base.Threads.@threads for k in eachindex(map.active_indices)
            target_idx = map.active_indices[k]
            rays = map.backward_rays[k]
            val_accum = 0.0

            for x_dep in rays
                indices, weights = compute_conservative_weights(x_dep, map.grid_meta, map.source_mask)
                for m in 1:16
                    s_idx = indices[m]
                    w = weights[m]
                    if w > 0
                        req = map.demand_map[s_idx]
                        # Scale by 1/demand if over-requested (Mass conservation)
                        scale = req > 1.0 ? (1.0 / req) : 1.0

                        # Accumulate: Weight is shared (0.25 per sub-ray)
                        val_accum += 0.25 * w * scale * source_data[s_idx]
                    end
                end
            end
            target_data[target_idx] = val_accum
        end
    end

    # --- Pass 2: Forward Push (Leakage Correction) ---
    for (k, s_idx) in enumerate(map.leakage_indices)
        val = source_data[s_idx]
        if val != 0.0
            req = map.demand_map[s_idx]
            leftover = (1.0 - req) * val
            x_arr = map.leakage_rays[k]

            indices, weights = compute_conservative_weights(x_arr, map.grid_meta, map.target_mask)
            for m in 1:16
                target_data[indices[m]] += weights[m] * leftover
            end
        end
    end

    return target_data
end

"""
    advect!(target, source, map; type=:conservative)

Conservative CCISL advection for `CartesianMeshField` values.
"""
function advect!(
    target::CartesianMeshField,
    source::CartesianMeshField,
    map::TransportMap;
    type::Symbol=:conservative,
)
    type === :conservative || throw(ArgumentError(
        "TransportMap advect! supports type=:conservative; " *
        "pass a frozen velocity and dt for type=:intensive."))
    target.grid == map.grid_meta || throw(ArgumentError(
        "Conservative advect! target belongs to a different Cartesian grid."))
    source.grid == map.grid_meta || throw(ArgumentError(
        "Conservative advect! source belongs to a different Cartesian grid."))
    advect!(target.data, source.data, map; type=:conservative)
    return target
end

# Utilities

function trace_ray(x::Point{D,T}, velocity_at, dt) where {D,T}
    raw_v1 = velocity_at(x)
    v1 = VectorValue(Float64(raw_v1[1]), Float64(raw_v1[2]))
    all(isfinite, v1) || throw(ArgumentError(
        "Non-finite velocity encountered while tracing from $x."))

    x_mid = x + v1 * dt
    raw_v2 = velocity_at(x_mid)
    v2 = VectorValue(Float64(raw_v2[1]), Float64(raw_v2[2]))
    all(isfinite, v2) || throw(ArgumentError(
        "Non-finite velocity encountered while tracing through $x_mid."))

    x_new = x + 0.5 * (v1 + v2) * dt
    return x_new
end

@inline function compute_conservative_weights!(
    indices_buffer::MVector{16,Int}, weights_buffer::MVector{16,Float64},
    x::Point{2,T}, grid::CartesianGridInfo, allowed=nothing,
) where {T}
    ox, oy = grid.origin
    dx, dy = grid.spacing
    nx, ny = grid.dims

    # Compute normalized coordinates
    ix_raw = 1 + (x[1] - ox) / dx
    iy_raw = 1 + (x[2] - oy) / dy

    # Keep the interpolation point at its TRUE location (guard only NaN, and clamp
    # to the node range [1, nx] so floor() stays sane). We deliberately do NOT
    # relocate near-wall points into the interior band [2, nx-2]: that funnels
    # every near-boundary characteristic onto the second node ring, inflating its
    # conservation demand and starving the boundary ring — producing a spurious
    # ρ≈0 rim one node thick along the domain wall (independent of velocity).
    #
    # Instead we follow Lentine, Grétarsson & Fedkiw (2011), "An unconditionally
    # stable fully conservative semi-Lagrangian method": stencil points that fall
    # outside the domain are "not visible" across the wall, get weight 0, and the
    # remaining (visible) weights are scaled up so Σ w = 1. Because both the
    # backward pull and the forward leakage cast go through this function, this
    # fixes conservation on both passes (w_ij and f_ij in the paper's notation).
    ix_float = clamp(isnan(ix_raw) ? 1.0 : ix_raw, 1.0, Float64(nx))
    iy_float = clamp(isnan(iy_raw) ? 1.0 : iy_raw, 1.0, Float64(ny))

    i = floor(Int, ix_float)
    j = floor(Int, iy_float)

    α = ix_float - i
    β = iy_float - j

    wL_x = quadratic_interpolation_weights(α, :left)
    wR_x = quadratic_interpolation_weights(α, :right)

    wL_y = quadratic_interpolation_weights(β, :left)
    wR_y = quadratic_interpolation_weights(β, :right)

    Wx = SVector(0.5 * wL_x[1], 0.5 * (wL_x[2] + wR_x[1]), 0.5 * (wL_x[3] + wR_x[2]), 0.5 * wR_x[3])
    Wy = SVector(0.5 * wL_y[1], 0.5 * (wL_y[2] + wR_y[1]), 0.5 * (wL_y[3] + wR_y[2]), 0.5 * wR_y[3])

    idx = 1
    sum_w = 0.0
    for (ny_local, wy_val) in enumerate(Wy)
        current_j = j - 1 + (ny_local - 1)
        visible_j = 1 <= current_j <= ny
        clamped_j = clamp(current_j, 1, ny)   # a valid buffer index; weight is 0 when not visible
        for (nx_local, wx_val) in enumerate(Wx)
            current_i = i - 1 + (nx_local - 1)
            visible_i = 1 <= current_i <= nx

            clamped_i = clamp(current_i, 1, nx)
            lin_idx = clamped_i + (clamped_j - 1) * nx
            indices_buffer[idx] = lin_idx

            # Not visible across the domain wall or outside the allowed support →
            # weight 0. Negative quadratic lobes are clipped.
            support_visible = allowed === nothing || allowed[lin_idx]
            w_vis = (visible_i && visible_j && support_visible) ?
                max(0.0, wx_val * wy_val) : 0.0
            weights_buffer[idx] = w_vis
            sum_w += w_vis
            idx += 1
        end
    end

    # Scale up the remaining visible weights so they sum to 1 (Lentine et al.).
    if sum_w > 0.0
        for k in 1:16
            weights_buffer[k] /= sum_w
        end
    else
        fill!(weights_buffer, 0.0)
    end

    return nothing
end

@inline function compute_conservative_weights(x::Point{2,T}, grid::CartesianGridInfo,
                                              allowed=nothing) where {T}
    indices = MVector{16,Int}(undef)
    weights = MVector{16,Float64}(undef)
    compute_conservative_weights!(indices, weights, x, grid, allowed)
    return SVector(indices), SVector(weights)
end

end # module
