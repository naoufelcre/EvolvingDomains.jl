module SemiLagrangian

using Gridap
using Gridap.Geometry: get_node_coordinates, num_nodes
using Gridap.TensorValues
using LinearAlgebra
using StaticArrays

using ...Geometric: CartesianMeshField, CartesianGridInfo, EvolvingDiscreteGeometry,
    get_active_indices, quadratic_interpolation_weights, grid_info

using ..Kinematic: AbstractVelocitySource, get_velocity

export TransportMap, advect!

# One-time warning flag for NaN velocities — avoids flooding logs in long simulations.
const _NAN_VELOCITY_WARNED = Ref(false)

#This module follows the method described in the paper

#An unconditionally stable fully conservative semi-Lagrangian method
#Lentine, Grétarsson, Fedkiw 2011
# https://doi.org/10.1016/j.jcp.2010.12.036
#
# The important details is that we integrate over the whole neigborhood of the backward ray destination using a supersampling.
# This avoids making too much mass dissapear in the first phase due to diverging characteristics and being obligated to compensate in the phase 𝟚
# Which could lead to unphysical results.

"""
    TransportMap

A discretized representation of the flow between two time steps.
Contains all geometric and kinematic information required for conservative advection.
This object is field-independent and should be reused for all fields advecting
with the same velocity (e.g., components of a strain tensor).

The velocity is a frozen spatial field over the step. `get_velocity` is therefore
evaluated with a dummy time of zero; sample or wrap the velocity at the desired time
before constructing the map.

The geometry must carry both time levels. Materialize its current cut before updating
the level set so that the update can preserve it as `prev_cut`.

the transport map is constructed by tracing rays and calculating conservation weights.

"""
struct TransportMap
    # --- Backward (Pull Phase) ---
    # target_active_idx[k] -> bundle of 4 departure_points
    active_indices::Vector{Int}
    backward_rays::Vector{SVector{4,Point{2,Float64}}}

    # --- Conservation ---
    # source_idx -> total_weight_pulled (Used to scale weights for conservation)
    demand_map::Vector{Float64}

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

function TransportMap(geom::EvolvingDiscreteGeometry, velocity::AbstractVelocitySource, dt::Real)
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
            # OLD: rays_buffer[j] = trace_ray(coords[i] + offsets[j], velocity, -dt)
            x_departure = coords[i] + offsets[j]
            rays_buffer[j] = trace_ray(x_departure, velocity, -dt)
            stationary &= rays_buffer[j] == x_departure
        end
        backward_rays[k] = SVector(rays_buffer)
    end

    if stationary
        for i in eachindex(source_mask)
            source_mask[i] && (stationary &= trace_ray(coords[i], velocity, dt) == coords[i])
        end
    end
    stationary &= source_mask == target_mask

    # 2. Conservation Demand (How much mass each source node 'owes' to the targets)
    demand = zeros(Float64, n_nodes)
    for rays in backward_rays
        for x_dep in rays
            # OLD: indices, weights = compute_conservative_weights(x_dep, meta)
            indices, weights = compute_conservative_weights(x_dep, meta, source_mask)
            for m in 1:16
                # OLD: demand[indices[m]] += 0.25 * weights[m]
                s_idx = indices[m]
                source_mask[s_idx] && (demand[s_idx] += 0.25 * weights[m])
            end
        end
    end

    # 3. Leakage Map (Forward rays for mass not 'pulled' by Pass 1)
    leak_idx = Int[]
    leak_rays = Point{2,Float64}[]
    # OLD: for i in 1:n_nodes
    for i in 1:n_nodes
        source_mask[i] || continue
        # If demand < 1.0, some mass at this source node might be left behind
        if demand[i] < 1.0 - 1e-12
            push!(leak_idx, i)
            push!(leak_rays, trace_ray(coords[i], velocity, dt))
        end
    end

    return TransportMap(active_current, backward_rays, demand, leak_idx, leak_rays, meta,
        source_mask, target_mask, stationary)
end


"""
    advect!(target_data::Vector{Float64}, source_data::Vector{Float64}, map::TransportMap)

Apply the transport operator defined by `map` to `source_data` and store in `target_data`.
Zero allocations in the inner loops.
"""
function advect!(target_data::Vector{Float64}, source_data::Vector{Float64}, map::TransportMap)
    if length(target_data) != length(source_data)
        error("advect!: vector dimension mismatch — target has $(length(target_data)) elements, source has $(length(source_data)).")
    end

    if map.is_identity
        fill!(target_data, 0.0)
        target_data[map.active_indices] .= source_data[map.active_indices]
        return target_data
    end

    fill!(target_data, 0.0)

    # --- Pass 1: Backward Pull (with 2x2 Supersampling) ---
    Base.Threads.@threads for k in eachindex(map.active_indices)
        target_idx = map.active_indices[k]
        rays = map.backward_rays[k]
        val_accum = 0.0

        for x_dep in rays
            # OLD: indices, weights = compute_conservative_weights(x_dep, map.grid_meta)
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

    # --- Pass 2: Forward Push (Leakage Correction) ---
    for (k, s_idx) in enumerate(map.leakage_indices)
        val = source_data[s_idx]
        if abs(val) > 1e-15
            req = map.demand_map[s_idx]
            leftover = (1.0 - req) * val
            x_arr = map.leakage_rays[k]

            # OLD: indices, weights = compute_conservative_weights(x_arr, map.grid_meta)
            indices, weights = compute_conservative_weights(x_arr, map.grid_meta, map.target_mask)
            for m in 1:16
                target_data[indices[m]] += weights[m] * leftover
            end
        end
    end

    return target_data
end

"""
    advect!(target::CartesianMeshField, source::CartesianMeshField, map::TransportMap)

In-place field advection for CartesianMeshField types.
"""
function advect!(target::CartesianMeshField, source::CartesianMeshField, map::TransportMap)
    advect!(target.data, source.data, map)
    return target
end

# Utilities

function trace_ray(x::Point{D,T}, velocity::AbstractVelocitySource, dt) where {D,T}
    # TransportMap represents one frozen velocity snapshot, not v(x,t) evolution.
    v1 = get_velocity(velocity, x, 0.0)

    # Safety: If velocity is invalid, don't move.
    # Warn once so the user knows their velocity function may have a bug.
    if any(isnan, v1)
        if !_NAN_VELOCITY_WARNED[]
            _NAN_VELOCITY_WARNED[] = true
            @warn "trace_ray: NaN velocity detected at $x — treating node as stationary. " *
                  "Check your velocity function for division-by-zero or out-of-domain evaluations. " *
                  "(This warning fires only once per session.)"
        end
        return x
    end

    x_mid = x + v1 * dt
    v2 = get_velocity(velocity, x_mid, 0.0)

    # Safety: If midpoint velocity is invalid, return mid-point.
    if any(isnan, v2)
        if !_NAN_VELOCITY_WARNED[]
            _NAN_VELOCITY_WARNED[] = true
            @warn "trace_ray: NaN midpoint velocity detected at $x_mid — returning midpoint. " *
                  "Check your velocity function for division-by-zero or out-of-domain evaluations. " *
                  "(This warning fires only once per session.)"
        end
        return x_mid
    end

    x_new = x + 0.5 * (v1 + v2) * dt
    return x_new
end

function compute_conservative_weights(x::Point{2,T}, grid::CartesianGridInfo,
                                     allowed=nothing) where {T}
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

    weights_buffer = MVector{16,Float64}(undef)
    indices_buffer = MVector{16,Int}(undef)

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
    if sum_w > 1e-12
        inv_sum = 1.0 / sum_w
        for k in 1:16
            weights_buffer[k] *= inv_sum
        end
    else
        fill!(weights_buffer, 0.0)
    end

    return indices_buffer, weights_buffer
end

end # module
