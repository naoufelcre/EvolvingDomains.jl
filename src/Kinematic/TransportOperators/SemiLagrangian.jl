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
end

function TransportMap(geom::EvolvingDiscreteGeometry, velocity::AbstractVelocitySource, dt::Real)
    grid = geom.grid
    meta = grid_info(grid)
    coords = get_node_coordinates(grid)
    n_nodes = num_nodes(grid)

    # Get nodes belonging to current geometry (targets for the backward pull)
    active_current = get_active_indices(geom, :current)

    # 1. Backward Flow (Where active nodes come from)
    dx, dy = meta.spacing
    offsets = SVector(
        VectorValue(-0.25 * dx, -0.25 * dy), VectorValue(0.25 * dx, -0.25 * dy),
        VectorValue(-0.25 * dx, 0.25 * dy), VectorValue(0.25 * dx, 0.25 * dy)
    )

    backward_rays = Vector{SVector{4,Point{2,Float64}}}(undef, length(active_current))
    for (k, i) in enumerate(active_current)
        rays_buffer = MVector{4,Point{2,Float64}}(undef)
        for j in 1:4
            rays_buffer[j] = trace_ray(coords[i] + offsets[j], velocity, -dt)
        end
        backward_rays[k] = SVector(rays_buffer)
    end

    # 2. Conservation Demand (How much mass each source node 'owes' to the targets)
    demand = zeros(Float64, n_nodes)
    for rays in backward_rays
        for x_dep in rays
            indices, weights = compute_conservative_weights(x_dep, meta)
            for m in 1:16
                demand[indices[m]] += 0.25 * weights[m]
            end
        end
    end

    # 3. Leakage Map (Forward rays for mass not 'pulled' by Pass 1)
    leak_idx = Int[]
    leak_rays = Point{2,Float64}[]
    for i in 1:n_nodes
        # If demand < 1.0, some mass at this source node might be left behind
        if demand[i] < 1.0 - 1e-12
            push!(leak_idx, i)
            push!(leak_rays, trace_ray(coords[i], velocity, dt))
        end
    end

    return TransportMap(active_current, backward_rays, demand, leak_idx, leak_rays, meta)
end


"""
    advect!(target_data::Vector{Float64}, source_data::Vector{Float64}, map::TransportMap)

Apply the transport operator defined by `map` to `source_data` and store in `target_data`.
Zero allocations in the inner loops.
"""
function advect!(target_data::Vector{Float64}, source_data::Vector{Float64}, map::TransportMap)
    @assert length(target_data) == length(source_data) "Vector dimension mismatch"
    fill!(target_data, 0.0)

    # --- Pass 1: Backward Pull (with 2x2 Supersampling) ---
    Base.Threads.@threads for k in eachindex(map.active_indices)
        target_idx = map.active_indices[k]
        rays = map.backward_rays[k]
        val_accum = 0.0

        for x_dep in rays
            indices, weights = compute_conservative_weights(x_dep, map.grid_meta)
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
        if val > 1e-15 # Only push if there is mass to push
            req = map.demand_map[s_idx]
            leftover = (1.0 - req) * val
            x_arr = map.leakage_rays[k]

            indices, weights = compute_conservative_weights(x_arr, map.grid_meta)
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
    t_dummy = 0.0
    v1 = get_velocity(velocity, x, t_dummy)

    # Safety: If velocity is invalid, don't move
    if any(isnan, v1)
        return x
    end

    x_mid = x + v1 * dt
    v2 = get_velocity(velocity, x_mid, t_dummy)

    # Safety: If midpoint velocity is invalid, return mid-point or start
    if any(isnan, v2)
        return x_mid
    end

    x_new = x + 0.5 * (v1 + v2) * dt
    return x_new
end

function compute_conservative_weights(x::Point{2,T}, grid::CartesianGridInfo) where {T}
    ox, oy = grid.origin
    dx, dy = grid.spacing
    nx, ny = grid.dims

    # Compute normalized coordinates
    ix_raw = 1 + (x[1] - ox) / dx
    iy_raw = 1 + (x[2] - oy) / dy

    # Clamp coordinates to the region where the 4x4 stencil [i-1, i+2] is valid.
    # For nx nodes, valid i ranges from 2 to nx-2.
    # We clamp ix_float to [2.0, nx-2.0] to ensure floor(ix_float) is in that range
    # and α = ix_float - i is always in [0, 1].
    ix_float = clamp(isnan(ix_raw) ? 2.0 : ix_raw, 2.0, Float64(nx) - 2.0)
    iy_float = clamp(isnan(iy_raw) ? 2.0 : iy_raw, 2.0, Float64(ny) - 2.0)

    i = floor(Int, ix_float)
    j = floor(Int, iy_float)

    α = ix_float - i
    β = iy_float - j

    # Ensure i,j are strictly within bounds (should be redundant with clamp above)
    i = clamp(i, 2, nx - 2)
    j = clamp(j, 2, ny - 2)

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
        clamped_j = clamp(current_j, 1, ny)
        for (nx_local, wx_val) in enumerate(Wx)
            current_i = i - 1 + (nx_local - 1)
            clamped_i = clamp(current_i, 1, nx)

            lin_idx = clamped_i + (clamped_j - 1) * nx
            indices_buffer[idx] = lin_idx

            w_raw = wx_val * wy_val
            w_clipped = max(0.0, w_raw)
            weights_buffer[idx] = w_clipped
            sum_w += w_clipped
            idx += 1
        end
    end

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
