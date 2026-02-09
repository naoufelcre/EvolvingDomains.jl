module SemiLagrangian

using Gridap
using Gridap.Geometry: get_node_coordinates
using Gridap.TensorValues
using Gridap: CartesianDiscreteModel
using LinearAlgebra
using StaticArrays

using ...Geometric
using ...Geometric: CartesianMeshField, CartesianGridInfo, EvolvingDiscreteGeometry, get_active_indices
using ...Geometric: quadratic_interpolation_weights

using ..Kinematic: AbstractVelocitySource, get_velocity

export RayMaps, compute_raymaps, advect_conservative

"""
    struct RayMaps{T}

Stores the backward and forward flow maps.
"""
struct RayMaps{T}
    # Maps index k (in active_current) -> x_prev (departure point)
    backward::Vector{Point{2,T}} 
end

"""
    compute_raymaps(geom::EvolvingDiscreteGeometry, velocity::AbstractVelocitySource, dt::Real) -> RayMaps
"""
function compute_raymaps(geom::EvolvingDiscreteGeometry, velocity::AbstractVelocitySource, dt::Real)
    grid = geom.grid
    active_current = get_active_indices(geom, :current)
    coords = get_node_coordinates(grid)
    
    backward_map = map(active_current) do I
        x_target = coords[I]
        trace_ray(x_target, velocity, -dt) 
    end

    return RayMaps(backward_map)
end

function trace_ray(x::Point{D,T}, velocity::AbstractVelocitySource, dt) where {D,T}
    t_dummy = 0.0 
    v1 = get_velocity(velocity, x, t_dummy)
    x_mid = x + v1 * dt
    v2 = get_velocity(velocity, x_mid, t_dummy)
    x_new = x + 0.5 * (v1 + v2) * dt
    return x_new
end

# =============================================================================
# Conservative Weight Calculation
# =============================================================================

function compute_conservative_weights(x::Point{2,T}, grid::CartesianGridInfo) where T
    ox, oy = grid.origin
    dx, dy = grid.spacing
    nx, ny = grid.dims
    
    ix_float = 1 + (x[1] - ox) / dx
    iy_float = 1 + (x[2] - oy) / dy
    
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
    
    weights_buffer = MVector{16, Float64}(undef)
    indices_buffer = MVector{16, Int}(undef)
    
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
        for k in 1:16; weights_buffer[k] *= inv_sum; end
    else
        fill!(weights_buffer, 0.0)
    end
    
    return indices_buffer, weights_buffer
end

# =============================================================================
# Main Advection Routine
# =============================================================================

function advect_conservative(geom::EvolvingDiscreteGeometry, field::CartesianMeshField, velocity::AbstractVelocitySource, dt::Real)
    grid_meta = field.grid
    coords = get_node_coordinates(geom.grid)
    n_nodes = length(field.data)
    
    new_data = zeros(Float64, n_nodes)
    weights_requested = zeros(Float64, n_nodes)
    
    # Precompute Raymaps for the active domain (Targets)
    raymaps = compute_raymaps(geom, velocity, dt)
    active_current = get_active_indices(geom, :current)
    
    # Pass 1: Accumulate Demand (Pull phase)
    # We iterate over TARGETS to see what they want to pull from SOURCES
    for (k, target_idx) in enumerate(active_current)
        x_dep = raymaps.backward[k]
        indices, weights = compute_conservative_weights(x_dep, grid_meta)
        for m in 1:16
            weights_requested[indices[m]] += weights[m]
        end
    end
    
    # Pass 2: Main Transport (Pull)
    for (k, target_idx) in enumerate(active_current)
        x_dep = raymaps.backward[k]
        indices, weights = compute_conservative_weights(x_dep, grid_meta)
        val_accum = 0.0
        for m in 1:16
            src_idx = indices[m]
            w = weights[m]
            if w > 0
                req = weights_requested[src_idx]
                # Scale by 1/req if req > 1 (Mass conservation constraint)
                scale_factor = req > 1.0 ? (1.0 / req) : 1.0
                val_accum += w * scale_factor * field.data[src_idx]
            end
        end
        new_data[target_idx] = val_accum
    end
    
    # Pass 3: Forward Correction (Push Leftovers)
    # CRITICAL: We iterate over ALL nodes that have mass, not just the active ones.
    # This prevents mass from "leaking" out of the active set and being forgotten.
    for src_idx in 1:n_nodes
        val = field.data[src_idx]
        if val > 1e-15 # Mass present
            req = weights_requested[src_idx]
            if req < 1.0 - 1e-12
                # Some mass was not pulled, so we push it forward
                leftover_mass = (1.0 - req) * val
                x_source = coords[src_idx]
                x_arrival = trace_ray(x_source, velocity, dt)
                
                indices, weights = compute_conservative_weights(x_arrival, grid_meta)
                for m in 1:16
                    dest_idx = indices[m]
                    new_data[dest_idx] += weights[m] * leftover_mass
                end
            end
        end
    end
    
    return CartesianMeshField(new_data, grid_meta)
end

end # module