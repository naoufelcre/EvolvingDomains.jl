# =============================================================================
# TensorDerivatives.jl
# =============================================================================

using EvolvingDomains: sample_velocity, AbstractVelocitySource, CartesianGridInfo
using StaticArrays: SMatrix, SVector

"""
    velocity_gradient!(out::Vector{<:AbstractMatrix}, u, info; method=:central)

Compute the velocity gradient tensor L = ∇u at all grid nodes.
out[i] will contain [∂u_x/∂x  ∂u_x/∂y; ∂u_y/∂x  ∂u_y/∂y].
"""
function velocity_gradient!(out::Vector{<:AbstractMatrix}, u::AbstractVelocitySource, 
                           info::CartesianGridInfo; method::Symbol=:central)
    nx, ny = info.dims
    hx, hy = info.spacing
    ox, oy = info.origin
    
    # 2x2 identity for size checking
    @assert length(out) == nx * ny
    
    @inbounds for j in 1:ny
        for i in 1:nx
            idx = i + (j-1)*nx
            x = ox + (i-1)*hx
            y = oy + (j-1)*hy
            
            # Central difference for velocity gradient
            # Handle boundaries by one-sided differences
            
            # u(x+h)
            ux_p, uy_p = if i < nx
                sample_velocity(u, (x+hx, y), 0.0)
            else
                sample_velocity(u, (x, y), 0.0) # Fallback (should improve)
            end

            # u(x-h)
            ux_m, uy_m = if i > 1
                sample_velocity(u, (x-hx, y), 0.0)
            else
                sample_velocity(u, (x, y), 0.0)
            end
            
            # u(y+h)
            vx_p, vy_p = if j < ny
                sample_velocity(u, (x, y+hy), 0.0)
            else
                sample_velocity(u, (x, y), 0.0)
            end

            # u(y-h)
            vx_m, vy_m = if j > 1
                sample_velocity(u, (x, y-hy), 0.0)
            else
                sample_velocity(u, (x, y), 0.0)
            end
            
            # Gradients
            # If boundary, use 1st order, else 2nd order central
            dx = (i==1 || i==nx) ? hx : 2*hx
            dux_dx = (ux_p - ux_m) / dx
            duy_dx = (uy_p - uy_m) / dx
            
            dy = (j==1 || j==ny) ? hy : 2*hy
            dux_dy = (vx_p - vx_m) / dy
            duy_dy = (vy_p - vy_m) / dy
            
            # Store in output tensor (row-major: ∇u_ij = ∂u_i/∂x_j)
            # StaticArrays are column-major: constructor takes (col1..., col2...)
            # Col 1: [∂u_x/∂x, ∂u_y/∂x]
            # Col 2: [∂u_x/∂y, ∂u_y/∂y]
            out[idx] = SMatrix{2,2,Float64}(dux_dx, duy_dx, dux_dy, duy_dy)
        end
    end
    return out
end

"""
    upper_convected_derivative!(out, C, C_prev, dt, u, ∇u, info)

Compute the Upper Convected Derivative:
    C_upper = DC/Dt - (∇u)ᵀ⋅C - C⋅(∇u)
    
Where DC/Dt is the material derivative.
"""
function upper_convected_derivative!(out::Vector{<:AbstractMatrix}, 
                                   C::Vector{<:AbstractMatrix}, 
                                   C_prev::Vector{<:AbstractMatrix},
                                   dt::Real,
                                   u::AbstractVelocitySource,
                                   ∇u::Vector{<:AbstractMatrix}, 
                                   info::CartesianGridInfo)
    
    # 1. Compute Material Derivative DC/Dt -> out
    # Reuse 'out' as temporary storage for DC/Dt
    material_derivative!(out, C, C_prev, dt, u, info)
    
    # 2. Add convective terms: - (∇u)ᵀ⋅C - C⋅(∇u)
    @inbounds @simd for i in eachindex(out)
        L = ∇u[i]       # L = ∇u
        C_val = C[i]    # Current C
        
        # Upper convected terms
        # Note: L' is transpose(L)
        # term = - L' * C - C * L
        
        out[i] -= (L' * C_val + C_val * L)
    end
    return out
end

# Helper for just the advective part of Upper Convected (for explicit schemes)
# C_upper_adv = (u⋅∇)C - (∇u)ᵀ⋅C - C⋅(∇u)
function upper_convected_advection!(out::Vector{<:AbstractMatrix},
                                  C::Vector{<:AbstractMatrix},
                                  u::AbstractVelocitySource,
                                  ∇u::Vector{<:AbstractMatrix}, 
                                  info::CartesianGridInfo)

    # 1. Compute Advective Derivative (u⋅∇)C -> out
    advective_derivative!(out, C, u, info)
    
    # 2. Add deformation terms
    @inbounds @simd for i in eachindex(out)
        L = ∇u[i]
        C_val = C[i]
        out[i] -= (L' * C_val + C_val * L)
    end
    return out
end
