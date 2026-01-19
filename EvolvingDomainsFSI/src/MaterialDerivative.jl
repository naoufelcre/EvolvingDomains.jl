# =============================================================================
# MaterialDerivative.jl
# =============================================================================

using EvolvingDomains: sample_velocity, AbstractVelocitySource, CartesianGridInfo

"""
    advective_derivative!(out, f, u, info; method=:upwind)

Compute the advective term (u ⋅ ∇)f and store it in `out`.
Supports both scalar (Vector{Float64}) and tensor (Vector{Matrix}) fields.

# Methods
- `:upwind` (default): First-order upwind (stable, diffusive)
- `:central`: Second-order central (less diffusive, potentially unstable)
"""
function advective_derivative!(out::AbstractVector, f::AbstractVector, 
                              u::AbstractVelocitySource, info::CartesianGridInfo;
                              method::Symbol=:upwind)
    nx, ny = info.dims
    hx, hy = info.spacing
    ox, oy = info.origin
    
    # Pre-compute strides for efficiency
    sx = 1
    sy = nx
    
    @inbounds for j in 1:ny
        for i in 1:nx
            idx = i + (j-1)*nx
            
            # Grid node position
            x = ox + (i-1)*hx
            y = oy + (j-1)*hy
            
            # Velocity at node
            vx, vy = sample_velocity(u, (x, y), 0.0)
            
            # X-derivative
            df_dx = if method == :upwind
                if vx > 0
                    i > 1 ? (f[idx] - f[idx-sx]) / hx : (f[idx+sx] - f[idx]) / hx
                else
                    i < nx ? (f[idx+sx] - f[idx]) / hx : (f[idx] - f[idx-sx]) / hx
                end
            else # :central
                if i > 1 && i < nx
                    (f[idx+sx] - f[idx-sx]) / (2*hx)
                elseif i == 1
                    (f[idx+sx] - f[idx]) / hx
                else
                    (f[idx] - f[idx-sx]) / hx
                end
            end
            
            # Y-derivative
            df_dy = if method == :upwind
                if vy > 0
                    j > 1 ? (f[idx] - f[idx-sy]) / hy : (f[idx+sy] - f[idx]) / hy
                else
                    j < ny ? (f[idx+sy] - f[idx]) / hy : (f[idx] - f[idx-sy]) / hy
                end
            else # :central
                if j > 1 && j < ny
                    (f[idx+sy] - f[idx-sy]) / (2*hy)
                elseif j == 1
                    (f[idx+sy] - f[idx]) / hy
                else
                    (f[idx] - f[idx-sy]) / hy
                end
            end
            
            out[idx] = vx * df_dx + vy * df_dy
        end
    end
    return out
end

"""
    material_derivative!(out, f, f_prev, dt, u, info; method=:upwind)

Compute the full material derivative Df/Dt = (f - f_prev)/dt + (u ⋅ ∇)f.
"""
function material_derivative!(out::AbstractVector, f::AbstractVector, f_prev::AbstractVector, 
                             dt::Real, u::AbstractVelocitySource, info::CartesianGridInfo;
                             method::Symbol=:upwind)
    # Compute advection term in-place
    advective_derivative!(out, f, u, info; method=method)
    
    # Add time derivative
    @inbounds @simd for i in eachindex(out)
        out[i] += (f[i] - f_prev[i]) / dt
    end
    return out
end
