using Gridap.TensorValues: VectorValue
using Gridap: FESpace, ReferenceFE, lagrangian, FEFunction, gradient, CellField, evaluate
using LinearAlgebra: norm
using StaticArrays: SVector

abstract type AbstractExtensionOperator end

"""
    ClosestPointExtension <: AbstractExtensionOperator

Uses the level set gradient to perform a closest-point extension from the 
bulk domain to the void.
"""
struct ClosestPointExtension <: AbstractExtensionOperator
    grid_info::CartesianGridInfo
    ϕ_values::Vector{Float64}            # Signed distance values on grid
    grad_ϕ::Vector{VectorValue{2,Float64}} # Pre-computed gradients at nodes
end

"""
    ClosestPointExtension(info::CartesianGridInfo, ϕ::Vector{Float64})

Convenience constructor that computes gradients automatically via central FD.
"""
function ClosestPointExtension(info::CartesianGridInfo, ϕ::Vector{Float64})
    grad_ϕ = _compute_fd_gradient(info, ϕ)
    return ClosestPointExtension(info, ϕ, grad_ϕ)
end

"""
Compute gradient via central finite differences on Cartesian grid.
"""
function _compute_fd_gradient(info::CartesianGridInfo, ϕ::Vector{Float64})
    nx, ny = info.dims
    dx, dy = info.spacing
    grad = Vector{VectorValue{2,Float64}}(undef, nx*ny)
    
    Base.Threads.@threads for j in 1:ny
        for i in 1:nx
            idx = i + (j-1)*nx
            # Central differences with one-sided at boundaries
            im, ip = max(1, i-1), min(nx, i+1)
            jm, jp = max(1, j-1), min(ny, j+1)
            
            dϕdx = (ϕ[ip + (j-1)*nx] - ϕ[im + (j-1)*nx]) / ((ip-im)*dx)
            dϕdy = (ϕ[i + (jp-1)*nx] - ϕ[i + (jm-1)*nx]) / ((jp-jm)*dy)
            grad[idx] = VectorValue(dϕdx, dϕdy)
        end
    end
    return grad
end

"""
    extend_field(op::ClosestPointExtension, field::CartesianMeshField)

Extend a field from the bulk to the void using the closest-point mapping:
u_ext(x) = u(x - ϕ(x) * ∇ϕ(x) / |∇ϕ|)

For interior points (ϕ <= 0), it returns the original value.
"""
function extend_field(op::ClosestPointExtension, field::CartesianMeshField)
    nx, ny = op.grid_info.dims
    origin = op.grid_info.origin
    spacing = op.grid_info.spacing
    
    itp = get_interpolator(field)
    u_interpolator = x -> itp(x[1], x[2])
    
    new_data = similar(field.data)
    Base.Threads.@threads for j in 1:ny
        for i in 1:nx
            idx = i + (j-1)*nx
            dist = op.ϕ_values[idx]
            
            if dist <= 0
                # Interior (Bulk): Keep original value
                new_data[idx] = field.data[idx]
            else
                # Exterior (Void): Map to interface and sample
                x_phys = SVector{2,Float64}(origin[1] + (i-1)*spacing[1], origin[2] + (j-1)*spacing[2])
                grad = op.grad_ϕ[idx]
                
                # Normalize gradient for mapping
                gnorm = norm(grad)
                if gnorm > 1e-12
                    # x_cp = x - dist * (grad/|grad|)
                    grad_hat = SVector{2,Float64}(grad[1] / gnorm, grad[2] / gnorm)
                    x_cp = x_phys .- dist .* grad_hat
                else
                    x_cp = x_phys
                end
                
                new_data[idx] = u_interpolator(x_cp)
            end
        end
    end
    
    return CartesianMeshField(new_data, op.grid_info)
end

export AbstractExtensionOperator, ClosestPointExtension, extend_field
