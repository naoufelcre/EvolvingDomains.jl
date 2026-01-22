using Gridap.TensorValues: VectorValue
using Gridap: FESpace, ReferenceFE, lagrangian, FEFunction, gradient, get_model, CellField, evaluate

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
    extend_field(op::ClosestPointExtension, u_grid::AbstractArray{T,2}) where T

Extend a field `u_grid` from the bulk to the void using the closest-point mapping:
u_ext(x) = u(x - ϕ(x) * ∇ϕ(x) / |∇ϕ|)

For interior points (ϕ <= 0), it returns the original value.
"""
function extend_field(op::ClosestPointExtension, u_grid::AbstractArray{T,2}) where T
    nx, ny = op.grid_info.dims
    origin = op.grid_info.origin
    spacing = op.grid_info.spacing
    
    # interpolator for the mesh-derived field on the grid
    u_interpolator = x -> bilinear_interpolate_scalar(
        vec(u_grid), origin, spacing, op.grid_info.dims, (x[1], x[2])
    )
    
    result = similar(u_grid)
    for j in 1:ny, i in 1:nx
        idx = i + (j-1)*nx
        dist = op.ϕ_values[idx]
        
        if dist <= 0
            # Interior (Bulk): Use the prolonged field directly
            result[i,j] = u_grid[i,j]
        else
            # Exterior (Void): Map to interface and sample
            x_phys = SVector{2,Float64}(origin[1] + (i-1)*spacing[1], origin[2] + (j-1)*spacing[2])
            grad = op.grad_ϕ[idx]
            
            # Normalize gradient for mapping (must handle zero gradient)
            gnorm = norm(grad)
            if gnorm > 1e-12
                # x_cp = x - dist * (grad/|grad|)
                grad_hat = SVector{2,Float64}(grad[1] / gnorm, grad[2] / gnorm)
                x_cp = x_phys .- dist .* grad_hat
            else
                x_cp = x_phys
            end
            
            result[i,j] = u_interpolator(x_cp)
        end
    end
    
    return result
end

export AbstractExtensionOperator, ClosestPointExtension, extend_field
