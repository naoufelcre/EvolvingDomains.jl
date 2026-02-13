module WENO5

using ...Geometric: CartesianGridInfo, CartesianMeshField
using ...Geometric: weno5⁻, weno5⁺
using Gridap.TensorValues: VectorValue
using StaticArrays

export weno5_step!

#The weno 5 time stepping method.
# For stencils see Geometric/Stencils.jl

function weno5_rhs!(rhs::Vector{Float64}, phi::CartesianMeshField, velocity::Vector{VectorValue{2,Float64}})
    nx, ny = phi.grid.dims

    # Iterate over internal grid points
    @inbounds Base.Threads.@threads for j in 1:ny
        for i in 1:nx
            I = CartesianIndex(i, j)
            idx = i + (j - 1) * nx

            # 1. Reconstruct Gradients
            dx_L = weno5⁻(phi, I, 1)
            dx_R = weno5⁺(phi, I, 1)
            dy_L = weno5⁻(phi, I, 2)
            dy_R = weno5⁺(phi, I, 2)

            # 2. Get Velocity (Pre-sampled)
            v = velocity[idx]

            # 3. Upwind Flux (Hamilton-Jacobi)
            grad_x = (v[1] > 0) ? dx_L : dx_R
            grad_y = (v[2] > 0) ? dy_L : dy_R

            # RHS = - (v ⋅ ∇ϕ)
            rhs[idx] = -(v[1] * grad_x + v[2] * grad_y)
        end
    end
end

function weno5_step!(phi::Vector{Float64}, grid::CartesianGridInfo, velocity::Vector{VectorValue{2,Float64}}, dt::Float64)
    # Forward Euler Time Integration
    n = length(phi)

    phi_wrapper = CartesianMeshField(phi, grid)
    rhs = zeros(Float64, n)

    weno5_rhs!(rhs, phi_wrapper, velocity)

    # Update in place: phi = phi + dt * rhs
    @. phi += dt * rhs

    return phi
end

end # module
