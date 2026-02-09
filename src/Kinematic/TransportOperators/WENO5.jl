module WENO5

using ...Geometric: CartesianGridInfo, CartesianMeshField
using ...Geometric: D⁺, D⁻
using Gridap.TensorValues: VectorValue
using StaticArrays

export weno5_step!

# =============================================================================
# Reconstruction Kernels
# =============================================================================

@inline function _weno5_weights(v1, v2, v3, v4, v5)
    ϵ = 1e-6
    # Smoothness indicators (Jiang & Peng)
    S1 = 13/12 * (v1 - 2*v2 + v3)^2 + 1/4 * (v1 - 4*v2 + 3*v3)^2
    S2 = 13/12 * (v2 - 2*v3 + v4)^2 + 1/4 * (v2 - v4)^2
    S3 = 13/12 * (v3 - 2*v4 + v5)^2 + 1/4 * (3*v3 - 4*v4 + v5)^2

    # Non-linear weights
    α1 = 0.1 / (S1 + ϵ)^2
    α2 = 0.6 / (S2 + ϵ)^2
    α3 = 0.3 / (S3 + ϵ)^2
    sum_α = α1 + α2 + α3

    return α1/sum_α, α2/sum_α, α3/sum_α
end

function weno5⁻(ϕ::CartesianMeshField, I, dim)
    # Left-biased stencil (for v > 0)
    offset = CartesianIndex(ntuple(k -> k == dim ? 1 : 0, 2))

    v1 = D⁻(ϕ, I - 2*offset, dim)
    v2 = D⁻(ϕ, I - offset,   dim)
    v3 = D⁻(ϕ, I,            dim)
    v4 = D⁻(ϕ, I + offset,   dim)
    v5 = D⁻(ϕ, I + 2*offset, dim)

    w1, w2, w3 = _weno5_weights(v1, v2, v3, v4, v5)

    q1 = v1/3 - 7*v2/6 + 11*v3/6
    q2 = -v2/6 + 5*v3/6 + v4/3
    q3 = v3/3 + 5*v4/6 - v5/6

    return w1*q1 + w2*q2 + w3*q3
end

function weno5⁺(ϕ::CartesianMeshField, I, dim)
    # Right-biased stencil (for v < 0)
    offset = CartesianIndex(ntuple(k -> k == dim ? 1 : 0, 2))

    v1 = D⁺(ϕ, I + 2*offset, dim)
    v2 = D⁺(ϕ, I + offset,   dim)
    v3 = D⁺(ϕ, I,            dim)
    v4 = D⁺(ϕ, I - offset,   dim)
    v5 = D⁺(ϕ, I - 2*offset, dim)

    w1, w2, w3 = _weno5_weights(v1, v2, v3, v4, v5)

    q1 = v1/3 - 7*v2/6 + 11*v3/6
    q2 = -v2/6 + 5*v3/6 + v4/3
    q3 = v3/3 + 5*v4/6 - v5/6

    return w1*q1 + w2*q2 + w3*q3
end

# =============================================================================
# Time Stepping
# =============================================================================

function weno5_rhs!(rhs::Vector{Float64}, phi::CartesianMeshField, velocity::Vector{VectorValue{2,Float64}})
    nx, ny = phi.grid.dims

    # Iterate over internal grid points
    @inbounds for j in 1:ny
        for i in 1:nx
            I = CartesianIndex(i, j)
            idx = i + (j-1)*nx

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
            rhs[idx] = - (v[1]*grad_x + v[2]*grad_y)
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
