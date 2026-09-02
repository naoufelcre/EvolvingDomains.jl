module WENO5

using ...Geometric: CartesianGridInfo, CartesianMeshField, WENO5Cache
using ...Geometric: weno5⁻, weno5⁺
using Gridap.TensorValues: VectorValue

export weno5_step!

# The weno 5 time stepping method.
# For stencils see Geometric/Stencils.jl

# Module-level singleton cache (fallback for direct weno5_step! calls)
const WENO_CACHE = WENO5Cache()

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

function weno5_step!(phi::Vector{Float64}, grid::CartesianGridInfo, 
                     velocity::Vector{VectorValue{2,Float64}}, dt::Float64,
                     cache::Union{WENO5Cache,Nothing}=nothing)
    # Use provided cache or module-level cache
    wcache = isnothing(cache) ? WENO_CACHE : cache
    
    n = length(phi)
    
    # Lazy allocation of temporaries (reallocate if size mismatch)
    if isnothing(wcache.rhs) || length(wcache.rhs) != n
        wcache.rhs = zeros(Float64, n)
        wcache.stage = zeros(Float64, n)
        wcache.phi0 = zeros(Float64, n)
    end
    
    rhs = wcache.rhs
    stage = wcache.stage
    phi0 = wcache.phi0
    
    # Store original phi
    @. phi0 = phi
    
    # RK3 SSP (Shu-Osher form)
    # Stage 1: u1 = u^n + dt*L(u^n)
    phi_wrapper = CartesianMeshField(phi, grid)
    weno5_rhs!(rhs, phi_wrapper, velocity)
    @. stage = phi + dt * rhs  # u1 stored in stage
    
    # Stage 2: u2 = 3/4*u^n + 1/4*u1 + 1/4*dt*L(u1)
    @. phi = stage  # u1 now in phi for L evaluation
    weno5_rhs!(rhs, phi_wrapper, velocity)
    @. stage = 0.75 * phi0 + 0.25 * stage + 0.25 * dt * rhs  # u2 stored in stage
    
    # Stage 3: u^{n+1} = 1/3*u^n + 2/3*u2 + 2/3*dt*L(u2)
    @. phi = stage  # u2 now in phi for L evaluation
    weno5_rhs!(rhs, phi_wrapper, velocity)
    @. phi = (1.0/3.0) * phi0 + (2.0/3.0) * stage + (2.0/3.0) * dt * rhs
    
    return phi
end

end # module
