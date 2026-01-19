using Test
using EvolvingDomains
using EvolvingDomainsFSI
using Gridap
using StaticArrays
using LinearAlgebra

@testset "Tensor Advection" begin
    domain = (0.0, 1.0, 0.0, 1.0)
    n = 20
    model = CartesianDiscreteModel(domain, (n, n))
    geom = EvolvingDiscreteGeometry(model, x -> norm(x .- (0.5, 0.5)) - 0.2)
    info = grid_info(geom)
    
    # 1. Velocity Gradient Test
    # Velocity field u = (y, 0) (Simple Shear)
    # ∇u = [0 1; 0 0]
    u_shear = StaticFunctionVelocity(x -> (x[2], 0.0))
    
    grad_u = Vector{SMatrix{2,2,Float64,4}}(undef, length(geom.coords_tuples))
    velocity_gradient!(grad_u, u_shear, info)
    
    mid_idx = length(grad_u) ÷ 2
    L = grad_u[mid_idx]
    
    @test L[1,1] ≈ 0.0 atol=1e-5 # dux/dx
    @test L[1,2] ≈ 1.0 atol=1e-5 # dux/dy
    @test L[2,1] ≈ 0.0 atol=1e-5 # duy/dx
    @test L[2,2] ≈ 0.0 atol=1e-5 # duy/dy

    # 2. Upper Convected Derivative Test
    # For Upper Convected Maxwell fluid in simple shear:
    # If steady state -> C_upper = 0
    # But let's test the operator consistency.
    # C = I (Identity everywhere)
    # u = (y, 0)
    # D C / Dt = 0 (stationary)
    # Upper Convective term = -L'C - CL = -L' - L = -[0 0; 1 0] - [0 1; 0 0] = -[0 1; 1 0]
    # So C_upper should be -[0 1; 1 0]
    
    C = [SMatrix{2,2,Float64}(1.0, 0.0, 0.0, 1.0) for _ in 1:length(geom.coords_tuples)]
    C_prev = copy(C)
    dt = 0.1
    
    out_ucd = zeros(SMatrix{2,2,Float64,4}, length(C))
    upper_convected_derivative!(out_ucd, C, C_prev, dt, u_shear, grad_u, info)
    
    res = out_ucd[mid_idx]
    
    # Expected: -[0 1; 1 0]
    @test res[1,1] ≈ 0.0 atol=1e-5
    @test res[1,2] ≈ -1.0 atol=1e-5
    @test res[2,1] ≈ -1.0 atol=1e-5
    @test res[2,2] ≈ 0.0 atol=1e-5
end
