using Test
using EvolvingDomains
using EvolvingDomainsFSI
using Gridap
using LinearAlgebra

@testset "Scalar Advection" begin
    
    # Setup simple geometry and velocity
    domain = (0.0, 1.0, 0.0, 1.0)
    n = 20
    model = CartesianDiscreteModel(domain, (n, n))
    geom = EvolvingDiscreteGeometry(model, x -> norm(x .- (0.5, 0.5)) - 0.2)
    info = grid_info(geom)
    
    # Constant velocity: u = (1, 0)
    u_const = StaticFunctionVelocity(x -> (1.0, 0.0))
    set_velocity!(geom, u_const)
    
    # 1. Test Advective Derivative (u ⋅ ∇)f
    # Field f(x,y) = x
    # Analytical: (u⋅∇)f = 1*1 + 0*0 = 1
    f_x = [c[1] for c in geom.coords_tuples]
    out = zeros(length(f_x))
    
    advective_derivative!(out, f_x, u_const, info)
    
    # Interior nodes should be exactly 1.0 (upwind with dx=1 is exact for linear)
    # Boundaries might differ due to one-sided diff
    # Check center node
    center_idx = length(out) ÷ 2 + n ÷ 2 # approximate center
    @test out[center_idx] ≈ 1.0 atol=1e-5
    
    # 2. Test Material Derivative Df/Dt
    # If the field is just being advected, Df/Dt = 0 implies ∂f/∂t = -(u⋅∇)f
    # Let's test consistent transport:
    # If f(x, t) = x - t (travelling wave), then ∂f/∂t = -1, (u⋅∇)f = 1 => Df/Dt = 0
    
    dt = 0.1
    f_prev = [c[1]     for c in geom.coords_tuples] # f at t=0
    f_curr = [c[1] - dt for c in geom.coords_tuples] # f at t=dt
    
    out_md = zeros(length(f_curr))
    material_derivative!(out_md, f_curr, f_prev, dt, u_const, info)
    
    # Should be zero
    @test abs(out_md[center_idx]) < 1e-10
end
