# =============================================================================
# Test: Time-Dependent and FE-Coupled Velocity Sources
# =============================================================================
# This test verifies that velocity sources work correctly:
# 1. Static velocity (backwards compatible)
# 2. Time-dependent function velocity
# 3. FE-coupled velocity (mocked)

using Test
using Gridap
using GridapEmbedded
using LinearAlgebra
using StaticArrays

# Add the local module
push!(LOAD_PATH, joinpath(@__DIR__, "..", "src"))
using EvolvingDomains

# Try to load LevelSetMethods (required for this test)
lsm_available = try
    using LevelSetMethods
    true
catch
    @warn "LevelSetMethods.jl not available, skipping full tests"
    false
end

println("=" ^ 60)
println("Test: Velocity Sources")
println("=" ^ 60)

# Enable geometry visualization (comment out to disable)
setup_test_viz!(enabled=true)

# =============================================================================
# Test 1: VelocitySource Abstract Interface
# =============================================================================
@testset "VelocitySource Interface" begin
    println("\n--- Test 1: VelocitySource Interface ---")
    
    # Static velocity
    u_static(x) = (1.0, 0.0)
    vs_static = StaticFunctionVelocity(u_static)
    @test is_time_dependent(vs_static) == false
    
    # Time-dependent velocity
    u_timedep(x, t) = (-x[2] * cos(t), x[1] * cos(t))
    vs_timedep = TimeDependentVelocity(u_timedep)
    @test is_time_dependent(vs_timedep) == true
    
    # Sample at some coordinates
    coords = [(0.0, 0.0), (1.0, 0.0), (0.0, 1.0), (1.0, 1.0)]
    
    # Static sampling
    vel_static = sample_velocity(vs_static, coords, 0.0)
    @test length(vel_static) == 4
    @test all(v -> v ≈ SVector(1.0, 0.0), vel_static)
    
    # Time-dependent sampling at t=0
    vel_t0 = sample_velocity(vs_timedep, coords, 0.0)
    @test length(vel_t0) == 4
    # At (1,0), t=0: u = (-0, 1) * cos(0) = (0, 1)
    @test vel_t0[2] ≈ SVector(0.0, 1.0)
    
    # Time-dependent sampling at t=π/2
    vel_tpi2 = sample_velocity(vs_timedep, coords, π/2)
    # cos(π/2) ≈ 0, so all velocities should be ≈ 0
    @test all(v -> norm(v) < 1e-10, vel_tpi2)
    
    println("  VelocitySource interface tests passed!")
end

# =============================================================================
# Test 2: Velocity Extension Strategies
# =============================================================================
@testset "Velocity Extension" begin
    println("\n--- Test 2: Velocity Extension Strategies ---")
    
    vel = SVector(1.0, 2.0)
    x_inside = (0.5, 0.5)
    x_outside = (2.0, 2.0)
    domain(x) = x[1]^2 + x[2]^2 < 1.0  # Unit circle
    
    # NoExtension (default)
    ext_none = NoExtension()
    @test extend_velocity(ext_none, vel, x_inside, domain) == vel
    @test extend_velocity(ext_none, vel, x_outside, domain) == vel  # Still uses vel
    
    # ZeroExtension
    ext_zero = ZeroExtension()
    @test extend_velocity(ext_zero, vel, x_inside, domain) == vel
    @test extend_velocity(ext_zero, vel, x_outside, domain) == zero(vel)
    
    # ConstantExtension
    ext_const = ConstantExtension(SVector(0.5, 0.5))
    @test extend_velocity(ext_const, vel, x_inside, domain) == vel
    @test extend_velocity(ext_const, vel, x_outside, domain) == SVector(0.5, 0.5)
    
    println("  Velocity extension tests passed!")
end

# =============================================================================
# Test 3: FE Velocity Source (with mock FE function)
# =============================================================================
@testset "FEVelocitySource" begin
    println("\n--- Test 3: FEVelocitySource ---")
    
    # Create a simple model and FE function
    domain = (0.0, 1.0, 0.0, 1.0)
    partition = (5, 5)
    model = CartesianDiscreteModel(domain, partition)
    
    # Vector-valued FE space
    reffe = ReferenceFE(lagrangian, VectorValue{2,Float64}, 1)
    V = FESpace(model, reffe)
    
    # Interpolate a known velocity field: u = (x, y)
    u_ana(x) = VectorValue(x[1], x[2])
    uh = interpolate_everywhere(u_ana, V)
    
    # Create FE velocity source
    vs_fe = FEVelocitySource(uh, model)
    @test is_time_dependent(vs_fe) == true
    
    # Sample at some coordinates
    coords = [(0.0, 0.0), (0.5, 0.5), (1.0, 1.0)]
    vel = sample_velocity(vs_fe, coords, 0.0)
    
    @test length(vel) == 3
    @test vel[1] ≈ SVector(0.0, 0.0) atol=1e-10
    @test vel[2] ≈ SVector(0.5, 0.5) atol=1e-10
    @test vel[3] ≈ SVector(1.0, 1.0) atol=1e-10
    
    # Test update
    u_new(x) = VectorValue(2*x[1], 2*x[2])
    uh_new = interpolate_everywhere(u_new, V)
    update_velocity!(vs_fe, uh_new)
    
    vel_new = sample_velocity(vs_fe, coords, 0.0)
    @test vel_new[2] ≈ SVector(1.0, 1.0) atol=1e-10
    
    println("  FEVelocitySource tests passed!")
end

# =============================================================================
# Test 4: Integration with LevelSetMethodsEvolver (if available)
# =============================================================================
if lsm_available
    @testset "LevelSetMethodsEvolver Integration" begin
        println("\n--- Test 4: LevelSetMethodsEvolver Integration ---")
        
        # Setup
        domain = (0.0, 1.0, 0.0, 1.0)
        partition = (20, 20)
        model = CartesianDiscreteModel(domain, partition)
        
        # Initial level set (circle at center)
        R = 0.15
        center = (0.5, 0.5)
        ϕ0(x) = R - norm(x .- center)
        
        # Test with static velocity (backwards compatible)
        u_static(x) = (0.5, 0.0)
        evolver_static = LevelSetMethodsEvolver(;
            bg_model = model,
            initial_ls = ϕ0,
            velocity = u_static,
            integrator = :RK3,
            bc = :Neumann
        )
        
        eg_static = EvolvingDiscreteGeometry(evolver_static, model)
        t0 = EvolvingDomains.current_time(eg_static)
        advance!(eg_static, 0.1)
        t1 = EvolvingDomains.current_time(eg_static)
        @test t1 ≈ t0 + 0.1
        println("  Static velocity: OK")
        
        # Test with VelocitySource wrapper
        vs = StaticFunctionVelocity(u_static)
        evolver_vs = LevelSetMethodsEvolver(;
            bg_model = model,
            initial_ls = ϕ0,
            velocity = vs,
            integrator = :RK3,
            bc = :Neumann
        )
        
        eg_vs = EvolvingDiscreteGeometry(evolver_vs, model)
        advance!(eg_vs, 0.1)
        @test EvolvingDomains.current_time(eg_vs) ≈ 0.1
        println("  VelocitySource wrapper: OK")
        
        # Test with time-dependent velocity
        u_rot(x, t) = (-x[2] + 0.5, x[1] - 0.5)  # Rotation around center
        vs_rot = TimeDependentVelocity(u_rot)
        evolver_rot = LevelSetMethodsEvolver(;
            bg_model = model,
            initial_ls = ϕ0,
            velocity = vs_rot,
            integrator = :RK3,
            bc = :Neumann
        )
        
        eg_rot = EvolvingDiscreteGeometry(evolver_rot, model)
        # Advance multiple steps
        for _ in 1:5
            advance!(eg_rot, 0.1)
            show_geometry(eg_rot, "rotation")  # Quick visual sanity check
        end
        @test EvolvingDomains.current_time(eg_rot) ≈ 0.5 atol=1e-10
        println("  Time-dependent velocity: OK")
        
        # Test supports_velocity_update
        @test supports_velocity_update(evolver_rot) == true
        
        println("  LevelSetMethodsEvolver integration tests passed!")
    end
else
    println("\n--- Skipping LevelSetMethodsEvolver tests (LevelSetMethods.jl not available) ---")
    println("Install with: ] add LevelSetMethods")
end

println("\n" * "=" ^ 60)
println("All tests completed!")
println("=" ^ 60)
