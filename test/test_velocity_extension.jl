# =============================================================================
# Test: Velocity Extension for FEVelocitySource
# =============================================================================
#
# Verifies that velocity extension correctly extends bulk velocity to void.
#
# =============================================================================

using Test
using EvolvingDomains
using Gridap
using LinearAlgebra: norm

@testset "Velocity Extension" begin
    # Setup: circle geometry with known velocity field
    domain = (0.0, 2.0, 0.0, 2.0)
    n = 20
    model = CartesianDiscreteModel(domain, (n, n))
    
    # Circle centered at (1,1) with radius 0.5
    # ϕ < 0 inside circle (bulk), ϕ > 0 outside (void)
    geom = EvolvingDiscreteGeometry(model, 
        x -> norm(x .- (1.0, 1.0)) - 0.5;
        reinit_freq = 0
    )
    
    # Create a mock FE function that returns constant velocity in bulk
    # In real usage this would be from a Stokes solve
    bulk_velocity = VectorValue(1.0, 0.5)
    
    # Create a simple callable that mimics FEFunction behavior
    struct MockFEFunction
        value::VectorValue{2,Float64}
    end
    (f::MockFEFunction)(x) = f.value
    
    mock_fe = MockFEFunction(bulk_velocity)
    vel = FEVelocitySource(mock_fe)
    
    # Update with geometry to enable extension
    update_velocity!(vel, mock_fe, geom)
    
    @testset "Bulk sampling" begin
        # Point inside circle (bulk) - should return exact FE value
        x_bulk = (1.0, 1.0)  # Center of circle
        v = sample_velocity(vel, x_bulk, 0.0)
        @test v[1] ≈ 1.0
        @test v[2] ≈ 0.5
    end
    
    @testset "Void extension" begin
        # Point just outside interface - should have blended extension
        h = 2.0 / n
        x_near_void = (1.0 + 0.5 + 2*h, 1.0)  # Just outside circle
        v = sample_velocity(vel, x_near_void, 0.0)
        
        # Should be non-zero (extended) but blended
        @test v[1] > 0.0  # Extended velocity present
        @test v[1] < 1.0  # But blended (less than full)
    end
    
    @testset "Far void" begin
        # Point far from interface - should be zero
        x_far = (0.1, 0.1)  # Corner, far from circle
        v = sample_velocity(vel, x_far, 0.0)
        @test v[1] ≈ 0.0 atol=1e-8
        @test v[2] ≈ 0.0 atol=1e-8
    end
    
    @testset "Backward compatibility" begin
        # FEVelocitySource without update_velocity!(geom) should work
        vel_noext = FEVelocitySource(mock_fe)
        v = sample_velocity(vel_noext, (1.0, 1.0), 0.0)
        @test v[1] ≈ 1.0
        @test v[2] ≈ 0.5
    end
end

println("✅ All velocity extension tests passed")
