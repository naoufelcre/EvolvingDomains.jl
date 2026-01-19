# =============================================================================
# Test: FEVelocitySource Coupling with Real Gridap FEFunction
# =============================================================================
#
# End-to-end test verifying that FEVelocitySource correctly couples
# PDE-computed velocity fields to geometry evolution.
#
# =============================================================================

using Test
using EvolvingDomains
using Gridap
using LinearAlgebra: norm

@testset "FEVelocitySource with Gridap FEFunction" begin
    
    @testset "Basic coupling — uniform velocity" begin
        # Setup grid
        domain = (0.0, 2.0, 0.0, 2.0)
        n = 30
        model = CartesianDiscreteModel(domain, (n, n))
        
        # Create FE space for velocity (vector-valued)
        reffe = ReferenceFE(lagrangian, VectorValue{2,Float64}, 1)
        V = FESpace(model, reffe)
        
        # Create a known velocity field: uniform rightward flow
        u_h = interpolate(x -> VectorValue(1.0, 0.0), V)
        
        # Create geometry (circle at center)
        geom = EvolvingDiscreteGeometry(model, 
            x -> norm(x .- (1.0, 1.0)) - 0.3;
            reinit_freq = 5
        )
        
        # Create FEVelocitySource and couple
        vel = FEVelocitySource(u_h)
        update_velocity!(vel, u_h, geom)
        set_velocity!(geom, vel)
        
        # Record initial state
        ϕ_before = copy(current_levelset(geom))
        
        # Advance geometry
        advance!(geom, 0.1)
        
        # Verify geometry moved
        ϕ_after = current_levelset(geom)
        @test ϕ_before ≠ ϕ_after
        
        # More specific: the interface should have shifted right
        # Check that ϕ decreased on the right side of domain (geometry moved into void)
        info = grid_info(geom)
        nx, ny = info.dims
        
        # Sample a point on the right side of the original circle
        # At x=1.4 (right of center), ϕ should decrease as geometry moves right
        right_col = round(Int, 1.4 / info.spacing[1])  # Column index for x≈1.4
        mid_row = ny ÷ 2  # Middle row
        right_idx = right_col + (mid_row - 1) * nx
        
        @test ϕ_after[right_idx] < ϕ_before[right_idx]

        # Visualization
        if !isdir("output")
            mkdir("output")
        end
        writevtk(get_grid(model), "output/fe_coupling_geometry", nodaldata=["phi"=>current_levelset(geom)])
        writevtk(Triangulation(model), "output/fe_coupling_velocity", cellfields=["u"=>u_h])
    end
    
    @testset "Velocity sampling consistency" begin
        # Verify that FEVelocitySource samples correctly from FEFunction
        domain = (0.0, 1.0, 0.0, 1.0)
        model = CartesianDiscreteModel(domain, (10, 10))
        
        reffe = ReferenceFE(lagrangian, VectorValue{2,Float64}, 1)
        V = FESpace(model, reffe)
        
        # Linear velocity field: v = (x, y)
        u_h = interpolate(x -> VectorValue(x[1], x[2]), V)
        
        geom = EvolvingDiscreteGeometry(model, 
            x -> norm(x .- (0.5, 0.5)) - 0.2;
            reinit_freq = 0
        )
        
        vel = FEVelocitySource(u_h)
        update_velocity!(vel, u_h, geom)
        
        # Sample at (0.5, 0.5) — should return ≈ (0.5, 0.5)
        v = sample_velocity(vel, (0.5, 0.5), 0.0)
        @test v[1] ≈ 0.5 atol=0.05
        @test v[2] ≈ 0.5 atol=0.05
        
        # Sample at (0.5, 0.6) — inside circle (bulk), should return ≈ (0.5, 0.6)
        v2 = sample_velocity(vel, (0.5, 0.6), 0.0)
        @test v2[1] ≈ 0.5 atol=0.05
        @test v2[2] ≈ 0.6 atol=0.05
    end
    
    @testset "Update velocity mid-simulation" begin
        # Simulate changing velocity field (as in FSI)
        domain = (0.0, 2.0, 0.0, 2.0)
        model = CartesianDiscreteModel(domain, (20, 20))
        
        reffe = ReferenceFE(lagrangian, VectorValue{2,Float64}, 1)
        V = FESpace(model, reffe)
        
        # Initial: rightward
        u_h1 = interpolate(x -> VectorValue(1.0, 0.0), V)
        
        geom = EvolvingDiscreteGeometry(model, 
            x -> norm(x .- (1.0, 1.0)) - 0.25;
            reinit_freq = 5
        )
        
        vel = FEVelocitySource(u_h1)
        update_velocity!(vel, u_h1, geom)
        set_velocity!(geom, vel)
        
        # First advance (rightward)
        advance!(geom, 0.05)
        ϕ_after_right = copy(current_levelset(geom))
        
        # Change to upward velocity
        u_h2 = interpolate(x -> VectorValue(0.0, 1.0), V)
        update_velocity!(vel, u_h2, geom)
        set_velocity!(geom, vel)
        
        # Second advance (upward)
        advance!(geom, 0.05)
        ϕ_after_up = current_levelset(geom)
        
        # Geometry should have moved — both states different
        @test ϕ_after_right ≠ ϕ_after_up
    end
    
    @testset "Type validation" begin
        domain = (0.0, 1.0, 0.0, 1.0)
        model = CartesianDiscreteModel(domain, (10, 10))
        
        reffe = ReferenceFE(lagrangian, VectorValue{2,Float64}, 1)
        V = FESpace(model, reffe)
        u_h = interpolate(x -> VectorValue(1.0, 0.0), V)
        
        geom = EvolvingDiscreteGeometry(model, 
            x -> norm(x .- (0.5, 0.5)) - 0.2;
            reinit_freq = 0
        )
        
        vel = FEVelocitySource(u_h)
        
        # Valid update should work
        @test_nowarn update_velocity!(vel, u_h, geom)
        
        # Invalid types should throw (after we add validation)
        # These are currently commented out until validation is added:
        # @test_throws ArgumentError update_velocity!(vel, [1.0, 0.0], geom)
        # @test_throws ArgumentError update_velocity!(vel, rand(2, 10), geom)
    end
end

println("✅ All FE coupling tests passed")
