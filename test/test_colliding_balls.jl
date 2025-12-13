# =============================================================================
# Tests for Colliding Balls with Prescribed Velocity
# =============================================================================
# This test exercises:
# - Multi-body level set (union of two circles)
# - Topology change handling (collision and merge)
# - Time evolution with advection

# Import functions explicitly for standalone test runs
using EvolvingDomains: reinitialize!, current_time

@testset "Colliding Balls" begin
    # Skip if LevelSetMethods is not available
    if !isdefined(Main, :LevelSetMethods)
        @info "Skipping colliding balls tests - LevelSetMethods not loaded"
        @test_skip "LevelSetMethods required"
        return
    end
    
    # =========================================================================
    # Setup
    # =========================================================================
    domain = (0.0, 2.0, 0.0, 2.0)
    partition = (60, 60)  # Higher resolution for accuracy
    model = CartesianDiscreteModel(domain, partition)
    
    # Two balls close to y=1.0 for faster merge
    center_top = (1.0, 1.2)
    center_bottom = (1.0, 0.8)
    radius = 0.15
    
    # Gap between balls: 1.2 - 0.8 - 2*0.15 = 0.1 units
    
    # Combined level set: union (min) of two circles
    function ϕ0(x)
        d_top = sqrt((x[1] - center_top[1])^2 + (x[2] - center_top[2])^2) - radius
        d_bottom = sqrt((x[1] - center_bottom[1])^2 + (x[2] - center_bottom[2])^2) - radius
        return min(d_top, d_bottom)
    end
    
    # Constant velocity toward y=1.0
    function u_toward_center(x)
        if x[2] > 1.0
            return (0.0, -0.5)  # Top ball moves down
        else
            return (0.0, 0.5)   # Bottom ball moves up
        end
    end
    
    @testset "Initial geometry: two separate balls" begin
        evolver = LevelSetMethodsEvolver(;
            bg_model = model,
            initial_ls = ϕ0,
            velocity = u_toward_center,
            spatial_scheme = :WENO5,
            integrator = :RK3,
            bc = :Neumann
        )
        eg = EvolvingDiscreteGeometry(evolver, model)
        
        ϕ = current_levelset(eg)
        
        # Check initial level set at center (should be positive - gap)
        info = grid_info(eg)
        nx, ny = info.dims
        mid_j = div(ny, 2) + 1
        center_i = div(nx, 2) + 1
        center_idx = center_i + (mid_j - 1) * nx
        
        @test ϕ[center_idx] > 0  # Gap between balls
        
        # Check that we have nodes inside the domain
        inside_count = count(ϕ .< 0)
        @test inside_count > 50
        
        @info "Initial: ϕ at center = $(round(ϕ[center_idx], digits=3)), inside count = $inside_count"
    end
    
    @testset "Collision and merge" begin
        evolver = LevelSetMethodsEvolver(;
            bg_model = model,
            initial_ls = ϕ0,
            velocity = u_toward_center,
            spatial_scheme = :WENO5,
            integrator = :RK3,
            bc = :Neumann
        )
        eg = EvolvingDiscreteGeometry(evolver, model)
        
        # Get center index
        info = grid_info(eg)
        nx, ny = info.dims
        mid_j = div(ny, 2) + 1
        center_i = div(nx, 2) + 1
        center_idx = center_i + (mid_j - 1) * nx
        
        ϕ_initial_center = current_levelset(eg)[center_idx]
        
        # Advance simulation
        Δt = 0.01
        nsteps = 50  # t = 0.5s should be enough
        
        for step in 1:nsteps
            advance!(eg, Δt)
        end
        
        ϕ_final = current_levelset(eg)
        
        @info "Final: t=$(round(current_time(eg), digits=2)), ϕ at center = $(round(ϕ_final[center_idx], digits=3))"
        
        # Level set at center should decrease (balls approaching)
        @test ϕ_final[center_idx] < ϕ_initial_center
        
        # After enough time, balls should merge (negative at center)
        @test ϕ_final[center_idx] < 0.0
        
        # Geometry should still be valid
        @test any(ϕ_final .< 0)  # Has interior
        @test any(ϕ_final .> 0)  # Has exterior
    end
end
