# =============================================================================
# Tests for External Solver API (set_levelset!, set_values!)
# =============================================================================
# These tests require LevelSetMethods.jl to be loaded

@testset "External Solver API" begin
    # Skip if LevelSetMethods is not available
    if !isdefined(Main, :LevelSetMethods)
        @info "Skipping external solver API tests - LevelSetMethods not loaded"
        @test_skip "LevelSetMethods required"
        return
    end
    
    domain = (0.0, 2.0, 0.0, 2.0)
    partition = (20, 20)
    model = CartesianDiscreteModel(domain, partition)
    
    # Circle level set
    center = (1.0, 1.0)
    radius = 0.5
    ϕ0(x) = sqrt((x[1] - center[1])^2 + (x[2] - center[2])^2) - radius
    
    # Simple velocity
    u(x) = (0.1, 0.0)
    
    @testset "set_values! roundtrip" begin
        evolver = LevelSetMethodsEvolver(;
            bg_model = model,
            initial_ls = ϕ0,
            velocity = u
        )
        
        # Get original values
        ϕ_original = copy(current_values(evolver))
        
        # Modify values
        ϕ_modified = ϕ_original .+ 0.1
        set_values!(evolver, ϕ_modified)
        
        # Verify modification took effect
        ϕ_retrieved = current_values(evolver)
        @test ϕ_retrieved ≈ ϕ_modified
        @test !(ϕ_retrieved ≈ ϕ_original)
    end
    
    @testset "set_levelset! invalidates geometry" begin
        evolver = LevelSetMethodsEvolver(;
            bg_model = model,
            initial_ls = ϕ0,
            velocity = u
        )
        eg = EvolvingDiscreteGeometry(evolver, model)
        
        # Access geometry to ensure it's built
        geo1 = current_geometry(eg)
        
        # Modify level set
        ϕ_original = current_levelset(eg)
        ϕ_modified = ϕ_original .+ 0.2
        set_levelset!(eg, ϕ_modified)
        
        # Geometry should be dirty
        @test eg.geometry_dirty == true
        @test isnothing(eg.cached_cut)
        
        # After accessing geometry, it should be rebuilt
        geo2 = current_geometry(eg)
        @test eg.geometry_dirty == false
    end
    
    @testset "grid_info from EvolvingDiscreteGeometry" begin
        evolver = LevelSetMethodsEvolver(;
            bg_model = model,
            initial_ls = ϕ0,
            velocity = u
        )
        eg = EvolvingDiscreteGeometry(evolver, model)
        
        info = grid_info(eg)
        
        @test info isa CartesianGridInfo
        @test info.dims == (21, 21)
        @test info.cells == (20, 20)
        @test info.spacing ≈ (0.1, 0.1)
    end
end
