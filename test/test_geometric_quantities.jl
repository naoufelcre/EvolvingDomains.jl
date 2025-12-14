# =============================================================================
# Tests for Geometric Quantities (Curvature/Normals for CutFEM)
# =============================================================================

using Statistics: mean

@testset "Geometric Quantities" begin
    # Skip if LevelSetMethods is not available
    if !isdefined(Main, :LevelSetMethods)
        @info "Skipping geometric quantities tests - LevelSetMethods not loaded"
        @test_skip "LevelSetMethods required"
        return
    end
    
    # =========================================================================
    # Setup: Circle test case
    # =========================================================================
    domain = (0.0, 1.0, 0.0, 1.0)
    partition = (50, 50)
    model = CartesianDiscreteModel(domain, partition)
    
    R = 0.3
    center = (0.5, 0.5)
    ϕ(x) = sqrt((x[1]-center[1])^2 + (x[2]-center[2])^2) - R
    expected_κ = 1/R  # Curvature of circle with radius R
    
    evolver = LevelSetMethodsEvolver(;
        bg_model = model,
        initial_ls = ϕ,
        velocity = x -> (0.0, 0.0),
        spatial_scheme = :WENO5,
        bc = :Neumann
    )
    eg = EvolvingDiscreteGeometry(evolver, model)
    
    @testset "get_cut_cells" begin
        cut_geo = current_cut(eg)
        cut_cells = get_cut_cells(cut_geo)
        
        @test length(cut_cells) > 0
        @test length(cut_cells) < prod(partition)  # Not all cells are cut
        
        # CUT cells should be a reasonable fraction (3-50% for circle)
        ratio = length(cut_cells) / prod(partition)
        @test 0.03 < ratio < 0.5
    end
    
    @testset "expand_to_band" begin
        cut_geo = current_cut(eg)
        cut_cells = get_cut_cells(cut_geo)
        
        # Expand by 1 layer
        band1 = expand_to_band(cut_cells, 1, partition)
        @test length(band1) > length(cut_cells)
        
        # Expand by 2 layers
        band2 = expand_to_band(cut_cells, 2, partition)
        @test length(band2) > length(band1)
        
        # Original cells should be in all bands
        for cell in cut_cells
            @test cell in band1
            @test cell in band2
        end
    end
    
    @testset "cells_to_nodes" begin
        cut_geo = current_cut(eg)
        cut_cells = get_cut_cells(cut_geo)
        
        band_nodes = cells_to_nodes(cut_cells, partition)
        
        @test length(band_nodes) > 0
        @test length(band_nodes) < prod(partition .+ 1)  # Not all nodes
        
        # Each quad cell has 4 nodes, but nodes are shared
        # So: n_cells * 4 >= n_nodes >= n_cells
        @test length(band_nodes) >= length(cut_cells)
    end
    
    @testset "curvature_at_band" begin
        κ_values, band_nodes = curvature_at_band(eg)
        
        @test length(κ_values) == prod(grid_info(eg).dims)
        @test length(band_nodes) > 0
        
        # Band should be a fraction of total nodes
        ratio = length(band_nodes) / length(κ_values)
        @test 0.1 < ratio < 0.5
        
        # Non-zero curvature only in band
        nonzero_count = count(!=(0.0), κ_values)
        @test nonzero_count <= length(band_nodes)
    end
    
    @testset "Circle curvature accuracy" begin
        κ_values, _ = curvature_at_band(eg)
        info = grid_info(eg)
        ϕ_vals = current_levelset(eg)
        
        # Get curvature near interface (|ϕ| < threshold)
        threshold = 0.05
        near_interface = findall(i -> abs(ϕ_vals[i]) < threshold && κ_values[i] != 0, 1:length(κ_values))
        
        @test length(near_interface) > 10
        
        κ_near = κ_values[near_interface]
        mean_κ = mean(κ_near)
        
        # Check curvature accuracy (should be close to 1/R)
        rel_error = abs(mean_κ - expected_κ) / expected_κ
        @test rel_error < 0.05  # 5% tolerance
        
        @info "Circle curvature test: expected κ=$(round(expected_κ, digits=3)), got $(round(mean_κ, digits=3)) (error: $(round(100*rel_error, digits=2))%)"
    end
    
    @testset "Curvature caching" begin
        # First call computes curvature
        eg.curvature_dirty = true
        eg.cached_curvature = nothing
        
        κ1 = interface_curvature(eg)
        @test !eg.curvature_dirty
        @test !isnothing(eg.cached_curvature)
        
        # Second call returns cached
        κ2 = interface_curvature(eg)
        @test κ1 === eg.cached_curvature  # Same object (cached)
        
        # After advance, cache is invalidated
        advance!(eg, 0.001)
        @test eg.curvature_dirty
    end
    
    @testset "interface_normal" begin
        n_Γ = interface_normal(eg)
        # Should return a CellField (from GridapEmbedded)
        @test n_Γ !== nothing
    end
end
