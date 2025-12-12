# =============================================================================
# Tests for CartesianGridInfo and grid utilities
# =============================================================================

@testset "CartesianGridInfo" begin
    # Create a simple 2D Cartesian model
    domain = (0.0, 2.0, 0.0, 1.0)  # x ∈ [0,2], y ∈ [0,1]
    partition = (20, 10)
    model = CartesianDiscreteModel(domain, partition)
    
    @testset "Construction from model" begin
        info = CartesianGridInfo(model)
        
        @test info.origin == (0.0, 0.0)
        @test info.dims == (21, 11)  # nodes = cells + 1
        @test info.cells == (20, 10)
        @test info.spacing ≈ (0.1, 0.1)
    end
    
    @testset "Grid spacing calculation" begin
        info = CartesianGridInfo(model)
        
        # Verify: spacing = (corner - origin) / cells
        expected_dx = 2.0 / 20
        expected_dy = 1.0 / 10
        @test info.spacing[1] ≈ expected_dx
        @test info.spacing[2] ≈ expected_dy
    end
    
    @testset "Consistency with existing utilities" begin
        info = CartesianGridInfo(model)
        (origin, corner, partition) = cartesian_descriptor(model)
        
        @test info.origin == origin
        @test info.cells == partition
    end
end

@testset "Domain Masks" begin
    domain = (0.0, 2.0, 0.0, 2.0)
    partition = (20, 20)
    model = CartesianDiscreteModel(domain, partition)
    
    # Circle level set: ϕ(x) = |x - center| - radius
    center = (1.0, 1.0)
    radius = 0.5
    ϕ0(x) = sqrt((x[1] - center[1])^2 + (x[2] - center[2])^2) - radius
    
    @testset "domain_mask basic" begin
        coords = get_node_coords_as_vector(model)
        ϕ_values = [ϕ0(c) for c in coords]
        
        mask = domain_mask(ϕ_values)
        
        # Check known points
        # Center should be inside (ϕ < 0)
        # Corners should be outside (ϕ > 0)
        n_inside = count(mask)
        n_total = length(mask)
        
        @test n_inside > 0
        @test n_inside < n_total
        
        # Approximate area check (rough)
        Δx = 2.0 / 20
        approx_inside_area = n_inside * Δx^2
        true_area = π * radius^2
        @test abs(approx_inside_area - true_area) / true_area < 0.2  # 20% tolerance
    end
    
    @testset "narrow_band_mask basic" begin
        coords = get_node_coords_as_vector(model)
        ϕ_values = [ϕ0(c) for c in coords]
        
        bandwidth = 0.2
        band = narrow_band_mask(ϕ_values, bandwidth)
        
        # Band should contain fewer points than total
        @test count(band) < length(band)
        @test count(band) > 0
        
        # All band points should have |ϕ| < bandwidth
        for (i, in_band) in enumerate(band)
            if in_band
                @test abs(ϕ_values[i]) < bandwidth
            end
        end
    end
end
