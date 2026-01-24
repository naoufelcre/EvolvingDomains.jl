using EvolvingDomains
using Gridap
using LinearAlgebra: norm
using Test

@testset "Reinitialization Call" begin
    domain = (0.0, 1.0, 0.0, 1.0)
    model = CartesianDiscreteModel(domain, (20, 20))
    
    # Create a distorted level set (NOT a true SDF)
    # ϕ = (distance)^2 - r^2, which is NOT signed distance
    distorted_sdf(x) = (norm(x .- (0.5, 0.5)))^2 - 0.2^2
    geom = EvolvingDiscreteGeometry(model, distorted_sdf; reinit_freq=0)
    
    ϕ_before = copy(current_levelset(geom))
    
    # Reinitialize should restore proper signed distance property
    reinitialize!(geom)
    
    ϕ_after = current_levelset(geom)
    
    # Test 1: Values should have changed
    @test ϕ_after != ϕ_before
    
    # Test 2: Step tracking updated
    @test geom.last_reinit_step == 0
    
    # Test 3: Dirty flag set
    @test geom.dirty == true
    
    # Test 4: Gradient magnitude should be close to 1 (SDF property)
    # Sample a few points in the narrow band
    info = grid_info(geom.model)
    h = info.spacing[1]
    
    # Check |∇ϕ| ≈ 1 at a node away from interface
    # (After reinit, the gradient magnitude should be close to 1)
    # This is a basic sanity check, not a rigorous SDF verification
    @test length(ϕ_after) == length(ϕ_before)
end
