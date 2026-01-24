using EvolvingDomains
using Gridap
using GridapEmbedded
using Test

@testset "Incremental Safety Check" begin
    # 1. Setup
    n = 10
    model = CartesianDiscreteModel((0, 1, 0, 1), (n, n))
    
    # Static Circle
    sdf(x, t) = norm(VectorValue(x...) - VectorValue(0.5, 0.5)) - 0.25
    geom = EvolvingDiscreteGeometry(model, x->sdf(x, 0.0))
    integrator = IncrementalIntegrator(model, 2)
    
    # 2. Advance Geometry
    # This increments geom.step, but integrator.last_updated_step remains -1 (or 0)
    advance!(geom, 0.1) 
    
    # 3. Intentional Error
    println("Attempting integration without update_integrator!...")
    
    # Should throw ErrorException
    @test_throws ErrorException measure_Ω(integrator, geom)
    @test_throws ErrorException measure_Γ(integrator, geom)
    
    # 4. Fix
    update_integrator!(integrator, geom)
    
    # Should Pass
    dΩ = measure_Ω(integrator, geom)
    @test dΩ !== nothing
    println("Safety check passed: Error thrown as expected when out of sync.")
end
