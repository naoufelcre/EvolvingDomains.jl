using EvolvingDomains
using Gridap
using Test

@testset "Reinitialization Call" begin
    domain = (0.0, 1.0, 0.0, 1.0)
    model = CartesianDiscreteModel(domain, (10, 10))
    geom = EvolvingDiscreteGeometry(model, x -> norm(x) - 0.5)

    # test_reinit_env.jl: Just verify it doesn't crash 
    # (reinitialization logic is inside LevelSetMethods)
    @test_nowarn reinitialize!(geom)
    @test geom.last_reinit_step == 0
end
