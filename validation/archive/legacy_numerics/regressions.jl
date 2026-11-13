using Test
using EvolvingDomains

# Run separately against the active dev package; preserve the established checks.
@testset "EvolvingDomains numerical regressions" begin
    include("TestBilinearInterpolation.jl")
    include("TestTransportWeightCache.jl")
    include("TestConservativeTransport.jl")
    include("TestIntensiveTransport.jl")
    include("TestGeometryEvolution.jl")
    geom = TestGeometryEvolution.run_zalesak()
    @test length(TestGeometryEvolution.current_levelset(geom)) == 101 * 101
    @test all(isfinite, TestGeometryEvolution.current_levelset(geom))
    include("TestAggregation.jl")
    include("TestMultiComponentTransfer.jl")
    include("TestVisualTransfer.jl")
    include("TestReinitialization.jl")
    include("TestCurvature.jl")
    include("TestWENO5.jl")
end
