using Test
using EvolvingDomains
using Gridap

@testset "all-cut component aggregation" begin
    grid = CartesianDiscreteModel((0.0, 1.0, 0.0, 1.0), (1, 1))
    geom = EvolvingDiscreteGeometry([-1.0, 1.0, 1.0, 1.0], grid)
    roots = @test_logs (:warn, r"no interior cell") aggregate_cut_cells(ensure_cut!(geom))
    @test roots == Int32[1]
end
