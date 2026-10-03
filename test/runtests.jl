using Test

@testset "EvolvingDomains Tests" begin
    @testset "Bilinear interpolation" begin
        include("TestBilinearInterpolation.jl")
    end
    @testset "Transport weight cache" begin
        include("TestTransportWeightCache.jl")
    end
    @testset "TerminalPlot" begin
        include("TestTerminalPlot.jl")
    end
    @testset "TestConservativeTransport" begin
        include("TestConservativeTransport.jl")
    end
    @testset "TestIntensiveTransport" begin
        include("TestIntensiveTransport.jl")
    end
    @testset "TestGeometryEvolution" begin
        include("TestGeometryEvolution.jl")
    end
    @testset "TestAggregation" begin
        include("TestAggregation.jl")
    end
    @testset "TestMultiComponentTransfer" begin
        include("TestMultiComponentTransfer.jl")
    end
    @testset "TestVisualTransfer" begin
        include("TestVisualTransfer.jl")
    end

    @testset "TestReinitialization" begin
        include("TestReinitialization.jl")
    end

    @testset "TestCurvature" begin
        include("TestCurvature.jl")
    end

    @testset "WENO5 sign selection" begin
        include("TestWENO5.jl")
    end
end
