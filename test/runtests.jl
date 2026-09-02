using Test

@testset "EvolvingDomains Tests" begin
    @testset "TerminalPlot" begin
        include("TestTerminalPlot.jl")
    end
    @testset "TestConservativeTransport" begin
        include("TestConservativeTransport.jl")
    end
    @testset "TestGeometryEvolution" begin
        include("TestGeometryEvolution.jl")
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
end
