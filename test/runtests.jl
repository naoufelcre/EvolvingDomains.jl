using Test

@testset "EvolvingDomains Tests" begin
    @testset "TestConservativeTransport" begin
        include("TestConservativeTransport.jl")
    end
    @testset "TestDumbellParabolic" begin
        include("TestDumbellParabolic.jl")
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
end
