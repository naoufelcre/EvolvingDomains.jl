using Test
using EvolvingDomains

@testset "EvolvingDomains.jl Tests" begin
    @testset "Reinitialization" begin
        include("test_reinit_env.jl")
    end

    @testset "Transfer Cache" begin
        include("test_cache_transfer.jl")
    end

    @testset "Velocity Extension" begin
        include("test_extension_physics.jl")
    end
end

