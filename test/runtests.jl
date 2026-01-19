using Test
using EvolvingDomains

@testset "EvolvingDomains.jl Tests" begin
    @testset "FE Coupling" begin
        include("test_fe_coupling.jl")
    end

    @testset "Mixed Grid" begin
        include("test_mixed_grid.jl")
    end

    @testset "Velocity Extension" begin
        include("test_velocity_extension.jl")
    end




end
