using Test

@testset "EvolvingDomains.jl" begin
    @testset "Evolving Domains" begin
        include("test_evolving_domains.jl")
    end
    
    @testset "Poisson Evolving" begin
        include("test_poisson_evolving.jl")
    end
    
    @testset "Zalesak Disk" begin
        include("test_zalesak.jl")
    end
end
