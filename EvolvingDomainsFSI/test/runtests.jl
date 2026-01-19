using Test

@testset "EvolvingDomainsFSI.jl" begin
    include("test_scalar_advection.jl")
    include("test_tensor_advection.jl")
end
