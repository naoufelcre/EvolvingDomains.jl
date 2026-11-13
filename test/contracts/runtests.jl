using Test

# Lean contract suite: no native CutCells library, no plotting, no long simulations.
if !isdefined(Main, :EvolvingDomains)
    include(joinpath(@__DIR__, "..", "..", "src", "EvolvingDomains.jl"))
end
isdefined(Main.EvolvingDomains, :Field) ||
    error("run this runner in a fresh Julia process; the loaded EvolvingDomains is not the redesign")

@testset "Redesign contracts" begin
    include("TestDataModel.jl")
    include("TestMultilevel.jl")
    include("TestCutData.jl")
    include("TestCutData3D.jl")
    include("TestCutGeometry.jl")
    include("TestCutTransport.jl")
    include("TestCutTransport3D.jl")
    include("TestFieldEvolution.jl")
    include("TestFieldEvolution3D.jl")
    include("TestFieldCIP.jl")
    include("TestFieldCIP3D.jl")
    include("TestCutAggregation.jl")
    include("TestCutAggregation3D.jl")
    include("TestReinitialization3D.jl")
    include("TestClosestPointExtension.jl")
    include("TestExtendCutValues.jl")
    include("TestTopologyFilter.jl")
end
