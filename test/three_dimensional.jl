using Test

# Focused 3D contracts and conformance; canonical/advanced cases live in validation/.
if !isdefined(Main, :EvolvingDomains)
    include(joinpath(@__DIR__, "..", "src", "EvolvingDomains.jl"))
end
isdefined(Main.EvolvingDomains, :Field) ||
    error("run this runner in a fresh Julia process; the loaded EvolvingDomains is not the redesign")

@testset "Three-dimensional contracts and conformance" begin
    # Defines `_checkSubcellIndex`, which `TestCutData3D.jl` calls; same order as
    # test/contracts/runtests.jl.
    include("contracts/TestCutData.jl")
    include("contracts/TestCutData3D.jl")
    include("conformance/TestCutCells3DABI.jl")
    include("conformance/TestNativeCuts3D.jl")
    include("contracts/TestFieldEvolution3D.jl")
    include("contracts/TestFieldCIP3D.jl")
    include("contracts/TestCutAggregation3D.jl")
    include("contracts/TestReinitialization3D.jl")
    include("contracts/TestCutTransport3D.jl")
    include("conformance/TestCutWorkflow3D.jl")
    include("conformance/TestUnresolvedFilterCut3D.jl")
end
