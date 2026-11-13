using Test

# Conformance/integration suite. Requires an instantiated environment and the native
# compatibility library; it is intentionally separate from Pkg.test.
if !isdefined(Main, :EvolvingDomains)
    include(joinpath(@__DIR__, "..", "src", "EvolvingDomains.jl"))
end
isdefined(Main.EvolvingDomains, :Field) ||
    error("run this runner in a fresh Julia process with the redesign source")

@testset "EvolvingDomains conformance" begin
    include("conformance/TestNativeCuts.jl")
    include("conformance/TestNativeCuts3D.jl")
    include("conformance/TestCutWorkflow.jl")
    include("conformance/TestCutWorkflow3D.jl")
    include("conformance/TestSubResolutionFilterCut.jl")
end
