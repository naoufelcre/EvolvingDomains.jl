# Standalone native checks need only Julia stdlibs and the built C ABI library (2D and 3D).
module EvolvingDomains
include(joinpath(@__DIR__, "..", "src", "DataModel.jl"))
include(joinpath(@__DIR__, "..", "src", "CutData.jl"))
include(joinpath(@__DIR__, "..", "src", "Bindings", "CutCells.jl"))
cut(field::Field, level_id::Symbol; kwargs...) = CutCellsBinding.fromField(field, level_id; kwargs...)
end

include("conformance/TestCutCellsABI.jl")
include("conformance/TestCutCells3DABI.jl")
include("conformance/TestNativeCuts.jl")
