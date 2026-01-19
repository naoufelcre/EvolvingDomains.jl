module EvolvingDomainsFSI

using EvolvingDomains
using StaticArrays
using LinearAlgebra

include("MaterialDerivative.jl")
include("TensorDerivatives.jl")

export advective_derivative!, material_derivative!
export velocity_gradient!, upper_convected_derivative!

end
