module Kinematic

using ..EvolvingDomains: Field

# Include the implementation
include("TransportOperators/WENO5.jl")
include("TransportOperators/SemiLagrangian.jl")
include("TransportOperators/CIP.jl")
include("FieldCIP.jl")

using .WENO5: WENO5Cache, weno5Step!
using .SemiLagrangian: TransportMap
import .SemiLagrangian: lagrangianTransport!
using .CIP: CIPCache, _cipLagrangianTransport!
using .FieldCIP: _fieldCipLagrangianTransport!
include("FieldEvolution.jl")


"""
    lagrangianTransport!(target, source, velocity, dt; level_id=nothing, cache=nothing)

Provided a source field, a nodal vector field and a time step, mutate the target Field by
transporting the values in the Lagrangian sense, at the given level. Source and target must share
grid identity and axes, and a departure that leaves the box, crosses its neighbour or has a
non-positive characteristic Jacobian raises `DomainError`.

Algorithm: directionally split cubic-interpolated propagation (CIP).

Note: `cache` can carry the interpolation profile of a previous time step, so that an
unchanged source does not recompute it.

Reference: Takewaki, Nishiguchi & Yabe (1985), "Cubic interpolated pseudo-particle method (CIP)
for solving hyperbolic-type equations", Journal of Computational Physics 61(2) 261–268
(unverified).
"""
function lagrangianTransport!(
    target::Field,
    source::Field,
    velocity,
    dt::Real;
    level_id::Union{Nothing,Symbol}=nothing,
    cache::Union{Nothing,CIPCache}=nothing,
)
    _, id = _fieldLevel(target, level_id)
    return _fieldCipLagrangianTransport!(target, source, velocity, dt;
        level_id=id, cache=cache)
end


end # module
