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
using .CIP: CIPCache, cipLagrangianTransport!
using .FieldCIP: fieldCipLagrangianTransport!
include("FieldEvolution.jl")

# `eulerianTransport!` evolves one Field level in place; `lagrangianTransport!` transports
# between a source and target with CIP or a conservative map.

"""
    lagrangianTransport!(target::Field, source::Field, velocity, dt; level_id=:init)

Transport an intensive scalar Field with directionally split CIP. `velocity` is a frozen
nodal vector field. A `CIPCache` carries the interpolation profile between calls.

CIP builds a cubic Hermite profile from one-sided derivatives, so it needs **at least three
nodes per axis**; on two nodes it degenerates to a linear profile between them.

The argument shape selects the operator: pass a `TransportMap` for conservative Cut
transport, or a velocity and timestep for intensive CIP.

Omitting `level_id` selects `:init` and warns once when the Field carries other levels,
exactly like `eulerianTransport!`.
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
    return fieldCipLagrangianTransport!(target, source, velocity, dt;
        level_id=id, cache=cache)
end


end # module
