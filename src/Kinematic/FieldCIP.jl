module FieldCIP

using ...Geometric: _GridMetrics, _gridMetrics, _checklevel, sameid
using ..CIP: CIPCache, cipLagrangianTransport!

using StaticArrays: SVector

"""
    fieldCipLagrangianTransport!(target, source, velocity, dt; level_id=:init, cache=nothing)

Adapt the level-indexed DataModel Field interface to the full-grid uniform
Cartesian CIP implementation. `velocity` is a nodal SVector vector; no geometry mask
is applied.
"""
function _fieldCipLagrangianTransportVector!(target, source, velocity, dt::Real;
        level_id::Symbol=:init,
        cache::Union{Nothing,CIPCache}=nothing)
    td = target.discretization
    sd = source.discretization
    tg = td.levels[level_id]
    sg = sd.levels[level_id]
    (sameid(tg, sg) && tg.axes == sg.axes) || throw(ArgumentError(
        "CIP source and target must use the same grid identity and axes."))
    info = _gridMetrics(sg)
    source_data = _checklevel(source, level_id, info, "Source Field values")
    source_data isa Array{Float64,length(info.dims)} || throw(ArgumentError(
        "Source Field values must be a contiguous Array; materialize views first."))
    target_data = _checklevel(target, level_id, info, "Target Field values")
    target_data isa Array{Float64,length(info.dims)} || throw(ArgumentError(
        "Target Field values must be a contiguous Array; materialize views first."))
    n = prod(info.dims)
    length(velocity) == n || throw(DimensionMismatch(
        "CIP velocity has $(length(velocity)) values; grid requires $n."))
    length(velocity[1]) == length(info.dims) || throw(ArgumentError(
        "CIP velocity dimension must match the Field grid."))
    Base.mightalias(source_data, target_data) && throw(ArgumentError(
        "CIP source and target values must not alias."))
    source_flat = vec(source_data)
    target_flat = vec(target_data)
    cipLagrangianTransport!(target_flat, source_flat, velocity, info, dt; cache=cache)
    return target
end

fieldCipLagrangianTransport!(target, source, velocity::AbstractVector{<:SVector{N}}, dt::Real;
        kwargs...) where {N} =
    _fieldCipLagrangianTransportVector!(target, source, velocity, dt; kwargs...)


end # module FieldCIP
