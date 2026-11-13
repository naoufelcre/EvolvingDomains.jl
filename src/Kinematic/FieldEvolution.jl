# Field adapter for the existing 2D uniform-Cartesian WENO5 + SSP-RK3 kernel.

using ..Geometric: _GridMetrics, _gridMetrics, _checklevel
using ..EvolvingDomains: Field, _levelgrid
using StaticArrays: SVector

const _FIELD_LEVEL_WARNING = "Field evolution defaults to the stable :init level; pass level_id=... to select another level explicitly."

function _validateFieldStep(velocity, dt)
    step = Float64(dt)
    isfinite(step) && step >= 0 || throw(ArgumentError(
        "Field evolution requires a finite, nonnegative timestep."))
    all(v -> all(isfinite, v), velocity) || throw(ArgumentError(
        "Field evolution velocity values must be finite."))
    return step
end

function _fieldLevel(field, level_id)
    d = field.discretization
    id = isnothing(level_id) ? :init : level_id
    haskey(getproperty(d, :levels), id) || throw(ArgumentError(
        "Field discretization has no grid level $(repr(id))."))
    if isnothing(level_id) && any(id -> id != :init, keys(getproperty(d, :levels)))
        @warn _FIELD_LEVEL_WARNING maxlog=1
    end
    return d, id
end

"""
    eulerianTransport!(field, velocity, dt; level_id=nothing, cache=nothing) -> field

Advance the selected uniform-Cartesian Float64 Field level (2D or 3D) with WENO5 +
SSP-RK3. Omitting `level_id` selects `:init` and warns when other levels are present. An
optional `WENO5Cache` reuses RK3 workspace. `velocity` must be a sampled nodal vector in node
order (first axis fastest), one `SVector` per node.

The scheme is a five-point reconstruction, so it needs **at least five nodes per axis**.
Shorter meshes are accepted — stencil indices clamp to the edge — but the result is no
longer WENO5, only a finite number computed from a clamped stencil. Use at least five nodes
per axis on every level you transport.
"""
# Exact container types are intentional: the kernel indexes this vector directly;
# views, SubArrays, and Float32 must be materialized by the caller first.
function eulerianTransport!(field::Field, velocity::Vector{SVector{N,Float64}}, dt::Real;
    level_id::Union{Nothing,Symbol}=nothing,
    cache::Union{Nothing,WENO5Cache}=nothing) where {N}
    d, id = _fieldLevel(field, level_id)
    info = _gridMetrics(_levelgrid(d, id))
    length(info.dims) == N || throw(ArgumentError(
        "Field evolution requires a $N-dimensional grid."))
    u = _checklevel(field, id, info, "Level-set Field values")
    u isa Array{Float64,N} || throw(ArgumentError(
        "Field evolution requires a contiguous Array{Float64,$N}; materialize views first."))
    length(velocity) == length(u) || throw(DimensionMismatch(
        "Velocity has $(length(velocity)) nodes; expected $(length(u))."))
    step = _validateFieldStep(velocity, dt)
    iszero(step) && return field

    WENO5.weno5Step!(vec(u), info, velocity, step, cache)
    return field
end
