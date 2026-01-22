# =============================================================================
# VelocitySource.jl — Velocity Abstractions for Level Set Evolution
# =============================================================================

using StaticArrays: SVector
using Gridap.TensorValues: VectorValue

"""
    AbstractVelocitySource

Base type for velocity field representations.
"""
abstract type AbstractVelocitySource end

"""
    sample_velocity(source::AbstractVelocitySource, x, t) -> Tuple{Float64, Float64}

Sample velocity at position `x` and time `t`.
"""
function sample_velocity end

"""
    is_time_dependent(source::AbstractVelocitySource) -> Bool

Return whether velocity varies with time.
"""
function is_time_dependent end

# =============================================================================
# Static Function Velocity: v(x)
# =============================================================================

"""
    StaticFunctionVelocity <: AbstractVelocitySource

Velocity defined by a static function v(x).

# Example
```julia
vel = StaticFunctionVelocity(x -> (1.0, 0.0))
```
"""
struct StaticFunctionVelocity{F} <: AbstractVelocitySource
    func::F
end

@inline sample_velocity(v::StaticFunctionVelocity, x, t) = v.func(x)
@inline is_time_dependent(::StaticFunctionVelocity) = false

# =============================================================================
# Time-Dependent Velocity: v(x, t)
# =============================================================================

"""
    TimeDependentVelocity <: AbstractVelocitySource

Velocity defined by a time-dependent function v(x, t).

# Example
```julia
vel = TimeDependentVelocity((x, t) -> (-x[2], x[1]))
```
"""
struct TimeDependentVelocity{F} <: AbstractVelocitySource
    func::F
end

@inline sample_velocity(v::TimeDependentVelocity, x, t) = v.func(x, t)
@inline is_time_dependent(::TimeDependentVelocity) = true
