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

function get_velocity end

function is_time_dependent end


struct StaticFunctionVelocity{F} <: AbstractVelocitySource
    func::F
end

@inline get_velocity(v::StaticFunctionVelocity, x, t) = v.func(x)
@inline is_time_dependent(::StaticFunctionVelocity) = false

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

@inline get_velocity(v::TimeDependentVelocity, x, t) = v.func(x, t)
@inline is_time_dependent(::TimeDependentVelocity) = true

"""
    sample_velocity(vel::AbstractVelocitySource, grid, t)

Sample the velocity field onto the grid nodes.
"""
function sample_velocity(vel::AbstractVelocitySource, grid, t)
    nx, ny = grid.dims
    x0, y0 = grid.origin
    dx, dy = grid.spacing

    # We return a flat vector of VectorValue to match grid indexing
    v_field = Vector{VectorValue{2,Float64}}(undef, nx*ny)

    @inbounds for j in 1:ny
        for i in 1:nx
            px = x0 + (i-1)*dx
            py = y0 + (j-1)*dy
            v = get_velocity(vel, SVector(px, py), t)

            # Convert to VectorValue if not already (assuming get_velocity returns iterable)
            v_field[i + (j-1)*nx] = VectorValue(v[1], v[2])
        end
    end
    return v_field
end
