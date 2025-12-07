# =============================================================================
# Velocity Source Abstraction
# =============================================================================
# This module provides a flexible interface for specifying velocity fields
# used in level set advection. It supports static functions, time-dependent
# functions, and FE-coupled velocities.

using StaticArrays
using Gridap
using Gridap.CellData

# =============================================================================
# Abstract Interface
# =============================================================================

"""
    AbstractVelocitySource

Abstract type for velocity sources used in level set advection.

Implementations must define:
- `sample_velocity(source, coords, t)`: Sample velocity at coordinates and time
- `is_time_dependent(source)`: Whether velocity varies with time
"""
abstract type AbstractVelocitySource end

"""
    sample_velocity(source::AbstractVelocitySource, coords, t) -> Vector{SVector{2,Float64}}

Sample velocity at grid coordinates and time t.
Returns a vector of 2D velocity vectors, one per coordinate.
"""
function sample_velocity(source::AbstractVelocitySource, coords, t)
    error("sample_velocity not implemented for $(typeof(source))")
end

"""
    is_time_dependent(source::AbstractVelocitySource) -> Bool

Return whether velocity varies with time. If false, velocity is sampled
once at construction and cached.
"""
function is_time_dependent(source::AbstractVelocitySource)
    error("is_time_dependent not implemented for $(typeof(source))")
end

# =============================================================================
# Velocity Extension Strategy (Modular)
# =============================================================================

"""
    AbstractVelocityExtension

Abstract type for velocity extension strategies.

Used when evaluating velocity at points outside the physical domain
(e.g., ghost points near cut cells). Implementations define how to
extrapolate or extend the velocity field.
"""
abstract type AbstractVelocityExtension end

"""
    extend_velocity(ext::AbstractVelocityExtension, vel, x, domain_indicator)

Extend velocity `vel` at point `x`. The `domain_indicator` function
returns whether a point is inside the physical domain.
"""
function extend_velocity(ext::AbstractVelocityExtension, vel, x, domain_indicator)
    error("extend_velocity not implemented for $(typeof(ext))")
end

"""
    NoExtension <: AbstractVelocityExtension

Default extension: use velocity as-is (assumes all sample points are valid).
"""
struct NoExtension <: AbstractVelocityExtension end

extend_velocity(::NoExtension, vel, x, domain_indicator) = vel

"""
    ZeroExtension <: AbstractVelocityExtension

Zero extension: return zero velocity for points outside domain.
"""
struct ZeroExtension <: AbstractVelocityExtension end

function extend_velocity(::ZeroExtension, vel, x, domain_indicator)
    return domain_indicator(x) ? vel : zero(vel)
end

"""
    ConstantExtension <: AbstractVelocityExtension

Constant extension: use a fixed velocity for points outside domain.
"""
struct ConstantExtension{V} <: AbstractVelocityExtension
    value::V
end

function extend_velocity(ext::ConstantExtension, vel, x, domain_indicator)
    return domain_indicator(x) ? vel : ext.value
end

"""
    NarrowBandExtension <: AbstractVelocityExtension

Osher-Sethian narrow band velocity extension.

Extends velocity from inside nodes to a narrow band around the interface
using constant extrapolation in the normal direction (∇u · ∇ϕ = 0).

# Fields
- `bandwidth`: Width of narrow band γ (should be k × Δx where k ≥ stencil width)
- `nx, ny`: Grid dimensions for 2D neighbor lookup
- `max_iters`: Maximum extension iterations (default: 50)

# Notes
Velocity outside the narrow band is set to zero, as it doesn't affect
level set evolution near the interface.
"""
struct NarrowBandExtension <: AbstractVelocityExtension
    bandwidth::Float64
    nx::Int
    ny::Int
    max_iters::Int
end

NarrowBandExtension(γ, nx, ny) = NarrowBandExtension(γ, nx, ny, 50)

"""
    extend_velocity_narrow_band!(vel, ϕ, inside_mask, ext::NarrowBandExtension)

Extend velocity from inside nodes to narrow band using neighbor propagation.

Uses iterative relaxation: each outside node in the narrow band copies
velocity from its neighbor with smallest |ϕ| that already has velocity.

This approximates the PDE extension ∇u · ∇ϕ = 0 (constant in normal direction).
"""
function extend_velocity_narrow_band!(vel::Vector{SVector{2,Float64}}, 
                                       ϕ::Vector{Float64},
                                       inside_mask::BitVector,
                                       ext::NarrowBandExtension)
    n = length(vel)
    γ = ext.bandwidth
    nx, ny = ext.nx, ext.ny
    
    # Identify nodes needing extension
    narrow_band = [abs(ϕ[i]) < γ for i in 1:n]
    needs_extension = [narrow_band[i] && !inside_mask[i] for i in 1:n]
    
    # 2D indexing helper (assumes row-major: x varies fastest)
    idx(i, j) = (j - 1) * nx + i
    
    for iter in 1:ext.max_iters
        changed = false
        
        for j in 1:ny, i in 1:nx
            k = idx(i, j)
            if !needs_extension[k]
                continue
            end
            
            # Collect valid neighbors
            neighbors = Int[]
            if i > 1 push!(neighbors, idx(i-1, j)) end
            if i < nx push!(neighbors, idx(i+1, j)) end
            if j > 1 push!(neighbors, idx(i, j-1)) end
            if j < ny push!(neighbors, idx(i, j+1)) end
            
            # Find neighbor with smallest |ϕ| that has velocity
            best_k = -1
            best_phi = Inf
            for nk in neighbors
                has_velocity = inside_mask[nk] || !needs_extension[nk]
                if has_velocity && abs(ϕ[nk]) < best_phi
                    best_phi = abs(ϕ[nk])
                    best_k = nk
                end
            end
            
            # Copy velocity from best neighbor
            if best_k > 0 && vel[k] != vel[best_k]
                vel[k] = vel[best_k]
                needs_extension[k] = false  # Mark as done
                changed = true
            end
        end
        
        if !changed
            break  # Converged
        end
    end
    
    # Zero velocity for nodes that are:
    # 1. Outside the narrow band (|ϕ| >= γ), AND
    # 2. Outside the physical domain (ϕ >= 0)
    # Inside nodes always keep their velocity (even if outside band)
    for i in 1:n
        if abs(ϕ[i]) >= γ && !inside_mask[i]
            vel[i] = SVector(0.0, 0.0)
        end
    end
    
    return vel
end

# =============================================================================
# Static Function Velocity
# =============================================================================

"""
    StaticFunctionVelocity{F, E<:AbstractVelocityExtension} <: AbstractVelocitySource

Velocity source from a static function `u(x) -> velocity`.
The velocity is time-independent and sampled once per call.

# Fields
- `func`: Function `x -> SVector{2,Float64}` or tuple
- `extension`: Velocity extension strategy

# Example
```julia
u(x) = (1.0, 0.0)  # Constant rightward flow
vel_source = StaticFunctionVelocity(u)
```
"""
struct StaticFunctionVelocity{F, E<:AbstractVelocityExtension} <: AbstractVelocitySource
    func::F
    extension::E
end

StaticFunctionVelocity(f) = StaticFunctionVelocity(f, NoExtension())

is_time_dependent(::StaticFunctionVelocity) = false

function sample_velocity(source::StaticFunctionVelocity, coords, t)
    f = source.func
    return [SVector{2,Float64}(f(c)...) for c in coords]
end

# =============================================================================
# Time-Dependent Function Velocity
# =============================================================================

"""
    TimeDependentVelocity{F, E<:AbstractVelocityExtension} <: AbstractVelocitySource

Velocity source from a time-dependent function `u(x, t) -> velocity`.

# Fields
- `func`: Function `(x, t) -> SVector{2,Float64}` or tuple
- `extension`: Velocity extension strategy

# Example
```julia
u(x, t) = (-x[2], x[1]) * cos(t)  # Oscillating rotation
vel_source = TimeDependentVelocity(u)
```
"""
struct TimeDependentVelocity{F, E<:AbstractVelocityExtension} <: AbstractVelocitySource
    func::F
    extension::E
end

TimeDependentVelocity(f) = TimeDependentVelocity(f, NoExtension())

is_time_dependent(::TimeDependentVelocity) = true

function sample_velocity(source::TimeDependentVelocity, coords, t)
    f = source.func
    return [SVector{2,Float64}(f(c, t)...) for c in coords]
end

# =============================================================================
# FE-Coupled Velocity Source (with Fast DOF Path)
# =============================================================================

"""
    FEVelocitySource{V, M, E<:AbstractVelocityExtension} <: AbstractVelocitySource

Velocity source wrapping a Gridap CellField (FEFunction) with optimized sampling.

Supports two sampling modes:
1. **Point evaluation** (fallback): Calls `fh(Point(x,y))` for each node
2. **Fast DOF path** (with NarrowBandExtension): Direct DOF array access + extension

# Fields
- `fe_function`: Current velocity CellField (mutable for in-place updates)
- `bg_model`: Background CartesianDiscreteModel
- `extension`: Velocity extension strategy
- `levelset_values`: Level set values at grid nodes (for narrow band)
- `inside_mask`: Nodes inside physical domain (ϕ < 0)

# Example
```julia
# Basic usage (point evaluation)
vel_source = FEVelocitySource(velocity_fh, model)

# Fast path with narrow band extension
nx, ny = partition .+ 1
γ = 6 * Δx  # 6 cells bandwidth for WENO5
ext = NarrowBandExtension(γ, nx, ny)
vel_source = FEVelocitySource(velocity_fh, model, ext)

# Update level set values (needed for narrow band)
update_levelset!(vel_source, ϕ_values)

# When velocity is updated from new FE solve
update_velocity!(vel_source, new_velocity_fh)
```
"""
mutable struct FEVelocitySource{V, M, E<:AbstractVelocityExtension} <: AbstractVelocitySource
    fe_function::V
    bg_model::M
    extension::E
    levelset_values::Union{Nothing, Vector{Float64}}
    inside_mask::Union{Nothing, BitVector}
end

# Constructors
FEVelocitySource(fh, model) = FEVelocitySource(fh, model, NoExtension(), nothing, nothing)
FEVelocitySource(fh, model, ext) = FEVelocitySource(fh, model, ext, nothing, nothing)

is_time_dependent(::FEVelocitySource) = true  # Always re-sample when called

"""
    update_levelset!(source::FEVelocitySource, ϕ::Vector{Float64})

Update the level set values for the narrow band extension.
Also recomputes the inside_mask.
"""
function update_levelset!(source::FEVelocitySource, ϕ::Vector{Float64})
    source.levelset_values = ϕ
    source.inside_mask = BitVector(ϕ .< 0)
    return source
end

"""
    sample_velocity(source::FEVelocitySource, coords, t)

Sample velocity at grid coordinates.

If `extension` is `NarrowBandExtension` and `levelset_values` is set,
uses the fast DOF path with narrow band extension.
Otherwise falls back to point evaluation.
"""
function sample_velocity(source::FEVelocitySource, coords, t)
    fh = source.fe_function
    ext = source.extension
    
    # Fast path: DOF + narrow band extension
    if ext isa NarrowBandExtension && !isnothing(source.levelset_values)
        return _sample_velocity_fast(source, coords)
    end
    
    # Fallback: point evaluation
    return _sample_velocity_point_eval(fh, coords)
end

"""
    _sample_velocity_fast(source::FEVelocitySource, coords)

Fast sampling using direct DOF access + narrow band extension.
DOF ordering verified to be interleaved: [u1_x, u1_y, u2_x, u2_y, ...]
"""
function _sample_velocity_fast(source::FEVelocitySource, coords)
    # Get DOF values directly (avoids point evaluation overhead)
    dofs = get_free_dof_values(source.fe_function)
    n = length(coords)
    
    # Extract velocity vectors from interleaved DOFs
    vel = Vector{SVector{2,Float64}}(undef, n)
    for i in 1:n
        vel[i] = SVector(dofs[2i-1], dofs[2i])
    end
    
    # Apply narrow band extension
    extend_velocity_narrow_band!(vel, source.levelset_values, source.inside_mask, 
                                  source.extension)
    
    return vel
end

"""
    _sample_velocity_point_eval(fh, coords)

Fallback sampling via point evaluation (slower but always works).
"""
function _sample_velocity_point_eval(fh, coords)
    return map(coords) do c
        pt = Point(c...)
        val = fh(pt)
        if val isa VectorValue
            SVector{2,Float64}(val[1], val[2])
        else
            SVector{2,Float64}(val...)
        end
    end
end

"""
    update_velocity!(source::FEVelocitySource, new_fh)

Update the wrapped FE velocity function with a new solution.
"""
function update_velocity!(source::FEVelocitySource, new_fh)
    source.fe_function = new_fh
    return source
end

