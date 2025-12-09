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
    # Pre-allocated buffers to avoid allocations per call
    narrow_band::BitVector
    needs_extension::BitVector
end

function NarrowBandExtension(γ, nx, ny, max_iters=50)
    n = nx * ny
    NarrowBandExtension(γ, nx, ny, max_iters, falses(n), falses(n))
end

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

    # Identify nodes needing extension (using pre-allocated buffers)
    narrow_band = ext.narrow_band
    needs_extension = ext.needs_extension
    @. narrow_band = abs(ϕ) < γ
    @. needs_extension = narrow_band & !inside_mask

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
    # Vectorized for performance
    outside_band = @. (abs(ϕ) >= γ) & !inside_mask
    vel[outside_band] .= Ref(SVector(0.0, 0.0))

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
uses the fast DOF-based path with narrow band extension.
Otherwise falls back to point evaluation.
"""
function sample_velocity(source::FEVelocitySource, coords, t)
    fh = source.fe_function
    ext = source.extension

    # Fast path: DOF-based sampling with narrow band extension
    if ext isa NarrowBandExtension && !isnothing(source.levelset_values)
        return _sample_velocity_fast(source, coords)
    end

    # Fallback: point evaluation with error handling for cut geometry
    return _sample_velocity_point_eval_safe(fh, coords, source.levelset_values)
end

"""    _sample_velocity_fast(source::FEVelocitySource, coords)

DOF-based velocity sampling with narrow band extension.

1. Extracts velocity at active triangulation nodes via point evaluation
2. Maps active nodes to background grid indices
3. Applies narrow band extension (∇u · ∇ϕ = 0) for outside nodes
4. Returns velocity array for full background grid

This is more efficient than evaluating at all grid points since it only
evaluates at active nodes (~10-20% of grid typically) and uses the extension
algorithm for the rest.
"""
function _sample_velocity_fast(source::FEVelocitySource, coords)
    fh = source.fe_function
    bg_model = source.bg_model
    ext = source.extension
    ϕ = source.levelset_values
    inside_mask = source.inside_mask
    
    n_bg_nodes = length(coords)
    vel = Vector{SVector{2,Float64}}(undef, n_bg_nodes)
    
    # Initialize all to zero
    fill!(vel, SVector(0.0, 0.0))
    
    # Get active triangulation info
    trian = Gridap.CellData.get_triangulation(fh)
    grid = Gridap.Geometry.get_grid(trian)
    active_nodes = Gridap.Geometry.get_node_coordinates(grid)
    
    # Get background grid parameters for coordinate mapping
    bg_desc = Gridap.Geometry.get_cartesian_descriptor(bg_model)
    origin = bg_desc.origin
    sizes = bg_desc.sizes
    partition = bg_desc.partition
    nx, ny = partition
    
    # Sample velocity at active nodes and map to background grid
    # Use level set values to skip nodes outside physical domain (avoids try-catch)
    for c in active_nodes
        # Convert coordinate to background grid index (column-major)
        gi = round(Int, (c[1] - origin[1]) / sizes[1]) + 1
        gj = round(Int, (c[2] - origin[2]) / sizes[2]) + 1
        bg_idx = gi + (gj - 1) * (nx + 1)
        
        # Bounds check
        if bg_idx < 1 || bg_idx > n_bg_nodes
            continue
        end
        
        # Skip nodes outside physical domain if we have level set info
        # (these nodes can't be evaluated safely in the FE space)
        if !isnothing(ϕ) && ϕ[bg_idx] >= 0
            continue
        end
        
        # Evaluate FE function at this node (should be safe now)
        pt = Point(c[1], c[2])
        val = fh(pt)
        vel[bg_idx] = SVector{2,Float64}(val[1], val[2])
    end
    
    # Apply narrow band extension (∇u · ∇ϕ = 0)
    if ext isa NarrowBandExtension && !isnothing(ϕ) && !isnothing(inside_mask)
        extend_velocity_narrow_band!(vel, ϕ, inside_mask, ext)
    end
    
    return vel
end



"""
    _sample_velocity_point_eval_safe(fh, coords, levelset_values)

Safe point evaluation for cut geometry FE functions.
Returns zero velocity for points outside the active triangulation.
If levelset_values is provided, uses it to skip evaluation for outside points.
"""
function _sample_velocity_point_eval_safe(fh, coords, levelset_values)
    n = length(coords)
    result = Vector{SVector{2,Float64}}(undef, n)
    
    for i in 1:n
        c = coords[i]
        
        # If we have level set values, use them to skip outside points
        if !isnothing(levelset_values) && levelset_values[i] >= 0
            # Outside the physical domain - use zero velocity
            result[i] = SVector{2,Float64}(0.0, 0.0)
            continue
        end
        
        # Try to evaluate at this point
        try
            pt = Point(c...)
            val = fh(pt)
            if val isa VectorValue
                result[i] = SVector{2,Float64}(val[1], val[2])
            else
                result[i] = SVector{2,Float64}(val...)
            end
        catch
            # Point is outside active triangulation - use zero velocity
            result[i] = SVector{2,Float64}(0.0, 0.0)
        end
    end
    
    return result
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
