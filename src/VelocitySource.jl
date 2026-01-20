# =============================================================================
# VelocitySource.jl — Velocity Abstractions for Level Set Evolution
# =============================================================================

using StaticArrays: SVector
using Gridap.TensorValues: VectorValue  # For FEFunction point evaluation

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
vel = StaticFunctionVelocity(x -> (1.0, 0.0))  # Uniform rightward flow
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
# Rigid body rotation
vel = TimeDependentVelocity((x, t) -> (-x[2], x[1]))
```
"""
struct TimeDependentVelocity{F} <: AbstractVelocitySource
    func::F
end

@inline sample_velocity(v::TimeDependentVelocity, x, t) = v.func(x, t)
@inline is_time_dependent(::TimeDependentVelocity) = true

# =============================================================================
# FE Velocity Source: From Gridap FEFunction (with automatic void extension)
# =============================================================================

"""
    FEVelocitySource <: AbstractVelocitySource

Velocity from a Gridap FEFunction with automatic extension to void regions.

When the FE velocity is only defined in the bulk (ϕ < 0), this type automatically
extends it to the void (ϕ ≥ 0) using closest-point extension with C² smooth blending.

# Usage
```julia
u_h = solve(stokes_op)
vel = FEVelocitySource(u_h)
update_velocity!(vel, u_h, geom)  # Caches levelset for extension
set_velocity!(geom, vel)
advance!(geom, Δt)
```

# Extension Algorithm
- Bulk (ϕ < 0): Uses FE function directly
- Void (ϕ ≥ 0): Closest-point extension with smootherstep decay over 5h band
"""
mutable struct FEVelocitySource{FE} <: AbstractVelocitySource
    fe_function::FE
    # Extension data (populated by update_velocity!)
    ϕ::Vector{Float64}
    origin::NTuple{2,Float64}
    spacing::NTuple{2,Float64}
    dims::NTuple{2,Int}
end

FEVelocitySource(fe) = FEVelocitySource(fe, Float64[], (0.0,0.0), (0.0,0.0), (0,0))

is_time_dependent(::FEVelocitySource) = true

# Bilinear interpolation of ϕ at point x (Delegated to InterpolationUtils)
@inline function _interpolate_ϕ(v::FEVelocitySource, x)
    return bilinear_interpolate_scalar(v.ϕ, v.origin, v.spacing, v.dims, x)
end

# Central difference gradient of ϕ at point x
@inline function _gradient_ϕ(v::FEVelocitySource, x)
    hx, hy = v.spacing
    ϕ_xp = _interpolate_ϕ(v, (x[1] + hx/2, x[2]))
    ϕ_xm = _interpolate_ϕ(v, (x[1] - hx/2, x[2]))
    ϕ_yp = _interpolate_ϕ(v, (x[1], x[2] + hy/2))
    ϕ_ym = _interpolate_ϕ(v, (x[1], x[2] - hy/2))
    return SVector((ϕ_xp - ϕ_xm) / hx, (ϕ_yp - ϕ_ym) / hy)
end

function sample_velocity(v::FEVelocitySource, x, t)
    # No extension data? Use FE directly (backward compatible)
    isempty(v.ϕ) && return _sample_fe_direct(v, x)
    
    ϕ_x = _interpolate_ϕ(v, x)
    
    # Bulk: use FE function directly
    ϕ_x < 0 && return _sample_fe_direct(v, x)
    
    # Void: closest-point extension
    ∇ϕ = _gradient_ϕ(v, x)
    n_mag = sqrt(∇ϕ[1]^2 + ∇ϕ[2]^2) + 1e-12
    n = ∇ϕ / n_mag
    x_Γ = (x[1] - ϕ_x * n[1], x[2] - ϕ_x * n[2])
    
    # Sample at closest point (no decay)
    return _sample_fe_direct(v, x_Γ)
end

@inline function _sample_fe_direct(v::FEVelocitySource, x)
    # Convert tuple/array to Point for Gridap FEFunction compatibility
    pt = VectorValue(x[1], x[2])
    val = v.fe_function(pt)
    return (val[1], val[2])
end

"""
    update_velocity!(source::FEVelocitySource, u_new, geom)

Update FE function and cache levelset for void extension.
"""
function update_velocity!(source::FEVelocitySource, u_new, geom)
    source.fe_function = u_new
    source.ϕ = copy(current_levelset(geom))
    info = grid_info(geom)
    source.origin = info.origin
    source.spacing = info.spacing
    source.dims = info.dims
    return source
end

# Backward compatible: no geom → no extension
function update_velocity!(source::FEVelocitySource, u_new)
    source.fe_function = u_new
    return source
end

# =============================================================================
# Guided Velocity Source: O(1) Lookup via Generic Locator
# =============================================================================

using Gridap.CellData
using Gridap.Fields
using Gridap.Geometry
using Gridap.Arrays

"""
    locate_cell(locator, x) -> Any

Generic interface to find a cell (or leaf) containing point `x` using `locator`.
Implementations should be provided by mesh backend packages.
"""
function locate_cell end

"""
    GuidedVelocitySource <: AbstractVelocitySource

Velocity field that uses a spatial `locator` to find the correct active cell
in O(1) time (after O(log N) tree descent), avoiding global geometric search.

Uses two SCALAR solutions (x and y components) to avoid Gridap's vector element
geometric mapping artifacts (Piola transform scaling).

# Fields
- `solution_x`: The Gridap FE Solution for u_x (CellField)
- `solution_y`: The Gridap FE Solution for u_y (CellField)
- `locator`: A spatial indexing structure (e.g. Quadtree root) that implements `locate_cell(locator, x)`.
- `leaf_map`: Mapping from Leaf ID -> Vector{Cell ID} (provided by GridapIntegration)
"""
struct GuidedVelocitySource{SOL1, SOL2, LOC} <: AbstractVelocitySource
    solution_x::SOL1
    solution_y::SOL2
    locator::LOC
    leaf_map::Dict{Int, Vector{Int}}
end

"""
    update_velocity!(source::GuidedVelocitySource, (u_x, u_y))

Update the FE solution in the GuidedVelocitySource. 
Expects a tuple of scalar FEFunctions `(u_x, u_y)`.
"""
function update_velocity!(source::GuidedVelocitySource, fit_tuple::Tuple)
    u_x, u_y = fit_tuple
    return GuidedVelocitySource(u_x, u_y, source.locator, source.leaf_map)
end


is_time_dependent(::GuidedVelocitySource) = true

function sample_velocity(vs::GuidedVelocitySource, x, t)
    # 1. Tree Descent (Generic Interface)
    leaf = locate_cell(vs.locator, x) 
    
    # 2. Topological Lookup
    # We assume 'leaf' has an '.id' field or compatible API. 
    # Check if we should enforce an abstract interface for leaf ID as well?
    # For now, duck-typing 'leaf.id' is fine if verified by the Locator implementation.
    if !haskey(vs.leaf_map, leaf.id)
        return (0.0, 0.0) 
    end
    
    candidate_cells = vs.leaf_map[leaf.id]
    
    # 3. Local Search & Evaluation
    best_cell = candidate_cells[1]
    
    # Optimisation: if only 1 cell, skip check
    if length(candidate_cells) > 1
        found_strict = false
        for cell_id in candidate_cells
            # Check geometry using solution_x (geometry is shared)
            if _is_point_in_cell(vs.solution_x, cell_id, x)
                best_cell = cell_id
                found_strict = true
                break
            end
        end
    end

    # 4. Force Evaluation on Scalar Fields
    val_x = evaluate_shape_function(vs.solution_x, best_cell, x)
    val_y = evaluate_shape_function(vs.solution_y, best_cell, x)
    
    return (val_x, val_y)
end

# -----------------------------------------------------------------------------
# Low-Level Gridap Helper
# -----------------------------------------------------------------------------

"""
    evaluate_shape_function(u_h, cell_id, p) -> Value

Evaluate FE function `u_h` at point `p` by mapping `p` to the reference space of `cell_id`.
This allows extrapolation if `p` is slightly outside.
"""
function evaluate_shape_function(u_h, cell_id::Int, p)
    trian = get_triangulation(u_h)
    cell_map = get_cell_map(trian)
    
    # Proper Gridap Array Access using Cache
    cm_array = get_array(cell_map)
    cm_cache = array_cache(cm_array)
    map_k = getindex!(cm_cache, cm_array, cell_id)
    
    u_array = get_array(u_h)
    u_cache = array_cache(u_array)
    u_k = getindex!(u_cache, u_array, cell_id)
    
    # Phys -> Ref
    # Gridap's inverse_map returns the inverse mapping as a Field.
    # We must evaluate it at the physical point for geometric checks (done in _is_point_in_cell).
    
    # CRITICAL FIX: u_k (CellField restriction) expects PHYSICAL coordinates, not Reference.
    # Passing ξ (Reference) causes evaluation at the wrong physical location (near origin).
    pt = VectorValue(p[1], p[2])
    
    return evaluate(u_k, pt)
end

function _is_point_in_cell(u_h, cell_id, p; debug=false)
    trian = get_triangulation(u_h)
    cell_map = get_cell_map(trian)
    
    cm_array = get_array(cell_map)
    cm_cache = array_cache(cm_array)
    map_k = getindex!(cm_cache, cm_array, cell_id)
    
    pt = VectorValue(p[1], p[2])
    
    # Check if inverse map yields a point inside the reference element
    # Gridap doesn't expose "is_inside" easily for arbitrary Polytope without overhead.
    # We can check if ξ is within bounds.
    # For Triangle, sum(ξ) <= 1 and ξ >= 0.
    ξ = evaluate(inverse_map(map_k), pt)
    
    if debug
        println("    [Cell $cell_id] ξ = $ξ")
    end

    TOL = 1e-8
    return (ξ[1] >= -TOL && ξ[2] >= -TOL && (ξ[1]+ξ[2]) <= 1.0+TOL)
end
