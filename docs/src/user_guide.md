# User Guide

A practical guide to implementing moving domain simulations with FE-coupled level set advection.

## Architecture Overview

```
┌─────────────────────────────────────────────────────────────────┐
│                     Your Simulation Code                        │
│  ┌──────────────┐  ┌──────────────┐  ┌──────────────────────┐  │
│  │ Physics      │  │ Geometry     │  │ Time Loop            │  │
│  │ (Stokes,     │→ │ Evolution    │→ │ advance!, solve,     │  │
│  │  Poisson)    │  │              │  │ repeat               │  │
│  └──────────────┘  └──────────────┘  └──────────────────────┘  │
└─────────────────────────────────────────────────────────────────┘
                              ↓
┌─────────────────────────────────────────────────────────────────┐
│                     EvolvingDomains Package                     │
│  ┌──────────────────────────────────────────────────────────┐  │
│  │ EvolvingDiscreteGeometry                                  │  │
│  │ • Manages geometry state                                  │  │
│  │ • advance!(eg, Δt) → moves interface                      │  │
│  │ • current_geometry(eg) → DiscreteGeometry for Gridap     │  │
│  └──────────────────────────────────────────────────────────┘  │
│                              ↓                                  │
│  ┌──────────────────────────────────────────────────────────┐  │
│  │ AbstractLevelSetEvolver (backend interface)               │  │
│  │ • LevelSetMethodsEvolver (uses LevelSetMethods.jl)       │  │
│  │ • evolve!, current_values, reinitialize!                 │  │
│  └──────────────────────────────────────────────────────────┘  │
│                              ↓                                  │
│  ┌──────────────────────────────────────────────────────────┐  │
│  │ VelocitySource (velocity abstraction)                     │  │
│  │ • StaticFunctionVelocity      - u(x)                     │  │
│  │ • TimeDependentVelocity       - u(x,t)                   │  │
│  │ • FEVelocitySource + NarrowBand - FE velocity (FAST)     │  │
│  └──────────────────────────────────────────────────────────┘  │
└─────────────────────────────────────────────────────────────────┘
                              ↓
┌─────────────────────────────────────────────────────────────────┐
│              External Dependencies                              │
│  • LevelSetMethods.jl (WENO5, RK3, reinitialization)           │
│  • Gridap + GridapEmbedded (FE spaces, cut meshes)             │
└─────────────────────────────────────────────────────────────────┘
```

## The Level Set Function

The `initial_ls` parameter expects a **callable function** `ϕ₀(x)` where:
- `x` is a coordinate tuple `(x[1], x[2])` for 2D
- Returns a `Float64` value representing the signed distance
- **ϕ < 0**: Inside the domain (physical region)
- **ϕ = 0**: Interface (boundary)
- **ϕ > 0**: Outside

## Initial Geometry: Analytical Functions

Define geometry using closed-form expressions. This is the most robust approach.

### Primitives

```julia
# Circle: center (cx, cy), radius R
ϕ_circle(x) = sqrt((x[1] - cx)^2 + (x[2] - cy)^2) - R

# Rectangle: centered at (cx, cy), half-widths (hw, hh)
function ϕ_rectangle(x)
    dx = abs(x[1] - cx) - hw
    dy = abs(x[2] - cy) - hh
    outside = sqrt(max(dx, 0)^2 + max(dy, 0)^2)
    inside = min(max(dx, dy), 0.0)
    return outside + inside
end

# Ellipse: center (cx, cy), semi-axes (a, b)
ϕ_ellipse(x) = sqrt(((x[1]-cx)/a)^2 + ((x[2]-cy)/b)^2) - 1.0

# Half-plane: normal (nx, ny), point (px, py)
ϕ_halfplane(x) = (x[1] - px) * nx + (x[2] - py) * ny
```

### CSG Operations

Combine primitives using min/max operations (Constructive Solid Geometry):

```julia
# Union: A ∪ B (inside either A or B)
ϕ_union(x) = min(ϕ_A(x), ϕ_B(x))

# Intersection: A ∩ B (inside both A and B)
ϕ_intersection(x) = max(ϕ_A(x), ϕ_B(x))

# Difference: A \ B (inside A, outside B)
ϕ_difference(x) = max(ϕ_A(x), -ϕ_B(x))

# Complement: ¬A (outside A)
ϕ_complement(x) = -ϕ_A(x)
```

### Example: Zalesak Disk

A classic test case: a slotted disk with a rectangular notch.

```julia
function zalesak_disk(x; center=(0.0, 0.25), radius=0.5, slot_width=0.12, slot_height=0.5)
    cx, cy = center
    
    # Distance to circle (negative inside)
    d_circle = sqrt((x[1] - cx)^2 + (x[2] - cy)^2) - radius
    
    # Slot rectangle
    xmin, xmax = cx - slot_width/2, cx + slot_width/2
    ymin, ymax = cy - radius, cy - radius + slot_height
    
    # Signed distance to rectangle
    dx = max(xmin - x[1], x[1] - xmax, 0.0)
    dy = max(ymin - x[2], x[2] - ymax, 0.0)
    inside_slot = (xmin ≤ x[1] ≤ xmax) && (ymin ≤ x[2] ≤ ymax)
    d_slot = inside_slot ? -min(x[1]-xmin, xmax-x[1], x[2]-ymin, ymax-x[2]) : sqrt(dx^2 + dy^2)
    
    # CSG: disk \ slot = max(-ϕ_disk, ϕ_slot)
    return max(-d_circle, d_slot)
end
```

## Code Structure Patterns

### Pattern 1: One-Way Coupling (Prescribed Velocity)

Use when velocity is known analytically or from external source.

```julia
# Setup
u(x, t) = (-x[2] + 0.5, x[1] - 0.5) * cos(t)  # Time-varying rotation
evolver = LevelSetMethodsEvolver(; velocity = u, ...)
eg = EvolvingDiscreteGeometry(evolver, model)

# Time loop - geometry doesn't affect velocity
for step in 1:N
    advance!(eg, Δt)
    geo = current_geometry(eg)
    # Optional: solve physics on current geometry
end
```

### Pattern 2: Two-Way Coupling (FE Velocity from Physics)

Use when velocity comes from solving equations on the moving domain.

```julia
# Setup with narrow band (FAST path)
nx, ny = partition .+ 1
Δx = 1.0 / partition[1]
γ = 6 * Δx  # WENO5 needs 6 cells
ext = NarrowBandExtension(γ, nx, ny)

# Initial FE velocity (placeholder)
V = FESpace(model, ReferenceFE(lagrangian, VectorValue{2,Float64}, 1))
u₀ = interpolate_everywhere(x -> VectorValue(0.0, 0.0), V)

vel_source = FEVelocitySource(u₀, model, ext)
evolver = LevelSetMethodsEvolver(; velocity = vel_source, ...)
eg = EvolvingDiscreteGeometry(evolver, model)

# Time loop - velocity depends on geometry
for step in 1:N
    # 1. Get current geometry
    geo = current_geometry(eg)
    cut_geo = cut(model, geo)
    
    # 2. Solve physics (e.g., Stokes)
    velocity_fh = solve_stokes(cut_geo)
    
    # 3. Update velocity source
    update_velocity!(vel_source, velocity_fh)
    update_levelset!(vel_source, current_values(eg.evolver))
    
    # 4. Advance geometry
    advance!(eg, Δt)
end
```

## Key Concepts

### Signed Distance Property

For accurate advection, ϕ should be a **signed distance function**:
- |∇ϕ| = 1 everywhere
- ϕ(x) = distance from x to interface

Use `reinitialize!` periodically to restore this property:

```julia
if step % 10 == 0
    reinitialize!(eg.evolver)
end
```

### CFL Condition

Time step must satisfy:

```
Δt < Δx / max(|u|)
```

For WENO5 + RK3, use CFL ≈ 0.5:

```julia
Δx = 1.0 / partition[1]
u_max = maximum(norm.(velocities))
Δt = 0.5 * Δx / u_max
```

### Narrow Band Concept

Velocity is only needed near the interface (|ϕ| < γ):
- **Inside band**: Velocity from FE solution or extended
- **Outside band**: Velocity = 0 (doesn't affect interface)

Bandwidth γ should cover the stencil:
```julia
γ = 6 * Δx  # WENO5 uses 6-point stencil
```

## Performance Optimization

### Do's ✓

```julia
# ✓ Use NarrowBandExtension for FE velocities (474× faster)
ext = NarrowBandExtension(γ, nx, ny)
vel_source = FEVelocitySource(fh, model, ext)
update_levelset!(vel_source, ϕ_values)

# ✓ Pre-allocate coordinate arrays
coords = [Tuple(c) for c in Gridap.Geometry.get_node_coordinates(trian)]

# ✓ Reuse FE spaces when possible
V = FESpace(model, reffe)  # Create once
for step in 1:N
    uh = interpolate_everywhere(u_func, V)  # Reuse V
end

# ✓ Use the cached cut geometry
cut_geo = current_cut(eg)  # Cached, fast
```

### Don'ts ✗

```julia
# ✗ Don't forget to update level set for narrow band
vel_source = FEVelocitySource(fh, model, ext)
# Missing: update_levelset!(vel_source, ϕ_values)
# Result: Falls back to slow point evaluation

# ✗ Don't recreate FE spaces each step
for step in 1:N
    V = FESpace(model, reffe)  # ✗ Expensive!
end

# ✗ Don't use tiny time steps unnecessarily
Δt = 1e-6  # ✗ Too small if CFL allows larger

# ✗ Don't call reinitialize! every step
for step in 1:N
    reinitialize!(evolver)  # ✗ Expensive, do every 10-20 steps
end
```

## Common Pitfalls

### Level Set Blows Up

**Symptom**: ϕ values become very large or NaN

**Fix**:
```julia
Δt = 0.5 * Δx / u_max  # Respect CFL
if step % 10 == 0
    reinitialize!(eg.evolver)
end
```

### Interface Disappears

**Symptom**: Geometry becomes empty or fills entire domain

**Fix**: Check your sign convention — ϕ < 0 is INSIDE.

### Slow Performance

**Fix**: Use narrow band extension:
```julia
ext = NarrowBandExtension(γ, nx, ny)
vel_source = FEVelocitySource(fh, model, ext)
update_levelset!(vel_source, current_values(evolver))
```
