# EvolvingDomains User Guide

A practical guide to implementing moving domain simulations with FE-coupled level set advection.

## Table of Contents

1. [Architecture Overview](#architecture-overview)
2. [Quick Start Example](#quick-start-example)
3. [Code Structure Patterns](#code-structure-patterns)
4. [Key Concepts](#key-concepts)
5. [Performance Optimization](#performance-optimization)
6. [Common Pitfalls](#common-pitfalls)
7. [API Reference](#api-reference)

---

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

---

## Quick Start Example

### Minimal Working Example

```julia
using EvolvingDomains
using Gridap
using GridapEmbedded
using LevelSetMethods

# 1. Define background grid
domain = (0.0, 1.0, 0.0, 1.0)
partition = (50, 50)
model = CartesianDiscreteModel(domain, partition)

# 2. Define initial geometry (circle)
R = 0.2
ϕ₀(x) = R - sqrt((x[1]-0.5)^2 + (x[2]-0.5)^2)

# 3. Define velocity field
u(x) = (1.0, 0.0)  # Constant rightward flow

# 4. Create evolver
evolver = LevelSetMethodsEvolver(;
    bg_model = model,
    initial_ls = ϕ₀,
    velocity = u,
    spatial_scheme = :WENO5,
    bc = :Neumann
)

# 5. Create evolving geometry
eg = EvolvingDiscreteGeometry(evolver, model)

# 6. Time stepping
for step in 1:100
    advance!(eg, 0.01)  # Move geometry
    
    # Use current geometry with Gridap
    geo = current_geometry(eg)
    cut_geo = cut(model, geo)
    # ... solve physics on cut_geo ...
end
```

---

## Initial Geometry Creation

This section covers how to define initial level set geometries, including **Constructive Solid Geometry (CSG)** operations for complex shapes.

### The Level Set Function

The `initial_ls` parameter expects a **callable function** `ϕ₀(x)` where:
- `x` is a coordinate tuple `(x[1], x[2])` for 2D
- Returns a `Float64` value representing the signed distance
- **ϕ < 0**: Inside the domain (physical region)
- **ϕ = 0**: Interface (boundary)
- **ϕ > 0**: Outside

### Approach 1: Analytical Functions (Recommended)

Define geometry using closed-form expressions. This is the most robust approach.

#### Primitives

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

#### CSG Operations

Combine primitives using min/max operations:

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

#### Example: Zalesak Disk

A classic test case: a slotted disk with a rectangular notch.

```julia
function zalesak_disk(x; center=(0.0, 0.25), radius=0.5, slot_width=0.12, slot_height=0.5)
    cx, cy = center
    
    # Distance to circle (negative inside)
    d_circle = sqrt((x[1] - cx)^2 + (x[2] - cy)^2) - radius
    
    # Slot rectangle: spans [cx - w/2, cx + w/2] × [cy - r, cy - r + h]
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

# Usage
evolver = LevelSetMethodsEvolver(;
    bg_model = model,
    initial_ls = zalesak_disk,
    velocity = u,
    integrator = :RK3
)
```

#### Example: Multi-Body Union

```julia
function two_circles(x)
    # Circle 1: center (0.3, 0.5), radius 0.15
    ϕ1 = sqrt((x[1] - 0.3)^2 + (x[2] - 0.5)^2) - 0.15
    
    # Circle 2: center (0.7, 0.5), radius 0.15
    ϕ2 = sqrt((x[1] - 0.7)^2 + (x[2] - 0.5)^2) - 0.15
    
    # Union
    return min(ϕ1, ϕ2)
end
```

### Approach 2: Using LevelSetMethods.jl Primitives

If you want to use `LevelSetMethods.jl`'s built-in geometry functions, you must create a compatible `CartesianGrid` first:

```julia
using LevelSetMethods

# Extract grid parameters from Gridap model
(origin, corner, partition) = EvolvingDomains.cartesian_descriptor(model)
node_partition = partition .+ 1  # LevelSetMethods uses node counts

# Create LevelSetMethods grid
lsm_grid = LevelSetMethods.CartesianGrid(origin, corner, node_partition)

# Now use LevelSetMethods primitives
disk = LevelSetMethods.circle(lsm_grid; center=(0.5, 0.5), radius=0.3)
rect = LevelSetMethods.rectangle(lsm_grid; center=(0.5, 0.3), width=(0.1, 0.3))

# CSG operations
shape = LevelSetMethods.setdiff(disk, rect)

# Convert back to function for EvolvingDomains
ϕ_vals = LevelSetMethods.values(shape)

# Use the LevelSetMethods object directly or define analytical equivalent
# (Analytical is recommended for signed distance accuracy)
```

> **⚠️ Important**: `LevelSetMethods.circle/rectangle` expect a `LevelSetMethods.CartesianGrid`, NOT a `Gridap.CartesianDiscreteModel`. Always convert first!

### Best Practices

1. **Use analytical functions** when possible — they maintain exact signed distance property
2. **Normalize your signed distance** — ensure |∇ϕ| ≈ 1 using `reinitialize!` periodically
3. **Keep CSG trees simple** — deep nesting can cause numerical artifacts at corners
4. **Test your geometry** — visualize before running simulations:

```julia
using EmbeddedViz
grid = GridInfo(model)
show_levelset(ϕ_values, grid)
```

---

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

### Pattern 3: Modular Physics Solver

Separate your physics into clean functions:

```julia
# physics.jl
function solve_stokes(cut_geo, model, μ, f)
    # Define FE spaces on cut geometry
    Ω = Triangulation(cut_geo, PHYSICAL_IN)
    Γ = EmbeddedBoundary(cut_geo)
    
    # Weak form, assemble, solve...
    return uh, ph
end

# main.jl
for step in 1:N
    geo = current_geometry(eg)
    cut_geo = cut(model, geo)
    
    uh, ph = solve_stokes(cut_geo, model, μ, f)
    
    update_velocity!(vel_source, uh)
    update_levelset!(vel_source, current_values(eg.evolver))
    advance!(eg, Δt)
end
```

---

## Key Concepts

### 1. Level Set Representation

The geometry is implicitly defined by a level set function ϕ:
- **ϕ < 0**: Inside the domain (physical region)
- **ϕ = 0**: Interface (boundary)
- **ϕ > 0**: Outside

```julia
# Circle of radius R centered at (cx, cy)
ϕ(x) = R - sqrt((x[1]-cx)^2 + (x[2]-cy)^2)

# Rectangle [a,b] × [c,d]
ϕ(x) = max(a-x[1], x[1]-b, c-x[2], x[2]-d)
```

### 2. Signed Distance Property

For accurate advection, ϕ should be a **signed distance function**:
- |∇ϕ| = 1 everywhere
- ϕ(x) = distance from x to interface

Use `reinitialize!` periodically to restore this property:

```julia
if step % 10 == 0
    reinitialize!(eg.evolver)
end
```

### 3. CFL Condition

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

### 4. Narrow Band Concept

Velocity is only needed near the interface (|ϕ| < γ):
- **Inside band**: Velocity from FE solution or extended
- **Outside band**: Velocity = 0 (doesn't affect interface)

Bandwidth γ should cover the stencil:
```julia
γ = 6 * Δx  # WENO5 uses 6-point stencil
```

### 5. DOF-to-Node Mapping

For Lagrangian P1 vector fields on Cartesian grids:
- Nodes numbered in row-major order (x varies fastest)
- DOFs are interleaved: `[u1_x, u1_y, u2_x, u2_y, ...]`
- DOF `2i-1` and `2i` correspond to node `i`

---

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

# ✓ Use WENO5 for accuracy, Upwind for speed
spatial_scheme = :WENO5  # More accurate
spatial_scheme = :Upwind # Faster, less accurate
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

---

## Common Pitfalls

### 1. Level Set Blows Up

**Symptom**: ϕ values become very large or NaN

**Causes**:
- Time step too large (CFL violation)
- Never reinitializing

**Fix**:
```julia
Δt = 0.5 * Δx / u_max  # Respect CFL
if step % 10 == 0
    reinitialize!(eg.evolver)
end
```

### 2. Interface Disappears

**Symptom**: Geometry becomes empty or fills entire domain

**Causes**:
- Velocity pointing wrong direction
- Sign convention mismatch

**Fix**: Check your sign convention:
```julia
# ϕ < 0 is INSIDE, so interface velocity should point OUTWARD
# for expansion, INWARD for contraction
```

### 3. Slow Performance

**Symptom**: Each time step takes seconds

**Causes**:
- Using point evaluation instead of DOF path
- Narrow band not configured

**Fix**:
```julia
ext = NarrowBandExtension(γ, nx, ny)
vel_source = FEVelocitySource(fh, model, ext)
update_levelset!(vel_source, current_values(evolver))
```

### 4. Wrong Velocity at Interface

**Symptom**: Interface moves in unexpected direction

**Causes**:
- FE velocity only defined inside, but level set needs it at interface
- Extension not working

**Fix**: Check that narrow band covers interface:
```julia
γ = 6 * Δx  # Must be > 0 to include interface nodes
```

---

## API Reference

### Core Types

| Type | Description |
|------|-------------|
| `EvolvingDiscreteGeometry` | Main container for evolving geometry |
| `LevelSetMethodsEvolver` | Level set advection backend |
| `FEVelocitySource` | Wraps Gridap FEFunction for velocity |
| `NarrowBandExtension` | Fast velocity extension strategy |

### Key Functions

| Function | Description |
|----------|-------------|
| `advance!(eg, Δt)` | Move geometry forward by Δt |
| `current_geometry(eg)` | Get DiscreteGeometry for Gridap |
| `current_values(evolver)` | Get level set values at nodes |
| `current_time(eg)` | Get current simulation time |
| `reinitialize!(evolver)` | Restore signed distance property |
| `update_velocity!(source, fh)` | Update FE velocity function |
| `update_levelset!(source, ϕ)` | Update level set for narrow band |

### Constructor Options

```julia
LevelSetMethodsEvolver(;
    bg_model,           # CartesianDiscreteModel
    initial_ls,         # ϕ₀(x) function
    velocity,           # Function or VelocitySource
    spatial_scheme = :WENO5,  # :WENO5 or :Upwind
    integrator = :RK3,        # :RK3, :RK2, :ForwardEuler
    bc = :Neumann,            # :Neumann, :Periodic, :Dirichlet
    reinit_freq = nothing     # Auto-reinitialize frequency
)
```

---

## Next Steps

1. **Start simple**: Use prescribed velocity to verify geometry evolution
2. **Add physics**: Solve equations on cut geometry
3. **Close the loop**: Use physics solution as velocity
4. **Optimize**: Enable narrow band extension for FE velocities
5. **Benchmark**: Compare with analytical solutions if available
