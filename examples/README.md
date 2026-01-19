# EvolvingDomains.jl — Examples

This directory contains working examples demonstrating the core functionality of `EvolvingDomains.jl`.

## Quick Start

Run any example from the project root:

```bash
julia --project=. examples/zalesak_minimal.jl
```

Or from the Julia REPL:

```julia
include("examples/zalesak_minimal.jl")
run_zalesak()
```

---

## Working Examples

### 1. `zalesak_minimal.jl` — Zalesak Disk Rotation

**Description**: Classic benchmark for level set advection. A slotted disk (Zalesak disk) rotates 360° around the origin.

**Demonstrates**:
- CSG operations (`setdiff` for creating the slot)
- WENO5 level set advection with reinitialization
- Rigid body rotation velocity field

**Physics**: One full revolution at constant angular velocity using `v(x) = (-y, x)`.

**Output**: `zalesak_minimal.gif`

---

### 2. `colliding_balls_minimal.jl` — Colliding Balls

**Description**: Two balls move toward each other and merge mid-simulation.

**Demonstrates**:
- Topological changes (domain merging)
- Handling of shock/rarefaction waves in level set methods
- Discontinuous velocity fields (different velocities above/below midline)

**Physics**: Top ball moves down, bottom ball moves up with constant velocities.

**Output**: `colliding_balls_minimal.gif`

---

### 3. `independent_motion_minimal.jl` — Independent Object Motion

**Description**: Two objects (Pacman and a dot) move independently through each other without merging.

**Demonstrates**:
- **Multi-level-set pattern**: Each object tracked by its own `EvolvingDiscreteGeometry`
- Independent velocity fields per object
- Manual level set union for visualization

**Key Pattern**:
```julia
# Create separate geometries for independent motion
geom_pacman = EvolvingDiscreteGeometry(model, φ_pacman)
geom_dot = EvolvingDiscreteGeometry(model, φ_dot)

# Independent velocities
set_velocity!(geom_pacman, StaticFunctionVelocity(x -> (0.5, 0.0)))
set_velocity!(geom_dot, StaticFunctionVelocity(x -> (-0.5, 0.0)))

# Advance separately
advance!(geom_pacman, dt)
advance!(geom_dot, dt)

# Combine for viz
φ_combined = min.(current_levelset(geom_pacman), current_levelset(geom_dot))
```

**Output**: `independent_motion_minimal.gif`

---

### 4. `csg.jl` — Constructive Solid Geometry

**Description**: Static visualization of CSG operations.

**Demonstrates**:
- Primitives: `Circle`, `Rectangle`
- Operations: `union`, `intersect`, `setdiff`
- Static geometry (no advection)

**Note**: Uses `GLMakie` (requires display). Outputs to `EvolvingDomainsDocumentation/assets/`.

---

## API Summary

| Function | Purpose |
|----------|---------|
| `EvolvingDiscreteGeometry(model, φ)` | Create evolving geometry from a Gridap model and level set function |
| `set_velocity!(geom, vel)` | Assign velocity source |
| `advance!(geom, dt)` | Evolve geometry by one time step |
| `current_levelset(geom)` | Get current level set array |
| `plot_levelset!(ax, geom)` | Visualize geometry (requires CairoMakie) |

### Velocity Sources

- `StaticFunctionVelocity(x -> ...)` — Time-independent velocity
- `TimeDependentVelocity((x, t) -> ...)` — Time-dependent velocity
- `FEVelocitySource(fespace, fefunc)` — Finite element velocity (for coupled simulations)
