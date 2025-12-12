# EvolvingDomains.jl

[![Build Status](https://github.com/naoufelcre/EvolvingDomains.jl/actions/workflows/CI.yml/badge.svg?branch=main)](https://github.com/naoufelcre/EvolvingDomains.jl/actions/workflows/CI.yml?query=branch%3Amain)

**A Julia package for evolving domain simulations with level set methods and cut-cell finite elements.**

EvolvingDomains.jl provides a clean, composable API for simulating physics on moving domains. It combines high-order level set advection ([LevelSetMethods.jl](https://github.com/maltezfaria/LevelSetMethods.jl)) with embedded finite element methods ([Gridap](https://github.com/gridap/Gridap.jl) + [GridapEmbedded](https://github.com/gridap/GridapEmbedded.jl)).

---

## Features

- 🚀 **High-order advection** — WENO5 spatial schemes with RK3 time stepping
- ⚡ **Fast FE coupling** — Narrow-band velocity extension (100x+ speedup)
- 🔌 **Pluggable backends** — Abstract interface for custom solvers
- 🧩 **External solver integration** — Clean API for operator-splitting workflows

---

## Installation

```julia
using Pkg
Pkg.add(url="https://github.com/naoufelcre/EvolvingDomains.jl")

# Also install the level set backend (not yet registered)
Pkg.add(url="https://github.com/maltezfaria/LevelSetMethods.jl")
```

---

## Quick Start

```julia
using EvolvingDomains
using Gridap, GridapEmbedded
using LevelSetMethods

# 1. Create a background Cartesian grid
model = CartesianDiscreteModel((0.0, 1.0, 0.0, 1.0), (50, 50))

# 2. Define initial geometry (circle) and velocity
ϕ₀(x) = 0.2 - sqrt((x[1]-0.5)^2 + (x[2]-0.5)^2)
u(x) = (0.5, 0.0)  # Rightward flow

# 3. Create the evolver and evolving geometry
evolver = LevelSetMethodsEvolver(;
    bg_model = model,
    initial_ls = ϕ₀,
    velocity = u,
    spatial_scheme = :WENO5
)
eg = EvolvingDiscreteGeometry(evolver, model)

# 4. Evolve!
for step in 1:100
    advance!(eg, 0.01)
    
    # Get current geometry for CutFEM
    cut_geo = current_cut(eg)
    # ... solve physics on cut_geo ...
end
```

📖 **See the [User Guide](docs/USER_GUIDE.md) for detailed examples and patterns.**

---

## API Overview

### Core Types

| Type | Description |
|------|-------------|
| `EvolvingDiscreteGeometry` | Main container — manages level set + Gridap geometry |
| `LevelSetMethodsEvolver` | Backend using LevelSetMethods.jl for advection |
| `FEVelocitySource` | Couples FE velocity solutions to level set evolution |

### Core Functions

| Function | Description |
|----------|-------------|
| `advance!(eg, Δt)` | Evolve the geometry by time Δt |
| `current_geometry(eg)` | Get `DiscreteGeometry` for use with `cut()` |
| `current_cut(eg)` | Get cached `EmbeddedDiscretization` |
| `current_levelset(eg)` | Get level set values as `Vector{Float64}` |
| `reinitialize!(eg)` | Restore signed distance property |

---

## External Solver Integration

EvolvingDomains.jl is designed for **operator-splitting workflows** where you might want to:
- Solve hyperbolic equations on the fixed Cartesian grid (external solver)
- Solve elliptic equations with CutFEM (Gridap)

### Getting Grid Metadata

```julia
info = grid_info(eg)

info.origin   # (0.0, 0.0) — grid origin
info.spacing  # (0.02, 0.02) — cell sizes Δx, Δy
info.dims     # (51, 51) — node counts
info.cells    # (50, 50) — cell counts
```

### Domain Masks

```julia
ϕ = current_levelset(eg)

# Which nodes are inside the domain?
inside = domain_mask(ϕ)  # BitVector: true where ϕ < 0

# Which nodes are near the interface?
γ = 6 * info.spacing[1]  # 6-cell bandwidth
near_interface = narrow_band_mask(ϕ, γ)  # BitVector: true where |ϕ| < γ
```

### Injecting External Solutions

If your external solver evolves the level set:

```julia
# Get current state
ϕ = current_levelset(eg)
info = grid_info(eg)
mask = domain_mask(ϕ)

# External solver computes new level set
ϕ_new = my_external_solver(info, ϕ, mask, Δt)

# Inject back into EvolvingDomains
set_levelset!(eg, ϕ_new)
reinitialize!(eg)  # Optional: restore signed distance

# Continue with CutFEM
cut_geo = current_cut(eg)
```

### Typical Operator-Splitting Loop

```julia
for step in 1:nsteps
    # 1. Hyperbolic step (external FD/FV solver)
    ϕ = current_levelset(eg)
    ϕ_new = hyperbolic_step(grid_info(eg), ϕ, u_data, Δt)
    set_levelset!(eg, ϕ_new)
    reinitialize!(eg)
    
    # 2. Elliptic step (CutFEM via Gridap)
    cut_geo = current_cut(eg)
    u_new = solve_elliptic(cut_geo)
    
    # 3. (Optional) Update velocity for next step
    update_velocity!(eg, FEVelocitySource(u_new, model))
end
```

---

## Velocity Sources

EvolvingDomains supports multiple velocity representations:

```julia
# Static function: u(x)
vel = StaticFunctionVelocity(x -> (1.0, 0.0))

# Time-dependent: u(x, t)
vel = TimeDependentVelocity((x, t) -> (-x[2], x[1]) * cos(t))

# FE-coupled (from Gridap solution)
vel = FEVelocitySource(velocity_fh, model)

# With narrow-band extension (FAST for cut geometries)
ext = NarrowBandExtension(γ, nx, ny)
vel = FEVelocitySource(velocity_fh, model, ext)
```

---

## Performance Tips

✅ **Do:**
- Use `NarrowBandExtension` for FE velocities (100x+ faster)
- Call `current_cut(eg)` instead of `cut(model, current_geometry(eg))` — it's cached
- Reinitialize every 10-20 steps, not every step

❌ **Avoid:**
- Very small time steps beyond CFL requirements
- Recreating FE spaces each time step
- Forgetting to call `update_levelset!` for narrow-band sources

---

## Documentation

- 📖 **[User Guide](docs/USER_GUIDE.md)** — Comprehensive guide with patterns and examples
- 🔬 **[API Reference](docs/USER_GUIDE.md#api-reference)** — Full API documentation

---

## Contributing

Contributions are welcome! Please open an issue to discuss proposed changes.

---

## License

MIT License — see [LICENSE](LICENSE) for details.
