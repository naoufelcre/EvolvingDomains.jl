# EvolvingDomains.jl Examples

This directory contains runnable example scripts demonstrating various features of EvolvingDomains.jl.

## Running Examples

From the package root directory:

```bash
# Activate the package environment
julia --project=.

# Run an example
julia --project=. examples/basic_advection.jl
```

Or from the Julia REPL:

```julia
# First activate the environment
using Pkg
Pkg.activate(".")

# Then include the example
include("examples/basic_advection.jl")
```

## Available Examples

### 1. Basic Advection (`basic_advection.jl`)

**Difficulty**: Beginner

Demonstrates the core workflow:
- Background Cartesian grid setup
- Level set definition (circle)
- Constant velocity field
- Time stepping with reinitialization
- Accessing geometry for CutFEM

This is the recommended starting point for new users.

### 2. Zalesak's Rotating Disk (`zalesak_disk.jl`)

**Difficulty**: Intermediate

Classic benchmark for level set methods:
- Non-convex slotted disk geometry
- Rigid body rotation velocity
- Error analysis after one revolution
- Tests corner and thin feature preservation

This example validates the numerical accuracy of the scheme.

## Creating Your Own Examples

Use this template as a starting point:

```julia
using EvolvingDomains
using Gridap, GridapEmbedded
using LevelSetMethods

# 1. Grid setup
domain = (x_min, x_max, y_min, y_max)
partition = (nx, ny)
model = CartesianDiscreteModel(domain, partition)

# 2. Initial level set (ϕ < 0 inside, ϕ > 0 outside)
ϕ₀(x) = your_signed_distance_function(x)

# 3. Velocity field
velocity(x) = (vx, vy)

# 4. Create evolver
evolver = LevelSetMethodsEvolver(;
    bg_model = model,
    initial_ls = ϕ₀,
    velocity = velocity,
    spatial_scheme = :WENO5,
    bc = :Neumann
)
eg = EvolvingDiscreteGeometry(evolver, model)

# 5. Time loop
for step in 1:nsteps
    advance!(eg, Δt)
    if step % reinit_freq == 0
        reinitialize!(eg)
    end
    
    # Access geometry for your physics solver
    cut_geo = current_cut(eg)
end
```

## Dependencies

These examples require:
- `EvolvingDomains.jl` (this package)
- `Gridap.jl`
- `GridapEmbedded.jl`
- `LevelSetMethods.jl`

All dependencies are listed in the main `Project.toml`.
