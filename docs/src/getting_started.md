# Getting Started

## Installation

```julia
using Pkg
Pkg.add(url="https://github.com/naoufelcre/EvolvingDomains.jl")

# Also install the level set backend (not yet registered)
Pkg.add(url="https://github.com/maltezfaria/LevelSetMethods.jl")
```

## Prerequisites

EvolvingDomains.jl builds on top of several Julia packages:

- **[Gridap.jl](https://github.com/gridap/Gridap.jl)**: Finite element framework
- **[GridapEmbedded.jl](https://github.com/gridap/GridapEmbedded.jl)**: Cut-cell methods for embedded boundaries
- **[LevelSetMethods.jl](https://github.com/maltezfaria/LevelSetMethods.jl)**: High-order level set advection

These will be installed automatically as dependencies.

## Basic Workflow

### 1. Create a Background Grid

```julia
using Gridap

domain = (0.0, 1.0, 0.0, 1.0)
partition = (50, 50)
model = CartesianDiscreteModel(domain, partition)
```

### 2. Define Initial Geometry

The geometry is defined implicitly via a **level set function** ``\phi(x)``:
- ``\phi < 0``: Inside the domain
- ``\phi = 0``: Interface
- ``\phi > 0``: Outside

```julia
# Circle centered at (0.5, 0.5) with radius 0.2
R = 0.2
ϕ₀(x) = R - sqrt((x[1]-0.5)^2 + (x[2]-0.5)^2)
```

### 3. Create the Evolver and Evolving Geometry

```julia
using EvolvingDomains
using LevelSetMethods

# Define velocity field
u(x) = (1.0, 0.0)  # Constant rightward flow

# Create evolver
evolver = LevelSetMethodsEvolver(;
    bg_model = model,
    initial_ls = ϕ₀,
    velocity = u,
    spatial_scheme = :WENO5,
    bc = :Neumann
)

# Create evolving geometry wrapper
eg = EvolvingDiscreteGeometry(evolver, model)
```

### 4. Time Stepping

```julia
using GridapEmbedded

for step in 1:100
    advance!(eg, 0.01)  # Move geometry

    # Use current geometry with Gridap
    geo = current_geometry(eg)
    cut_geo = cut(model, geo)

    # ... solve physics on cut_geo ...
end
```

## Next Steps

- [Examples](examples/index.md)
- [Visualization](visualization.md)
- [API Reference](api/core.md)

