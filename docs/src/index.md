# EvolvingDomains.jl

*A Julia package for evolving domain simulations with level set methods and cut-cell finite elements.*

[![Build Status](https://github.com/naoufelcre/EvolvingDomains.jl/actions/workflows/CI.yml/badge.svg?branch=main)](https://github.com/naoufelcre/EvolvingDomains.jl/actions/workflows/CI.yml?query=branch%3Amain)

## Overview

EvolvingDomains.jl provides a clean, composable API for simulating physics on moving domains. It combines high-order level set advection ([LevelSetMethods.jl](https://github.com/maltezfaria/LevelSetMethods.jl)) with embedded finite element methods ([Gridap](https://github.com/gridap/Gridap.jl) + [GridapEmbedded](https://github.com/gridap/GridapEmbedded.jl)).

## Features

- 🚀 **High-order advection** — WENO5 spatial schemes with RK3 time stepping
- ⚡ **Fast FE coupling** — Narrow-band velocity extension (100x+ speedup)
- 🔌 **Pluggable backends** — Abstract interface for custom solvers
- 🧩 **External solver integration** — Clean API for operator-splitting workflows

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

## Package Contents

```@index
```

## Manual Outline

```@contents
Pages = [
    "getting_started.md",
    "user_guide.md",
    "examples.md",
    "api/core.md",
    "api/velocity.md",
    "api/external.md",
]
Depth = 2
```
