# [Rotating Circle](@id rotating_circle)

A simple example advecting a circle with a rigid body rotation velocity field.

## Problem Setup

We set up a circular domain on a background Cartesian grid and evolve it under rigid body rotation about the origin.

- **Domain**: ``[-1, 1]^2``
- **Grid**: 100 × 100 cells
- **Initial shape**: Circle of radius 0.2 centered at (0.3, 0.5)
- **Velocity**: Rigid body rotation ``\mathbf{u}(x,y) = \omega(-y, x)`` with ``\omega = 2\pi``

## Code

```julia
using EvolvingDomains
using Gridap, GridapEmbedded
using LevelSetMethods

# Domain setup
domain = (-1.0, 1.0, -1.0, 1.0)
partition = (100, 100)
model = CartesianDiscreteModel(domain, partition)

# Initial circle
R = 0.2
center = (0.3, 0.5)
ϕ₀(x) = R - sqrt((x[1]-center[1])^2 + (x[2]-center[2])^2)

# Rigid body rotation about origin
ω = 2π  # Angular velocity (one revolution in T=1)
u(x) = (-ω * x[2], ω * x[1])

# Create evolver
evolver = LevelSetMethodsEvolver(;
    bg_model = model,
    initial_ls = ϕ₀,
    velocity = u,
    spatial_scheme = :WENO5,
    bc = :Neumann
)
eg = EvolvingDiscreteGeometry(evolver, model)

# Time stepping
Δt = 0.01
substeps = 2
nsteps = 50

for step in 1:nsteps
    for _ in 1:substeps
        EvolvingDomains.advance!(eg, Δt)
    end

    # Reinitialize periodically to maintain signed distance property
    if step % 5 == 0
        EvolvingDomains.reinitialize!(eg)
    end
end

println("Final time: ", EvolvingDomains.current_time(eg))
```

## Visualization

With GLMakie loaded, you can visualize the final state:

```julia
using GLMakie
fig = plot_levelset(eg; title = "Rotating Circle (Final)")
```

## Results

The circle rotates about the origin while maintaining its shape. After one complete revolution (T=1), the circle should return to its starting position.

![Rotating Circle Animation](../assets/circle_evolution_levelset.gif)
