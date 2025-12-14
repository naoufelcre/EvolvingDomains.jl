# [Colliding Balls](@id colliding_balls)

Two balls moving toward each other and merging.

## Code

```julia
using EvolvingDomains
using Gridap, GridapEmbedded
using LevelSetMethods

# Domain setup
domain = (0.0, 2.0, 0.0, 2.0)
partition = (80, 80)
model = CartesianDiscreteModel(domain, partition)

# Two balls above and below y=1.0
center_top = (1.0, 1.3)
center_bottom = (1.0, 0.7)
radius = 0.2

function ϕ0(x)
    # Positive inside convention: positive = inside ball
    d_top = radius - sqrt((x[1] - center_top[1])^2 + (x[2] - center_top[2])^2)
    d_bottom = radius - sqrt((x[1] - center_bottom[1])^2 + (x[2] - center_bottom[2])^2)
    return max(d_top, d_bottom)  # Union of two balls
end

# Velocity: move toward y=1.0
velocity(x) = x[2] > 1.0 ? (0.0, -0.5) : (0.0, 0.5)

# Create evolver
evolver = LevelSetMethodsEvolver(;
    bg_model = model,
    initial_ls = ϕ0,
    velocity = velocity,
    spatial_scheme = :WENO5,
    integrator = :RK3,
    bc = :Neumann
)
eg = EvolvingDiscreteGeometry(evolver, model)

# Time stepping
Δt = 0.01
substeps = 2
nsteps = 30

for step in 1:nsteps
    for _ in 1:substeps
        EvolvingDomains.advance!(eg, Δt)
    end

    # Reinitialize periodically
    if step % 10 == 0
        EvolvingDomains.reinitialize!(eg)
    end
end

println("Final time: ", EvolvingDomains.current_time(eg))
```

## Visualization

With GLMakie loaded:

```julia
using GLMakie

# Initial state
fig_initial = plot_levelset(eg; title = "Initial: Two Balls")

# Final state (after collision)
fig_final = plot_levelset(eg; title = "Final: Merged")
```

## Results

TODO
