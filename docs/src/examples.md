# Examples of moving domain construction

Simple construction of moving domains

## Example 1: Advecting a Circle

A simple example advecting a circle with constant velocity.

```julia
using EvolvingDomains
using Gridap, GridapEmbedded
using LevelSetMethods

# Domain setup
domain = (0.0, 2.0, 0.0, 1.0)
partition = (100, 50)
model = CartesianDiscreteModel(domain, partition)

# Initial circle at left side
R = 0.2
center = (0.3, 0.5)
ϕ₀(x) = R - sqrt((x[1]-center[1])^2 + (x[2]-center[2])^2)

# Constant rightward velocity
u(x) = (1.0, 0.0)

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
T_final = 1.0
nsteps = Int(T_final / Δt)

for step in 1:nsteps
    EvolvingDomains.advance!(eg, Δt)

    # Reinitialize periodically
    if step % 10 == 0
        EvolvingDomains.reinitialize!(eg)
    end

    # Access geometry if needed
    cut_geo = EvolvingDomains.current_cut(eg)
end
```


## Example 2: Rotating Zalesak Disk

```julia
using EvolvingDomains
using Gridap, GridapEmbedded
using LevelSetMethods

# Domain: [-1.5, 1.5]²
domain = (-1.5, 1.5, -1.5, 1.5)
partition = (100, 100)
model = CartesianDiscreteModel(domain, partition)

# Zalesak disk: circle with rectangular slot carved out
function zalesak_disk(x)
    center = (-0.75, 0.0)
    radius = 0.5

    # Circle SDF: positive inside, negative outside
    d_circle = radius - sqrt((x[1] - center[1])^2 + (x[2] - center[2])^2)

    # Rectangular slot: centered on circle, extends upward
    h = 1.0   # slot height
    w = 0.2   # slot width

    xmin = center[1] - w/2
    xmax = center[1] + w/2
    ymin = center[2]
    ymax = center[2] + h

    # Slot SDF: positive inside slot, negative outside
    dx = min(x[1] - xmin, xmax - x[1])
    dy = min(x[2] - ymin, ymax - x[2])
    d_slot = min(dx, dy)

    # Zalesak disk = circle - slot (positive inside convention)
    return min(d_circle, -d_slot)
end

# Rigid body rotation about origin
ω = 2π  # One revolution in T=1
u(x) = (-ω * x[2], ω * x[1])

# Create evolver
evolver = LevelSetMethodsEvolver(;
    bg_model = model,
    initial_ls = zalesak_disk,
    velocity = u,
    spatial_scheme = :WENO5,
    integrator = :RK3,
    bc = :Neumann
)
eg = EvolvingDiscreteGeometry(evolver, model)

# Time stepping: one full revolution
Δt = 0.005
substeps = 4
nsteps = 50

for step in 1:nsteps
    for _ in 1:substeps
        EvolvingDomains.advance!(eg, Δt)
    end

    if step % 5 == 0
        EvolvingDomains.reinitialize!(eg)
    end
end

println("Completed one revolution at t = ", EvolvingDomains.current_time(eg))
```

Results

![Contour](assets/evolution_levelset.gif)
![Grid](assets/evolution_grid.gif)
![Cut Mesh](assets/evolution_cutmesh.gif)


## Example 3: Colliding Balls

topological changes.

```julia
using EvolvingDomains
using Gridap, GridapEmbedded
using LevelSetMethods

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
    return max(d_top, d_bottom)  # Union of two balls (positive inside)
end

# Velocity: move toward y=1.0
velocity(x) = x[2] > 1.0 ? (0.0, -0.5) : (0.0, 0.5)

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
nsteps = 60

for step in 1:nsteps
    EvolvingDomains.advance!(eg, Δt)
    if step % 10 == 0
        EvolvingDomains.reinitialize!(eg)
    end
end
```
