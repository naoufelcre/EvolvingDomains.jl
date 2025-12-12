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
    advance!(eg, Δt)

    # Reinitialize periodically
    if step % 10 == 0
        reinitialize!(eg)
    end

    # Access geometry if needed
    cut_geo = current_cut(eg)
    # ... solve physics here ...
end

println("Final time: ", current_time(eg))
```

Results

![Contour](assets/circle_evolution_levelset.gif)
![Grid](assets/circle_evolution_grid.gif)
![Cut Mesh](assets/circle_evolution_cutmesh.gif)


## Example 2: Rotating Zalesak Disk

The classic benchmark for level set methods: a slotted disk under rigid body rotation.

```julia
using EvolvingDomains
using Gridap, GridapEmbedded
using LevelSetMethods

# Domain: [-0.5, 0.5]²
domain = (-0.5, 0.5, -0.5, 0.5)
partition = (100, 100)
model = CartesianDiscreteModel(domain, partition)

# Zalesak disk
function zalesak_disk(x)
    cx, cy = 0.0, 0.25
    radius = 0.15
    slot_width = 0.05
    slot_height = 0.25

    d_circle = sqrt((x[1] - cx)^2 + (x[2] - cy)^2) - radius

    xmin, xmax = cx - slot_width/2, cx + slot_width/2
    ymin, ymax = cy - radius, cy - radius + slot_height

    dx = max(xmin - x[1], x[1] - xmax, 0.0)
    dy = max(ymin - x[2], x[2] - ymax, 0.0)
    inside_slot = (xmin ≤ x[1] ≤ xmax) && (ymin ≤ x[2] ≤ ymax)
    d_slot = inside_slot ? -min(x[1]-xmin, xmax-x[1], x[2]-ymin, ymax-x[2]) : sqrt(dx^2 + dy^2)

    return max(-d_circle, d_slot)
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
Δt = 0.01
T_period = 1.0
nsteps = Int(T_period / Δt)

for step in 1:nsteps
    advance!(eg, Δt)

    if step % 10 == 0
        reinitialize!(eg)
    end
end

println("Completed one revolution at t = ", current_time(eg))
```

Results

![Contour](assets/evolution_levelset.gif)
![Grid](assets/evolution_grid.gif)
![Cut Mesh](assets/evolution_cutmesh.gif)
