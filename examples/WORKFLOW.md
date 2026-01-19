# EvolvingDomains.jl — Typical Workflow

This guide outlines the standard workflow for setting up and running a simulation with `EvolvingDomains.jl`. The workflow consists of five main steps:

1.  **Grid Setup**
2.  **Geometry Design & Initialization**
3.  **Velocity Definition**
4.  **Simulation Loop**
5.  **Visualization**

---

## 1. Grid Setup

The simulation domain is defined using a `CartesianDiscreteModel` from `Gridap`. This serves as the background mesh for the level set evolution.

```julia
using EvolvingDomains
using Gridap

# Define the physical domain limits (xmin, xmax, ymin, ymax)
domain = (-2.0, 2.0, -2.0, 2.0)

# Define the resolution (number of cells in x and y)
n_cells = (80, 80)

# Create the background model
model = CartesianDiscreteModel(domain, n_cells)
```

> [!TIP]
> Ensure your grid resolution is sufficient to capture the features of your geometry, especially if thin structures (like the slot in Zalesak's disk) are involved.

---

## 2. Geometry Design & Initialization

You can create complex initial geometries using Constructive Solid Geometry (CSG) primitives.

### Primitives
Available primitives in `EvolvingDomains`:
- `Circle(center, radius)`
- `Rectangle(pmin, pmax)`

### Operations
Combine primitives using standard set operations:
- `union(a, b)`: Combine shapes.
- `intersect(a, b)`: Overlapping region.
- `setdiff(a, b)`: Subtract shape `b` from `a`.

### Transformations
- `Translate(geometry, displacement_vector)`

### Example: Constructing a Shape
```julia
# 1. Create a base circle centered at origin
base_shape = Circle(VectorValue(0.0, 0.0), 0.8)

# 2. Create a "cutout" rectangle
cutout = Rectangle(VectorValue(-0.2, 0.0), VectorValue(0.2, 1.0))

# 3. Subtract the cutout from the circle
initial_shape = setdiff(base_shape, cutout)

# 4. Initialize the EvolvingDiscreteGeometry
# Convert the CSG object to a signed distance function (x -> Float64)
sdf(x) = initial_shape(VectorValue(x...))

geom = EvolvingDiscreteGeometry(model, sdf; reinit_freq=5)
```

> [!NOTE]
> `reinit_freq` controls how often the level set field is re-distanced. Frequent reinitialization preserves visual shape but can introduce small mass conservation errors.

---

## 3. Velocity Definition

The geometry evolves according to a velocity field. You can define this analytically or using Finite Elements.

### Option A: Analytic Velocity
Use `TimeDependentVelocity` or `StaticFunctionVelocity`.

```julia
# Rigid body rotation: v = (-y, x)
velocity_field = TimeDependentVelocity((x, t) -> VectorValue(-x[2], x[1]))

set_velocity!(geom, velocity_field)
```

### Option B: FE Velocity (FSI Coupling)
Coupling with a PDE solver often involves an `FEFunction`.

```julia
# Assuming u_h is a Gridap FEFunction from a Stokes/Navier-Stokes solve:
vel_source = FEVelocitySource(u_h)

# Update the source and assign to geometry
update_velocity!(vel_source, u_h, geom) 
set_velocity!(geom, vel_source)
```

---

## 4. Simulation Loop

The core simulation evolves the geometry in time steps `dt`.

```julia
t_end = 2π
dt = 0.05
n_steps = round(Int, t_end / dt)

for step in 1:n_steps
    # 1. Advance the interface
    advance!(geom, dt)
    
    # 2. Access the current state for physics solvers (if needed)
    # Get the current CutFEM discretization
    cut_mesh = current_cut(geom) 
    
    # ... Solve PDEs on cut_mesh ...
    
    # 3. (Optional) Visualization code goes here
end
```

---

## 5. Visualization

Use `CairoMakie` to visualize the evolving interface.

```julia
using CairoMakie

fig = Figure()
ax = Axis(fig[1, 1], aspect=1)

# Plot the zero-level set (the interface)
plot_levelset!(ax, geom; color=:black, linewidth=2)
```

### Animation Loop Example
```julia
record(fig, "simulation.gif", 1:n_steps; framerate=20) do step
    advance!(geom, dt)
    
    empty!(ax) # Clear previous frame
    plot_levelset!(ax, geom)
    ax.title = "Time: $(round(EvolvingDomains.current_time(geom), digits=2))"
end
```
