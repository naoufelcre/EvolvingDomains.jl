# Examples

This page contains complete, runnable examples demonstrating common use cases.

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

## Example 3: FE-Coupled Velocity

Two-way coupling where velocity comes from solving an FE problem.

```julia
using EvolvingDomains
using Gridap, GridapEmbedded
using LevelSetMethods

# Domain setup
domain = (0.0, 1.0, 0.0, 1.0)
partition = (50, 50)
model = CartesianDiscreteModel(domain, partition)

# Initial geometry
ϕ₀(x) = 0.2 - sqrt((x[1]-0.5)^2 + (x[2]-0.5)^2)

# Setup narrow band extension (crucial for performance)
nx, ny = partition .+ 1
Δx = 1.0 / partition[1]
γ = 6 * Δx  # WENO5 stencil width
ext = NarrowBandExtension(γ, nx, ny)

# Initial velocity field (placeholder)
V = FESpace(model, ReferenceFE(lagrangian, VectorValue{2,Float64}, 1))
u_init = interpolate_everywhere(x -> VectorValue(0.0, 0.0), V)

# Create velocity source with extension
vel_source = FEVelocitySource(u_init, model, ext)

# Create evolver
evolver = LevelSetMethodsEvolver(;
    bg_model = model,
    initial_ls = ϕ₀,
    velocity = vel_source,
    spatial_scheme = :WENO5,
    bc = :Neumann
)
eg = EvolvingDiscreteGeometry(evolver, model)

# Simulation loop
Δt = 0.01
for step in 1:100
    # 1. Get current geometry
    cut_geo = current_cut(eg)
    
    # 2. Your physics solver would go here
    #    velocity_fh = solve_stokes(cut_geo, ...)
    
    # For demo: use a simple prescribed velocity
    velocity_fh = interpolate_everywhere(x -> VectorValue(0.1, 0.0), V)
    
    # 3. Update velocity source (IMPORTANT!)
    update_velocity!(vel_source, velocity_fh)
    update_levelset!(vel_source, current_levelset(eg))
    
    # 4. Advance geometry
    advance!(eg, Δt)
    
    # 5. Periodic reinitialization
    if step % 10 == 0
        reinitialize!(eg)
    end
end
```

## Example 4: External Solver Integration

Using grid info for operator-splitting with external Cartesian-grid solvers.

```julia
using EvolvingDomains
using Gridap, GridapEmbedded
using LevelSetMethods

# Setup
domain = (0.0, 1.0, 0.0, 1.0)
partition = (50, 50)
model = CartesianDiscreteModel(domain, partition)

ϕ₀(x) = 0.2 - sqrt((x[1]-0.5)^2 + (x[2]-0.5)^2)
u(x) = (0.5, 0.0)

evolver = LevelSetMethodsEvolver(;
    bg_model = model,
    initial_ls = ϕ₀,
    velocity = u,
    spatial_scheme = :WENO5
)
eg = EvolvingDiscreteGeometry(evolver, model)

# Get grid metadata for external solver
info = grid_info(eg)
println("Grid origin: ", info.origin)
println("Grid spacing: ", info.spacing)
println("Node dimensions: ", info.dims)
println("Cell dimensions: ", info.cells)

# Example: extract level set data for external processing
ϕ = current_levelset(eg)
inside_mask = domain_mask(ϕ)
narrow_band = narrow_band_mask(ϕ, 6 * info.spacing[1])

println("Total nodes: ", length(ϕ))
println("Nodes inside domain: ", sum(inside_mask))
println("Nodes in narrow band: ", sum(narrow_band))

# Simulate external solver modifying the level set
function my_external_solver(ϕ, info, dt)
    # Example: simple shift (not a real solver!)
    return ϕ .- 0.01
end

# Inject external solution
ϕ_new = my_external_solver(ϕ, info, 0.01)
set_levelset!(eg, ϕ_new)
reinitialize!(eg)

# Continue with CutFEM
cut_geo = current_cut(eg)
```
