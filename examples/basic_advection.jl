# Basic Advection Example
# =======================
# 
# This example demonstrates the core workflow of EvolvingDomains.jl:
# 1. Setting up a background Cartesian grid
# 2. Defining an initial level set (circle)
# 3. Specifying a velocity field
# 4. Evolving the geometry over time
#
# To run this example:
#   julia --project=. examples/basic_advection.jl
#
# Or from the Julia REPL:
#   include("examples/basic_advection.jl")

using EvolvingDomains
using Gridap
using GridapEmbedded
using LevelSetMethods

# =============================================================================
# 1. Background Grid Setup
# =============================================================================
# Define the computational domain and discretization
domain = (0.0, 2.0, 0.0, 1.0)   # (x_min, x_max, y_min, y_max)
partition = (100, 50)           # Number of cells in each direction
model = CartesianDiscreteModel(domain, partition)

println("✓ Created background grid: $(partition[1])×$(partition[2]) cells")
println("  Domain: x ∈ [$(domain[1]), $(domain[2])], y ∈ [$(domain[3]), $(domain[4])]")

# =============================================================================
# 2. Initial Geometry (Level Set)
# =============================================================================
# Define a circle using a signed distance function
# ϕ(x) < 0 : inside the circle
# ϕ(x) > 0 : outside the circle
# ϕ(x) = 0 : on the boundary

R = 0.2                      # Circle radius
center = (0.3, 0.5)          # Circle center

# Signed distance function: negative inside, positive outside
ϕ₀(x) = sqrt((x[1] - center[1])^2 + (x[2] - center[2])^2) - R

println("\n✓ Initial geometry: Circle")
println("  Center: $center")
println("  Radius: $R")

# =============================================================================
# 3. Velocity Field
# =============================================================================
# Simple constant rightward velocity
velocity(x) = (0.5, 0.0)  # (u_x, u_y)

println("\n✓ Velocity field: Constant rightward")
println("  v = (0.5, 0.0)")

# =============================================================================
# 4. Create the Evolver
# =============================================================================
# The LevelSetMethodsEvolver handles the numerical scheme for level set advection
evolver = LevelSetMethodsEvolver(;
    bg_model = model,
    initial_ls = ϕ₀,
    velocity = velocity,
    spatial_scheme = :WENO5,    # High-order WENO scheme
    integrator = :RK3,          # Third-order Runge-Kutta
    bc = :Neumann               # Boundary conditions
)

# Wrap in EvolvingDiscreteGeometry for CutFEM integration
eg = EvolvingDiscreteGeometry(evolver, model)

println("\n✓ Evolver created")
println("  Spatial scheme: WENO5")
println("  Time integrator: RK3")

# =============================================================================
# 5. Time Evolution
# =============================================================================
Δt = 0.02                    # Time step
T_final = 1.0                # Final time
nsteps = Int(T_final / Δt)   # Number of steps
reinit_freq = 10             # Reinitialization frequency

println("\n✓ Starting time evolution...")
println("  Δt = $Δt")
println("  T_final = $T_final")
println("  Total steps: $nsteps")
println("-" ^ 50)

for step in 1:nsteps
    # Advance the level set by Δt
    advance!(eg, Δt)
    
    # Periodic reinitialization to maintain signed distance property
    if step % reinit_freq == 0
        reinitialize!(eg)
        
        # Print progress
        t = current_time(eg)
        ϕ = current_levelset(eg)
        ϕ_min, ϕ_max = extrema(ϕ)
        println("  Step $step / $nsteps (t = $(round(t, digits=3)))")
        println("    Level set range: [$(round(ϕ_min, digits=4)), $(round(ϕ_max, digits=4))]")
    end
end

# =============================================================================
# 6. Final Results
# =============================================================================
println("-" ^ 50)
println("\n✓ Simulation completed!")
println("  Final time: $(current_time(eg))")

# Get grid information
info = grid_info(eg)
println("\nGrid info:")
println("  Origin: $(info.origin)")
println("  Spacing: $(info.spacing)")
println("  Node dimensions: $(info.dims)")

# Get level set statistics
ϕ_final = current_levelset(eg)
inside_mask = domain_mask(ϕ_final)
println("\nGeometry statistics:")
println("  Total nodes: $(length(ϕ_final))")
println("  Nodes inside domain (ϕ < 0): $(sum(inside_mask))")

# Access the cut geometry for CutFEM
cut_geo = current_cut(eg)
println("\n✓ Ready for CutFEM: current_cut(eg) available")
