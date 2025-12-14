# =============================================================================
# Example: Using the Snapshot Visualization API
# =============================================================================
# Demonstrates the new snapshot-based visualization workflow.
#
# Usage:
#   julia --project=. examples/poc_live_animation.jl

using Pkg
Pkg.activate(dirname(@__DIR__))

using EvolvingDomains
using Gridap, GridapEmbedded
using LevelSetMethods
using GLMakie

# =============================================================================
# Setup
# =============================================================================
domain = (-1.0, 1.0, -1.0, 1.0)
partition = (50, 50)
model = CartesianDiscreteModel(domain, partition)

# Circle off-center so rotation is visible
R = 0.2
center = (0.5, 0.0)
ϕ₀(x) = R - sqrt((x[1]-center[1])^2 + (x[2]-center[2])^2)

# Rigid body rotation
ω = 2π
u(x) = (-ω * x[2], ω * x[1])

evolver = LevelSetMethodsEvolver(;
    bg_model = model,
    initial_ls = ϕ₀,
    velocity = u,
    spatial_scheme = :WENO5,
    bc = :Neumann
)
eg = EvolvingDiscreteGeometry(evolver, model)

# =============================================================================
# Run Simulation with Snapshotting
# =============================================================================
println("=" ^ 50)
println("Snapshot Visualization API Demo")
println("=" ^ 50)

frames = SimulationFrame[]
push!(frames, snapshot(eg))  # Initial state

Δt = 0.01
nsteps = 100
cache_every = 2

println("Running simulation...")
for step in 1:nsteps
    # Advance level set
    EvolvingDomains.advance!(eg, Δt)
    
    # Reinitialize periodically
    if step % 20 == 0
        EvolvingDomains.reinitialize!(eg)
    end
    
    # Snapshot when desired
    if step % cache_every == 0
        push!(frames, snapshot(eg))
    end
    
    if step % 25 == 0
        println("  Step $step/$nsteps")
    end
end

println("✓ Simulation complete: $(length(frames)) frames cached")

# =============================================================================
# Create SimulationResult and View
# =============================================================================
result = SimulationResult(grid_info(eg), frames)

println("\nOpening interactive viewer...")
println("Drag slider to navigate through time")

fig = viewer(result; wait=true)

println("Done!")
