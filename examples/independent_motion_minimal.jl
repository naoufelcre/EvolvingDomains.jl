# =============================================================================
# Independent Motion — Minimal Example
# =============================================================================
#
# Multi-level set pattern: Objects with independent velocities.
# Demonstrates: How to simulate objects passing through each other
# without forced merging (each has its own level set).
#
# =============================================================================

using EvolvingDomains
using Gridap
using LevelSetMethods
using CairoMakie

function run_independent_motion()
    println("═" ^ 60)
    println("  Independent Motion — EvolvingDomains.jl V2")
    println("═" ^ 60)
    
    # 1. Setup Grid
    domain = (0.0, 2.0, 0.0, 2.0)
    n = 60
    model = CartesianDiscreteModel(domain, (n, n))
    
    # 2. Two Objects with SEPARATE Level Sets
    # Pacman (left, moving right)
    pacman_disk = EvolvingDomains.Circle(VectorValue(0.0, 0.0), 0.3)
    pacman_mouth = EvolvingDomains.Circle(VectorValue(0.15, 0.0), 0.15)
    pacman = setdiff(pacman_disk, pacman_mouth)
    pacman_start = Translate(pacman, VectorValue(0.5, 1.0))
    
    # Dot (right, moving left)
    dot = EvolvingDomains.Circle(VectorValue(1.5, 1.0), 0.1)
    
    # 3. Create TWO Evolving Geometries (one per object)
    geom_pacman = EvolvingDiscreteGeometry(model, 
        x -> pacman_start(VectorValue(x...));
        reinit_freq = 5
    )
    
    geom_dot = EvolvingDiscreteGeometry(model, 
        x -> dot(VectorValue(x...));
        reinit_freq = 5
    )
    
    # 4. Independent Velocities
    set_velocity!(geom_pacman, StaticFunctionVelocity(x -> (0.5, 0.0)))  # Right
    set_velocity!(geom_dot, StaticFunctionVelocity(x -> (-0.5, 0.0)))     # Left
    
    # 5. Simulation Parameters
    t_end = 2.0
    dt = 0.04
    n_steps = round(Int, t_end / dt)
    
    println("Grid: $(n)×$(n), Δt = $dt, Steps = $n_steps")
    println("Pattern: Multi-Level Set (2 separate geometries)")
    println("")
    
    # 6. Record Animation
    fig = Figure(size=(600, 600))
    ax = Axis(fig[1, 1], aspect=1, limits=(0, 2, 0, 2),
              title="Independent Motion (Multi-LS)")
    
    output_path = joinpath(@__DIR__, "independent_motion_minimal.gif")
    
    record(fig, output_path, 1:n_steps; framerate=15) do step
        advance!(geom_pacman, dt)
        advance!(geom_dot, dt)
        
        empty!(ax)
        
        # Combine level sets for visualization (union)
        info = grid_info(geom_pacman)
        ϕ_pacman = current_levelset(geom_pacman)
        ϕ_dot = current_levelset(geom_dot)
        ϕ_combined = min.(ϕ_pacman, ϕ_dot)  # Union
        
        nx, ny = info.dims
        ϕ_2d = reshape(ϕ_combined, nx, ny)
        
        xs = range(info.origin[1], step=info.spacing[1], length=nx)
        ys = range(info.origin[2], step=info.spacing[2], length=ny)
        
        heatmap!(ax, xs, ys, ϕ_2d; colormap=:RdBu)
        contour!(ax, xs, ys, ϕ_2d; levels=[0.0], color=:black, linewidth=2)
        
        ax.title = "t = $(round(EvolvingDomains.current_time(geom_pacman), digits=2))"
    end
    
    println("✅ Animation saved to: $output_path")
    println("═" ^ 60)
end

# Run if executed directly
if abspath(PROGRAM_FILE) == @__FILE__
    run_independent_motion()
end
