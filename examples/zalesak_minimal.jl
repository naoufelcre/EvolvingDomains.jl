# =============================================================================
# Zalesak Disk Rotation — Minimal Example
# =============================================================================
#
# Classic benchmark: Zalesak disk (circle with slot) rotating 360°.
# Demonstrates: WENO5 advection, reinitialization, shape preservation.
#
# =============================================================================

using EvolvingDomains
using Gridap
using LevelSetMethods
using CairoMakie

function run_zalesak()
    println("═" ^ 60)
    println("  Zalesak Disk Rotation — EvolvingDomains.jl V2")
    println("═" ^ 60)
    
    # 1. Setup Grid
    domain = (-2.0, 2.0, -2.0, 2.0)
    n = 80
    model = CartesianDiscreteModel(domain, (n, n))
    
    # 2. Construct Zalesak Disk (CSG)
    # Circle centered at (0, 1) with slot cut out
    disk = EvolvingDomains.Circle(VectorValue(0.0, 0.0), 0.8)
    slot = EvolvingDomains.Rectangle(VectorValue(-0.15, 0.0), VectorValue(0.15, 1.0))
    zalesak = setdiff(disk, slot)
    zalesak_at_start = Translate(zalesak, VectorValue(0.0, 1.0))
    
    # 3. Create Evolving Geometry
    geom = EvolvingDiscreteGeometry(model, 
        x -> zalesak_at_start(VectorValue(x...));
        reinit_freq = 5
    )
    
    # 4. Set Velocity (Rigid Body Rotation)
    # One full revolution in T = 2π
    set_velocity!(geom, TimeDependentVelocity((x, t) -> (-x[2], x[1])))
    
    # 5. Simulation Parameters
    t_end = 2π
    dt = 0.05
    n_steps = round(Int, t_end / dt)
    
    println("Grid: $(n)×$(n), Δt = $dt, Steps = $n_steps")
    println("Target: Full 360° rotation")
    println("")
    
    # 6. Record Animation
    fig = Figure(size=(600, 600))
    ax = Axis(fig[1, 1], aspect=1, limits=(-2, 2, -2, 2),
              title="Zalesak Disk Rotation")
    
    output_path = joinpath(@__DIR__, "zalesak_minimal.gif")
    
    record(fig, output_path, 1:n_steps; framerate=20) do step
        advance!(geom, dt)
        empty!(ax)
        plot_levelset!(ax, geom)
        ax.title = "t = $(round(EvolvingDomains.current_time(geom), digits=2))"
    end
    
    println("✅ Animation saved to: $output_path")
    println("═" ^ 60)
end

# Run if executed directly
if abspath(PROGRAM_FILE) == @__FILE__
    run_zalesak()
end
