```@meta
EditURL = "../../../examples/zalesak_minimal.jl"
```

=============================================================================
Zalesak Disk Rotation — Minimal Example
=============================================================================

Classic benchmark: Zalesak disk (circle with slot) rotating 360°.
Demonstrates: WENO5 advection, reinitialization, shape preservation.

=============================================================================

````@example zalesak_minimal
using EvolvingDomains
using Gridap
using LevelSetMethods
using CairoMakie

println("═" ^ 60)
println("  Zalesak Disk Rotation — EvolvingDomains.jl V2")
println("═" ^ 60)
````

1. Setup Grid

````@example zalesak_minimal
domain = (-2.0, 2.0, -2.0, 2.0)
n = 80
model = CartesianDiscreteModel(domain, (n, n))
````

2. Construct Zalesak Disk (CSG)
Circle centered at (0, 1) with slot cut out

````@example zalesak_minimal
disk = EvolvingDomains.Circle(VectorValue(0.0, 0.0), 0.8)
slot = EvolvingDomains.Rectangle(VectorValue(-0.15, 0.0), VectorValue(0.15, 1.0))
zalesak = setdiff(disk, slot)
zalesak_at_start = Translate(zalesak, VectorValue(0.0, 1.0))
````

3. Create Evolving Geometry

````@example zalesak_minimal
geom = EvolvingDiscreteGeometry(model,
    x -> zalesak_at_start(VectorValue(x...));
    reinit_freq = 5
)
````

4. Set Velocity (Rigid Body Rotation)
One full revolution in T = 2π

````@example zalesak_minimal
set_velocity!(geom, TimeDependentVelocity((x, t) -> (-x[2], x[1])))
````

5. Simulation Parameters

````@example zalesak_minimal
t_end = 2π
dt = 0.05
n_steps = round(Int, t_end / dt)

println("Grid: $(n)×$(n), Δt = $dt, Steps = $n_steps")
println("Target: Full 360° rotation")
println("")
````

6. Record Animation

````@example zalesak_minimal
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
````

![Zalesak Disk](zalesak_minimal.gif)

---

*This page was generated using [Literate.jl](https://github.com/fredrikekre/Literate.jl).*

