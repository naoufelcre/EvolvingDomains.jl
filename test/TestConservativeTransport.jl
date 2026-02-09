using Test
using EvolvingDomains
using EvolvingDomains.Geometric
using EvolvingDomains.Kinematic
using EvolvingDomains.Kinematic.SemiLagrangian
using Gridap
using Gridap.TensorValues
using LinearAlgebra

# For visualization
using CairoMakie
include("GeometryEvolutionEngine/Visualization.jl")
using .Visualization

@testset "Rigorous Rotating Checkerboard" begin

    # 1. Setup Parameters
    # -------------------
    n = 400
    domain = (0, 1, 0, 1)
    grid = CartesianDiscreteModel(domain, (n, n))
    vol_elem = (1.0/n)^2

    # High-resolution time-stepping
    nTime = 100
    Δt = 0.01
    plot_interval = 10 # 100 frames total

    center = VectorValue(0.4, 0.4)
    width, height = 0.25, 0.5

    # 2. Output Setup
    # ---------------
    output_dir = joinpath(@__DIR__, "output_conservative")
    mkpath(output_dir)
    println("Saving plots to: $output_dir")

    # 3. Define Initial Geometry (Rectangle)
    # --------------------------------------
    function rectangle_sdf(x)
        dx = abs(x[1] - center[1]) - width/2
        dy = abs(x[2] - center[2]) - height/2
        return max(dx, dy)
    end

    pts = Gridap.Geometry.get_node_coordinates(grid)
    phi_init = map(rectangle_sdf, pts)
    cart_info = grid_info(grid)

    geom = EvolvingDiscreteGeometry(vec(collect(phi_init)), grid)
    get_active_indices(geom, :current)
    geom.cache.prev_cut = geom.cache.cut

    # 4. Define Scalar Field (Checkerboard Pattern)
    # ---------------------------------------------
    k = 2 * (2π / width)
    function checkerboard(x)
        dx, dy = x[1] - center[1], x[2] - center[2]
        return sin(k*dx) * sin(k*dy)
    end

    # Field = 1.0 + sin(...) inside, 0.0 outside
    data = [phi <= 0 ? 1.0 + checkerboard(p) : 0.0 for (phi, p) in zip(phi_init, pts)]
    initial_field = CartesianMeshField(vec(data), cart_info)

    mass_init = sum(initial_field.data) * vol_elem

    # 5. Velocity
    # -----------
    ω = 2π
    vel = StaticFunctionVelocity(x -> VectorValue(-ω * (x[2] - 0.5), ω * (x[1] - 0.5)))
    vel_sampled = sample_velocity(vel, cart_info, 0.0)

    # 6. Helper: Plotting
    # -------------------
    function plot_field!(ax, geom, field; title="Transport", mass=0.0)
        nx, ny = field.grid.dims
        xs = range(field.grid.origin[1], step=field.grid.spacing[1], length=nx)
        ys = range(field.grid.origin[2], step=field.grid.spacing[2], length=ny)
        field_2d = reshape(field.data, nx, ny)
        heatmap!(ax, xs, ys, field_2d; colormap=:viridis, colorrange=(0, 2))
        ax.title = "$title\nMass: $(round(mass, digits=8))"
    end

    fig = Figure(size=(800, 800))
    ax = Axis(fig[1, 1], aspect=1)

    # 7. Time Loop
    # ------------
    current_field = initial_field
    checkpoint_cut = geom.cache.cut

    for step in 1:nTime
        # Capture geom before update for diagnostics
        if (step-1) % plot_interval == 0
             get_active_indices(geom, :current)
             checkpoint_cut = geom.cache.cut
        end

        # Evolve
        advance!(geom, vel_sampled, Δt)
        current_field = advect_conservative(geom, current_field, vel, Δt)

        # Plotting
        if step % plot_interval == 0
            t = step * Δt
            current_mass = sum(current_field.data) * vol_elem
            mass_err = current_mass - mass_init
            println("Step $step (t=$(round(t, digits=3))) | Mass: $(round(current_mass, digits=10)) | Err: $mass_err")

            # Update cache history for visualization tools
            geom.cache.prev_cut = checkpoint_cut

            empty!(ax)
            plot_field!(ax, geom, current_field; title="Transport t=$(round(t, digits=3))", mass=current_mass)
            save(joinpath(output_dir, "frame_$(lpad(step, 6, '0')).png"), fig)
        end
    end

    # Final Verification
    mass_final = sum(current_field.data) * vol_elem
    @test abs(mass_final - mass_init) < 1e-12
    println("Final Mass Error: $(abs(mass_final - mass_init))")
end
