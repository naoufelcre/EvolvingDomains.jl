using Test
using EvolvingDomains
using EvolvingDomains.Geometric
using EvolvingDomains.Kinematic
using Gridap
using Gridap.TensorValues

# For visualization
using CairoMakie
# Extension loaded automatically

# Include local test engine tools
include("GeometryEvolutionEngine/Visualization.jl")
using .Visualization
using EvolvingDomains.Kinematic: compute_raymaps, RayMaps

@testset "Zalesak Disk" begin
    # Define Grid (100x100 on [0,1]x[0,1])
    n = 100
    grid = CartesianDiscreteModel((0, 1, 0, 1), (n, n))
    center = VectorValue(0.5, 0.5) # Center of the Grid

    # Create empty geometry container
    geom = EvolvingDiscreteGeometry(grid)

    # Define Geometry (Zalesak Disk)
    disk = EvolvingDomains.Geometric.Circle(VectorValue(0.5, 0.75), 0.15)
    slot = EvolvingDomains.Geometric.Rectangle(VectorValue(0.475, 0.60), VectorValue(0.525, 0.75))
    zalesak_geo = setdiff(disk, slot)

    # Initialize Level Set Field
    pts = Gridap.Geometry.get_node_coordinates(grid)
    phi_values = map(x -> zalesak_geo(x), pts)
    set_levelset!(geom, vec(collect(phi_values)))

    # Define Velocity as Solid Body Rotation.
    ω = 20π # Angular velocity
    function rigid_rotation(x)
        rx = x[1] - center[1]
        ry = x[2] - center[2]
        VectorValue(-ω * ry, ω * rx)
    end
    #Set up
    vel = StaticFunctionVelocity(rigid_rotation)

    # 5. Output Setup
    output_dir = joinpath(@__DIR__, "output")
    mkpath(output_dir)

    # Force initialization of the cut for the first step
    # This ensures 'prev_cut' is populated correctly when we advance! for the first time
    get_active_indices(geom, :current)

    # ====== Initial State Visualization =====
    # Classification and Level Set Values
    fig = Figure(size=(800, 800))
    ax = Axis(fig[1, 1], aspect=1, title="Zalesak Disk t=0.0")
    Visualization.plot_cell_classification!(ax, geom)
    save(joinpath(output_dir, "frame_000000_classification.png"), fig)
    EvolvingDomains.plot_levelset!(ax, geom; show_zero=true, linecolor=:black, linewidth=3)
    save(joinpath(output_dir, "frame_000000_levelset.png"), fig)

    # ====== Figure for velocity field. ======
    fig_vel = Figure(size=(800, 800))
    ax_vel = Axis(fig_vel[1, 1], aspect=1, title="Velocity Field (Rigid Rotation)")

    # Extract x and y components from Gridap points
    # pts was defined in step 3
    xs = [p[1] for p in pts]
    ys = [p[2] for p in pts]

    # Calculate velocity vectors at every point
    vs_vec = rigid_rotation.(pts)
    us = [v[1] for v in vs_vec]
    vs = [v[2] for v in vs_vec]

    # Calculate magnitude for coloring (optional but looks nice)
    magnitude = sqrt.(us .^ 2 .+ vs .^ 2)

    # ! IMPORTANT: 100x100 = 10,000 points. This is too dense for an arrow plot.
    # Create a mask or indices for subsampling.
    subsampling_index = 32
    indices = 1:subsampling_index:length(pts)

    arrows!(ax_vel,
        xs[indices], ys[indices],
        us[indices], vs[indices];
        arrowsize=10,
        lengthscale=1.0, # Actual velocity vectors
        linecolor=magnitude[indices],
        colormap=:viridis
    )

    x_grid_lines = unique(xs) # Returns the 101 unique x-coords
    y_grid_lines = unique(ys) # Returns the 101 unique y-coords

    vlines!(ax_vel, x_grid_lines, color=:black, linewidth=0.5, alpha=0.5)
    hlines!(ax_vel, y_grid_lines, color=:black, linewidth=0.5, alpha=0.5)
    save(joinpath(output_dir, "frame_000_velocity.png"), fig_vel)

    # ====== Fin Initial state vizualization ======

    # ====== Time Loop ======

    nTime = 1000000
    timeHorizon = 1
    Δt = timeHorizon / nTime

    # Diagnostic parameters
    diag_stride = 100
    Δt_diag = Δt * diag_stride # Time step for visualization (to make arrows visible)

    info = grid_info(grid)
    v_field = sample_velocity(vel, info, 0)

    # Run full simulation
    checkpoint_cut = nothing
    for i in 1:10000
        # At the very start of a diagnostic stride, capture the geometry
        if (i - 1) % diag_stride == 0
            get_active_indices(geom, :current)
            checkpoint_cut = geom.cache.cut
        end

        t = i * Δt

        # Evolve geometry
        advance!(geom, v_field, Δt)

        # Diagnostics
        if i % diag_stride == 0
            println("Step $i / $nTime (t=$t)")

            # Inject the checkpoint from the start of the stride into the cache
            # so that state=:prev correctly references the geometry from 10,000 steps ago.
            geom.cache.prev_cut = checkpoint_cut

            # 1. Geometry Figure
            fig_geo = Figure(size=(800, 800))
            ax_geo = Axis(fig_geo[1, 1], aspect=1, title="Zalesak Disk t=$(round(t, digits=3))")
            Visualization.plot_cell_classification!(ax_geo, geom)
            EvolvingDomains.plot_levelset!(ax_geo, geom; show_zero=true, linecolor=:black, linewidth=3)
            save(joinpath(output_dir, "frame_$(lpad(i, 8, '0'))_geometry.png"), fig_geo)
        end
    end
end
