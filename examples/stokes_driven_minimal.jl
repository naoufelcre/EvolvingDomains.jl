# =============================================================================
# Stokes-Driven Geometry — Minimal Example
# =============================================================================
#
# Demonstrates FEVelocitySource coupling: PDE-computed velocity drives
# geometry evolution. A circle is advected by a Stokes flow.
#
# This is the minimal working pattern for FSI-like simulations.
#
# =============================================================================

using EvolvingDomains
using Gridap
using GridapEmbedded
using LevelSetMethods
using CairoMakie
using LinearAlgebra: norm

# Inline plot_levelset! from extension to avoid loading issues
function plot_levelset!(ax, geom;
                        colormap=:RdBu,
                        levels::Int=20,
                        show_zero::Bool=true,
                        linewidth::Real=2,
                        linecolor=:black,
                        ϕ_buffer::Union{Nothing, Vector{Float64}}=nothing)
    # info = round_grid_info(geom) # Adjusted helper call if needed, or just extract directly
    # Original extension uses EvolvingDomains.grid_info(geom)
    # Since we are using .EvolvingDomains, we need to be careful with scope if copy-pasting
    
    # Just implement directly for this script:
    phi = EvolvingDomains.current_levelset(geom)
    info = EvolvingDomains.grid_info(geom)
    nx, ny = info.dims
    
    phi_2d = reshape(phi, nx, ny)
    
    ox, oy = info.origin
    dx, dy = info.spacing
    xs = range(ox, step=dx, length=nx)
    ys = range(oy, step=dy, length=ny)
    
    heatmap!(ax, xs, ys, phi_2d; colormap=colormap)
    if show_zero
        contour!(ax, xs, ys, phi_2d; levels=[0.0], color=linecolor, linewidth=linewidth)
    end
    return ax
end

function run_stokes_driven()
    println("═" ^ 60)
    println("  Stokes-Driven Geometry — EvolvingDomains.jl V2")
    println("  (With Quadtree + Guided Lookup Optimization)")
    println("═" ^ 60)
    
    # ==========================================================================
    # 1. Setup Quadtree Grid (Guided Lookup Requirement)
    # ==========================================================================
    println("[Grid] Generating Quadtree Mesh...")
    L = 6 # Resolution L=6
    
    # Minimal sizing: uniform refinement
    # To demonstrate GUIDED lookup, we actually just need A quadtree.
    # We will just refine everywhere to L=6 to mimic the previous 40x40 (which is close to 2^5=32 or 2^6=64)
    # 40x40 on length 2 => h = 0.05
    # L=5 => h = 1/32 = 0.031
    # L=6 => h = 1/64 = 0.015
    # Let's use L=5 to match ~40
    
    # We generate on [0,1] and map to [0,2]
    qmesh = generate_fine_mesh(1.0, 5) 
    balance!(qmesh)
    elements = pave_mesh(qmesh)
    
    model_unit, leaf_map = quadtree_to_discrete_model(elements)
    
    # Map to [0, 2]x[0, 2]
    # Gridap UnstructuredGrid allows node modification
    grid_unit = get_grid(model_unit)
    nodes = Gridap.Geometry.get_node_coordinates(grid_unit)
    new_nodes = [Gridap.Point(n[1]*2.0, n[2]*2.0) for n in nodes]
    
    grid_phys = Gridap.Geometry.UnstructuredGrid(new_nodes, 
                                 Gridap.Geometry.get_cell_node_ids(grid_unit), 
                                 Gridap.Geometry.get_reffes(grid_unit), 
                                 Gridap.Geometry.get_cell_type(grid_unit), 
                                 Gridap.Geometry.NonOriented())
    model = Gridap.Geometry.UnstructuredDiscreteModel(grid_phys)
    
    println("[Grid] Converted to Unstructured Model: $(num_cells(model)) cells")
    
    # We also need a scaling for the 'root' used in lookups?
    # GuidedVelocitySource uses 'root' which is in [0,1].
    # But checks point 'x' which is in [0,2].
    # PROBLEM: The Quadtree structure (in memory) thinks it covers [0,1].
    # If we pass physics point (1.5, 1.5), the tree search will fail or clamp.
    # SOLUTION: GuidedVelocitySource needs to know the domain map OR we map x -> x_unit before lookup.
    # Since GuidedVelocitySource is designed to be minimal, and `get_leaf_at` works on the tree coordinates...
    # We should probably modify GuidedVelocitySource to accept a coordinate transform?
    # OR simpler for this demo: Scale the problem to [0,1].
    # Let's scale the problem domain to [0,1].
    
    domain_len = 1.0 
    domain = (0.0, 1.0, 0.0, 1.0)
    
    # Re-generate Quadtree on [0,1] directly
    # qmesh acts on [0,1] by default size=1.0 center=(0.5,0.5)
    # So we don't need to remap nodes!
    model = model_unit # Just use the unit model directly
    
    println("[Grid] Domain is [0,1]x[0,1]")

    # ==========================================================================
    # 2. Create Evolving Geometry (circle)
    # ==========================================================================
    center = (0.5, 0.5) # Adjusted for [0,1] domain
    radius = 0.125     # Adjusted for [0,1] domain
    
    # Note: EvolvingDiscreteGeometry expects a Cartesian model?
    # NO. I checked `EvolvingDiscreteGeometry` constructor, it works with DiscreteModel 
    # IF the `lsm_model` (background for level set) is provided?
    # Actually, the minimal constructor `EvolvingDiscreteGeometry(model, dist)` 
    # might assume Cartesian if it builds the LSM grid automatically.
    # Let's check constructor signature... 
    # It usually creates a Cartesian background if not provided. Use explicit LSM grid for safety.
    
    lsm_model = CartesianDiscreteModel(domain, (60, 60)) # High res for geometry
    geom = EvolvingDiscreteGeometry(model, lsm_model, 
        x -> norm(x .- center) - radius;
        reinit_freq = 0
    )
    
    # ==========================================================================
    # 3. Create FE Velocity Field (Scalar Components)
    # ==========================================================================
    # We use Scalar FESpace to avoid Gridap Vector mapping artifacts
    reffe_scal = ReferenceFE(lagrangian, Float64, 1)
    V_scal = FESpace(model, reffe_scal)
    
    # Interpolate channel flow components: u_x = 4y(1-y), u_y = 0
    function flow_u(x)
        y = x[2]
        return 4.0 * y * (1.0 - y)
    end
    
    function flow_v(x)
        return 0.0
    end
    
    u_h_x = interpolate(flow_u, V_scal)
    u_h_y = interpolate(flow_v, V_scal)
    
    # ==========================================================================
    # 4. Guided Velocity Coupling
    # ==========================================================================
    # The Magic: GuidedVelocitySource with Lineage
    # Passing Tuple of Scalar Fields (x, y)
    vel = GuidedVelocitySource(u_h_x, u_h_y, qmesh.root, leaf_map)
    
    # Initial sync
    set_velocity!(geom, vel)
    
    println("Velocity: Parabolic channel flow (max magnitude ≈ 1.0)")
    
    # Verify non-zero velocity at interface
    v_sample = sample_velocity(vel, center, 0.0)
    println("Velocity at center: ($(round(v_sample[1], digits=3)), $(round(v_sample[2], digits=3)))")
    
    # ==========================================================================
    # 5. Simulation Parameters
    # ==========================================================================
    t_end = 0.75 # scaled time
    dt = 0.01
    n_steps = round(Int, t_end / dt)
    
    println("Steps: $n_steps, Δt = $dt")
    println("")
    
    # ==========================================================================
    # 6. Record Animation
    # ==========================================================================
    fig = Figure(size=(700, 600))
    ax = Axis(fig[1, 1], aspect=1, limits=(0, 1, 0, 1),
              xlabel="x", ylabel="y",
              title="Stokes-Driven (Guided Lookup) (t=0)")
    
    # Show velocity magnitude as background
    xs = range(0, 1, length=50)
    ys = range(0, 1, length=50)
    vel_mag = [norm(VectorValue(flow_u(VectorValue(x, y)), flow_v(VectorValue(x, y)))) for x in xs, y in ys]
    
    output_path = joinpath(@__DIR__, "stokes_driven_minimal.gif")
    
    record(fig, output_path, 1:n_steps; framerate=20) do step
        advance!(geom, dt)
        
        # In a real FSI, u_h changes. Here it is static.
        # But we demonstrate the update! API:
        # vel = update_velocity!(vel, u_h) # Returns NEW instance
        # set_velocity!(geom, vel)
        
        # Optimization: Since u_h didn't change, we don't strictly need to update,
        # but let's do it to show the pattern.
        # Pass Tuple (u_x, u_y)
        vel = update_velocity!(vel, (u_h_x, u_h_y))
        set_velocity!(geom, vel)
        
        empty!(ax)
        
        # Background: velocity magnitude
        heatmap!(ax, xs, ys, vel_mag; colormap=:blues, alpha=0.3)
        
        # Geometry contour
        plot_levelset!(ax, geom)
        
        t = EvolvingDomains.current_time(geom)
        ax.title = "Stokes-Driven (Guided Lookup) (t=$(round(t, digits=2)))"
    end
    
    println("✅ Animation saved to: $output_path")
    println("═" ^ 60)
    
    return geom
end

# Run if executed directly
if abspath(PROGRAM_FILE) == @__FILE__
    run_stokes_driven()
end
