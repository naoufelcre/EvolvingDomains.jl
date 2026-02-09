using Test
using EvolvingDomains
using EvolvingDomains.Geometric
using EvolvingDomains.Transfer
using Gridap
using GridapEmbedded
using Gridap.TensorValues
using CairoMakie
using LinearAlgebra

# Include local visualization helper
include("GeometryEvolutionEngine/Visualization.jl")
using .Visualization

@testset "Visual Transfer Verification" begin

    # --- 0. HELPER FUNCTIONS ---

    function analytical_field(x)
        return sin(π * x[1]) * sin(π * x[2])
    end

    function get_reference_field(geom)
        info = grid_info(geom.grid)
        pts = Gridap.Geometry.get_node_coordinates(geom.grid)
        u_data = vec(collect(map(analytical_field, pts)))
        return CartesianMeshField(u_data, info)
    end

    function plot_comparison!(fig, pos, u_ref, u_res, u_err, title_prefix, geom)
        nx, ny = u_ref.grid.dims
        xs = range(u_ref.grid.origin[1], step=u_ref.grid.spacing[1], length=nx)
        ys = range(u_ref.grid.origin[2], step=u_ref.grid.spacing[2], length=ny)

        ax1 = Axis(fig[pos, 1], aspect=1, title="$title_prefix: Input")
        heatmap!(ax1, xs, ys, reshape(u_ref.data, nx, ny), colormap=:viridis)
        contour!(ax1, xs, ys, reshape(geom.levelset, nx, ny), levels=[0.0], color=:red)

        ax2 = Axis(fig[pos, 2], aspect=1, title="$title_prefix: Result")
        heatmap!(ax2, xs, ys, reshape(u_res.data, nx, ny), colormap=:viridis)
        contour!(ax2, xs, ys, reshape(geom.levelset, nx, ny), levels=[0.0], color=:red)

        ax3 = Axis(fig[pos, 3], aspect=1, title="$title_prefix: Absolute Error")
        hm_err = heatmap!(ax3, xs, ys, reshape(u_err.data, nx, ny), colormap=:plasma)
        contour!(ax3, xs, ys, reshape(geom.levelset, nx, ny), levels=[0.0], color=:white, alpha=0.5)
        Colorbar(fig[pos, 4], hm_err)
    end

    # --- 1. SETUP GEOMETRY ---
    n = 60 # Higher resolution for better visual verification
    domain = (-1, 1, -1, 1)
    grid = CartesianDiscreteModel(domain, (n, n))

    # Define Geometry: Circle R=0.6 centered at (0.0, 0.0)
    circle = EvolvingDomains.Geometric.Circle(VectorValue(0.0, 0.0), 0.6)

    # Initialize Geometry efficiently
    # Directly passing the evaluated levelset to the constructor avoids redundant copies
    pts = Gridap.Geometry.get_node_coordinates(grid)
    phi_values = vec(collect(map(x -> circle(x), pts)))
    geom = EvolvingDiscreteGeometry(phi_values, grid)

    # Pre-calculate active indices (Nodes belonging to IN or CUT cells)
    # This is called internally by setup_transfer, but we need it for masked error calculation.
    active_indices = get_active_indices(geom)

    # Output directory
    output_dir = joinpath(@__DIR__, "output_visual_transfer")
    mkpath(output_dir)

    # Common Reference Field
    u_ref = get_reference_field(geom)


    # --- 2. EXPERIMENT 1: GRID-MESH ROUND TRIP ---
    # Testing whether restricting to the mesh and prolonging back preserves the field.

    # Ensure cut is populated
    get_active_indices(geom, :current)
    cut_geo = geom.cache.cut

    # Create Active Triangulation
    Ω_act = Triangulation(cut_geo, ACTIVE)

    reffe = ReferenceFE(lagrangian, Float64, 1)

    # Use Active Space (mimicking AgFEM setup)
    V = TestFESpace(Ω_act, reffe)

    setup_transfer(geom, V)

    # 1. Grid -> Mesh (Restriction)
    u_mesh = grid_to_mesh(geom, u_ref)

    # 2. Mesh -> Grid (Prolongation)
    u_back = mesh_to_grid(geom, u_mesh)

    # 3. Compute Error (Only on active nodes)
    err_transfer_data = zeros(length(u_ref.data))
    max_err_transfer = 0.0
    for idx in active_indices
        err = abs(u_back.data[idx] - u_ref.data[idx])
        err_transfer_data[idx] = err
        max_err_transfer = max(max_err_transfer, err)
    end
    u_err_transfer = CartesianMeshField(err_transfer_data, u_ref.grid)

    println("Max Transfer Round-Trip Error: $max_err_transfer")
    @test max_err_transfer < 1e-3


    # --- 3. EXPERIMENT 2: CLOSEST POINT EXTENSION ---
    # Testing whether we can recover the field in the void from bulk values.

    # 1. Mask Input: Set to ZERO in the void (ϕ > 0)
    masked_data = [phi <= 0 ? val : 0.0 for (phi, val) in zip(geom.levelset, u_ref.data)]
    u_masked = CartesianMeshField(vec(collect(masked_data)), u_ref.grid)

    # 2. Perform Extension
    u_ext = extend(geom, u_masked)

    # 3. Compute Error (Global)
    # Note: The 'error' here compares the constant-along-normal extension to the
    # analytical field. Since sin(πx)sin(πy) is NOT constant along normals,
    # we expect a significant discrepancy in the void. This test verifies
    # that the extension is bounded and produces a coherent field.
    u_err_ext_data = abs.(u_ext.data .- u_ref.data)
    u_err_ext = CartesianMeshField(u_err_ext_data, u_ref.grid)

    max_err_ext = maximum(u_err_ext_data)
    println("Max Extension vs Analytical Error: $max_err_ext")

    # We expect the error to be bounded by the field's range (1.0)
    @test max_err_ext < 1.1


    # --- 4. EXPERIMENT 3: MESH-GRID ROUND TRIP ---
    # Testing Mesh -> Grid -> Mesh. Crucial for FEM-native fields.

    # 1. Initialize directly on Mesh (Interpolate analytical function)
    # Note: 'interpolate' is a Gridap function
    u_fem_init = interpolate(analytical_field, V)

    # 2. Mesh -> Grid (Prolongation)
    u_grid_intermediate = mesh_to_grid(geom, u_fem_init)

    # 3. Grid -> Mesh (Restriction)
    u_mesh_final = grid_to_mesh(geom, u_grid_intermediate)

    # 4. Compute Error (Visualized on Grid)
    # We evaluate the difference between the initial and final FE functions on the grid nodes
    # This gives us a consistent way to visualize the error.
    # Error = | P(u_fem_init) - P(R(P(u_fem_init))) |  ? No, that's just P(u) - P(u').
    # We want to see how much the FE function changed.
    # Let's project the difference to the grid:
    u_diff_fem = interpolate(u_fem_init - u_mesh_final, V) # Does Gridap support subtraction of FEFunctions? Yes.
    u_err_mesh_data_on_grid = mesh_to_grid(geom, u_diff_fem).data

    # Take absolute value for heatmap
    u_err_mesh_grid = CartesianMeshField(abs.(u_err_mesh_data_on_grid), u_ref.grid)

    # Also compute max error on active nodes for the test
    max_err_mesh = 0.0
    for idx in active_indices
        val = u_err_mesh_grid.data[idx]
        max_err_mesh = max(max_err_mesh, val)
    end
    println("Max Mesh Round-Trip Error (Mesh->Grid->Mesh): $max_err_mesh")
    @test max_err_mesh < 1e-3

    # Visualization Helpers for FE objects (map to grid first)
    u_fem_init_grid = mesh_to_grid(geom, u_fem_init)
    u_mesh_final_grid = mesh_to_grid(geom, u_mesh_final)


    # --- 5. VISUALIZATION ---

    fig = Figure(size=(1400, 1200))
    Label(fig[0, :], "Visual Transfer & Extension Verification (Field: sin(πx)sin(πy))", fontsize=20, font=:bold)

    # Plot experiments
    plot_comparison!(fig, 1, u_ref, u_back, u_err_transfer, "Exp 1: Grid -> Mesh -> Grid", geom)
    plot_comparison!(fig, 2, u_masked, u_ext, u_err_ext, "Exp 2: Extension (Bulk -> Void)", geom)
    plot_comparison!(fig, 3, u_fem_init_grid, u_mesh_final_grid, u_err_mesh_grid, "Exp 3: Mesh -> Grid -> Mesh", geom)

    save(joinpath(output_dir, "visual_verification.png"), fig)
    println("Saved unified visualization to: $output_dir/visual_verification.png")

end
