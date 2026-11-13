# =============================================================================
# PublicationFigures.jl — Publication-ready B&W vector figures
# =============================================================================
# Generates individual PDF panels for three test cases using grayscale
# heatmaps, discrete contour lines, and 1-D cross-sections.
#
# Run: julia test/PublicationFigures.jl
# =============================================================================

using EvolvingDomains
using EvolvingDomains.Geometric
using EvolvingDomains.Kinematic
using EvolvingDomains.Kinematic.SemiLagrangian
using EvolvingDomains.Transfer
using Gridap
using Gridap.Geometry
using Gridap.TensorValues
using GridapEmbedded
using CairoMakie
using LinearAlgebra

# ============================================================
# OUTPUT SETUP
# ============================================================
const output_dir = joinpath(@__DIR__, "output_publication")
mkpath(output_dir)

# ============================================================
# HELPERS
# ============================================================
"""Bilinear interpolation of a CartesianMeshField at point (x,y)."""
function sample_field(field, x, y)
    info = field.grid
    nx, ny = info.dims
    ox, oy = info.origin
    dx, dy = info.spacing

    i = floor(Int, (x - ox) / dx) + 1
    j = floor(Int, (y - oy) / dy) + 1
    i = clamp(i, 1, nx - 1)
    j = clamp(j, 1, ny - 1)

    xi  = clamp((x - (ox + (i - 1) * dx)) / dx, 0.0, 1.0)
    eta = clamp((y - (oy + (j - 1) * dy)) / dy, 0.0, 1.0)

    idx(ii, jj) = ii + (jj - 1) * nx
    v11 = field.data[idx(i,   j)]
    v21 = field.data[idx(i+1, j)]
    v12 = field.data[idx(i,   j+1)]
    v22 = field.data[idx(i+1, j+1)]

    return (1 - xi) * (1 - eta) * v11 + xi * (1 - eta) * v21 +
           (1 - xi) * eta       * v12 + xi * eta       * v22
end

"""Interpolate a CartesianMeshField to a finer regular grid."""
function interpolate_to_fine_grid(field, fine_nx, fine_ny, bounds)
    info = field.grid
    ox, oy = info.origin
    dx, dy = info.spacing
    nx, ny = info.dims

    xmin, xmax, ymin, ymax = bounds
    fine_dx = (xmax - xmin) / (fine_nx - 1)
    fine_dy = (ymax - ymin) / (fine_ny - 1)

    fine_data = fill(NaN, fine_nx, fine_ny)

    for j in 1:fine_ny, i in 1:fine_nx
        x = xmin + (i - 1) * fine_dx
        y = ymin + (j - 1) * fine_dy

        ci = (x - ox) / dx + 1
        cj = (y - oy) / dy + 1
        i0 = floor(Int, ci)
        j0 = floor(Int, cj)

        if i0 < 1 || i0 >= nx || j0 < 1 || j0 >= ny
            continue
        end

        xi = ci - i0
        eta = cj - j0

        idx(ii, jj) = ii + (jj - 1) * nx
        v11 = field.data[idx(i0,   j0)]
        v21 = field.data[idx(i0+1, j0)]
        v12 = field.data[idx(i0,   j0+1)]
        v22 = field.data[idx(i0+1, j0+1)]

        if any(isnan, (v11, v21, v12, v22))
            continue
        end

        fine_data[i, j] = (1 - xi) * (1 - eta) * v11 +
                            xi       * (1 - eta) * v21 +
                            (1 - xi) * eta       * v12 +
                            xi       * eta       * v22
    end

    xs = range(xmin, step=fine_dx, length=fine_nx)
    ys = range(ymin, step=fine_dy, length=fine_ny)
    return xs, ys, fine_data
end


# ============================================================
# PART 1: CONSERVATIVE TRANSPORT  (Checkerboard)
# ============================================================
function generate_conservative_figures()
    n = 500
    domain = (0, 1, 0, 1)
    grid = CartesianDiscreteModel(domain, (n, n))
    vol_elem = (1.0 / n)^2

    nTime = 100
    Δt = 0.01

    center = VectorValue(0.4, 0.4)
    width, height = 0.25, 0.5

    function rectangle_sdf(x)
        dx = abs(x[1] - center[1]) - width / 2
        dy = abs(x[2] - center[2]) - height / 2
        return max(dx, dy)
    end

    pts = Gridap.Geometry.get_node_coordinates(grid)
    phi_init = map(rectangle_sdf, pts)
    cart_info = grid_info(grid)

    geom = EvolvingDiscreteGeometry(vec(collect(phi_init)), grid)
    get_active_indices(geom, :current)
    geom.cache.prev_cut = geom.cache.cut

    k = 2 * (2π / width)
    function checkerboard(x)
        dx, dy = x[1] - center[1], x[2] - center[2]
        return sin(k * dx) * sin(k * dy)
    end

    data = [phi <= 0 ? 1.0 + checkerboard(p) : 0.0 for (phi, p) in zip(phi_init, pts)]
    initial_field = CartesianMeshField(vec(data), cart_info)

    ω = 2π
    vel = StaticFunctionVelocity(x -> VectorValue(-ω * (x[2] - 0.5), ω * (x[1] - 0.5)))
    vel_sampled = sample_velocity(vel, cart_info, 0.0)

    current_field = initial_field
    mass_init = sum(initial_field.data) * vol_elem

    mass_history = Float64[mass_init]
    time_history = Float64[0.0]

    # Exact rotated-rectangle SDF (used for masking & contour)
    function exact_rect_sdf(x, y, cx, cy, w, h, θ)
        dx, dy = x - cx, y - cy
        rx =  dx * cos(θ) + dy * sin(θ)
        ry = -dx * sin(θ) + dy * cos(θ)
        return max(abs(rx) - w / 2, abs(ry) - h / 2)
    end

    # Snapshot helper ---------------------------------------------------------
    function save_conservative_snapshot(step, t, geom, field)
        θ = ω * t
        rel_x, rel_y = -0.1, -0.1
        c_x = 0.5 + rel_x * cos(θ) - rel_y * sin(θ)
        c_y = 0.5 + rel_x * sin(θ) + rel_y * cos(θ)

        # ---- 2-D contour plot ----
        fig_2d = Figure(size=(600, 600))
        ax_2d = Axis(fig_2d[1, 1], aspect=1,
                     title="Transport t = $t",
                     xlabel="x", ylabel="y",
                     xgridvisible=false, ygridvisible=false)

        nx, ny = field.grid.dims
        xs = range(field.grid.origin[1], step=field.grid.spacing[1], length=nx)
        ys = range(field.grid.origin[2], step=field.grid.spacing[2], length=ny)
        field_2d = reshape(field.data, nx, ny)

        # Mask with exact analytical SDF
        field_masked = copy(field_2d)
        for j in 1:ny, i in 1:nx
            if exact_rect_sdf(xs[i], ys[j], c_x, c_y, width, height, θ) > 0
                field_masked[i, j] = NaN
            end
        end

        # 4 discrete bands: 0–0.5, 0.5–1.0, 1.0–1.5, 1.5–2.0
        heatmap!(ax_2d, xs, ys, field_masked;
                 colormap=cgrad([:white, :gray85, :gray50, :gray20], 4, categorical=true),
                 colorrange=(0, 2))

        # Interface contour: use the *exact* analytical SDF evaluated on the grid
        # (the evolved numerical levelset phi_2d degenerates after t=0)
        exact_phi_2d = [exact_rect_sdf(x, y, c_x, c_y, width, height, θ) for x in xs, y in ys]
        contour!(ax_2d, xs, ys, exact_phi_2d;
                 levels=[0.0], color=:gray, linewidth=1.5, linestyle=:dash)

        xlims!(ax_2d, 0, 1)
        ylims!(ax_2d, 0, 1)

        save(joinpath(output_dir, "conservative_2d_t$(step).pdf"), fig_2d)

        # ---- 1-D slab in the rotating checkerboard frame (local y′ = 0.1) ----
        fig_1d = Figure(size=(700, 400))
        ax_1d = Axis(fig_1d[1, 1],
                     title="Cross-section (t = $t)",
                     xlabel="Local x′ (along width)",
                     ylabel="Field value",
                     xgridvisible=false, ygridvisible=false)

        n_samples = 600
        local_y_fixed = 0.1
        # Keep the sample line strictly inside the rectangle so that bilinear
        # interpolation never mixes interior values with the exterior zeros.
        margin = 2 * field.grid.spacing[1]
        local_xs = range(-width / 2 + margin, width / 2 - margin, length=n_samples)

        eulerian_xs = [c_x + lx * cos(θ) - local_y_fixed * sin(θ)
                       for lx in local_xs]
        eulerian_ys = [c_y + lx * sin(θ) + local_y_fixed * cos(θ)
                       for lx in local_xs]

        # Analytical reference (valid for all t in the local frame for this test)
        k = 2 * (2π / width)
        C = sin(k * local_y_fixed)
        exact_vals = [1.0 + sin(k * lx) * C for lx in local_xs]

        # Numerical sampling via bilinear interpolation
        field_vals = [sample_field(field, ex, ey)
                      for (ex, ey) in zip(eulerian_xs, eulerian_ys)]

        lines!(ax_1d, local_xs, field_vals, color=:black, linewidth=1.5, label="Numerical")
        lines!(ax_1d, local_xs, exact_vals, color=:gray, linewidth=0.8, linestyle=:dash, label="Analytical")

        # reference contour levels
        hlines!(ax_1d, [0.5, 1.0, 1.5], color=:gray,
                linewidth=0.5, linestyle=:dash)

        # rectangle edges
        vlines!(ax_1d, [-width/2, width/2], color=:gray,
                linewidth=0.8, linestyle=:dot)

        ylims!(ax_1d, -0.2, 2.2)
        xlims!(ax_1d, -width/2 - 0.03, width/2 + 0.03)

        save(joinpath(output_dir, "conservative_1d_t$(step).pdf"), fig_1d)
    end

    # t = 0
    save_conservative_snapshot(0, 0.0, geom, current_field)

    for step in 1:nTime
        advance!(geom, vel_sampled, Δt)
        k_map = TransportMap(geom, vel, Δt)
        new_data = zeros(Float64, length(current_field.data))
        advect!(new_data, current_field.data, k_map)
        current_field = CartesianMeshField(new_data, current_field.grid)

        current_mass = sum(current_field.data) * vol_elem
        push!(mass_history, current_mass)
        push!(time_history, step * Δt)

        if step in [33, 66, 100]
            save_conservative_snapshot(step, step * Δt, geom, current_field)
        end
    end

    # ---- Mass-evolution plot (relative error) ----
    mass_err = abs.(mass_history .- mass_init)
    fig_mass = Figure(size=(600, 400))
    ax_mass = Axis(fig_mass[1, 1],
                   title="Mass conservation error",
                   xlabel="Time", ylabel="|Mass − Mass₀|",
                   xgridvisible=false, ygridvisible=false)
    lines!(ax_mass, time_history, mass_err,
           color=:black, linewidth=1.5)
    # scientific notation often needs manual y-limits for tiny numbers
    ylims!(ax_mass, 0, max(maximum(mass_err) * 1.2, 1e-14))
    save(joinpath(output_dir, "conservative_mass.pdf"), fig_mass)
end


# ============================================================
# PART 2: DUMBBELL TEMPERATURE
# ============================================================
function save_dumbbell_frame(step, t, geom, u_grid, info)
    fig = Figure(size=(700, 600))
    ax = Axis(fig[1, 1], aspect=DataAspect(),
              title="Temperature t = $(round(t, digits=3))",
              xlabel="x", ylabel="y",
              xgridvisible=false, ygridvisible=false)

    # Fine grid for smooth rendering
    fine_n = 300
    bounds = (-2.0, 2.0, -2.0, 2.0)
    xs_fine, ys_fine, u_fine = interpolate_to_fine_grid(u_grid, fine_n, fine_n, bounds)

    # mask with analytical levelset evaluated on the fine grid
    for j in 1:fine_n, i in 1:fine_n
        x = xs_fine[i]
        y = ys_fine[j]
        phi = (t - 0.25) - ((x^2 - 0.3 * x^4) - y^2)
        if phi > 0
            u_fine[i, j] = NaN
        end
    end

    # grayscale: white → gray50  (high temp stays distinguishable from black contour)
    hm = heatmap!(ax, xs_fine, ys_fine, u_fine,
                  colormap=cgrad([:white, :gray50]), colorrange=(0, 1))

    # interface contour on the same fine grid
    phi_fine = [ (t - 0.25) - ((x^2 - 0.3 * x^4) - y^2)
                 for x in xs_fine, y in ys_fine ]
    contour!(ax, xs_fine, ys_fine, phi_fine,
             levels=[0.0], color=:black, linewidth=2)

    Colorbar(fig[1, 2], hm)

    save(joinpath(output_dir, "dumbbell_t$(step).pdf"), fig)
end

function generate_dumbbell_figures()
    domain = (-2, 2, -2.0, 2.0)
    n = 75
    grid = CartesianDiscreteModel(domain, (n, n))
    info = grid_info(grid)

    function analytical_levelset(t, x)
        x1, x2 = x[1], x[2]
        return (t - 0.25) - ((x1^2 - 0.3 * x1^4) - x2^2)
    end

    function u0_func(x)
        x1, x2 = x[1], x[2]
        return analytical_levelset(0, x)^2 + (x1 + x2 + 1)
    end

    function heat_varforms(u_prev, Δt, dΩ)
        a(u, v) = ∫((u * v) / Δt + (∇(u) ⋅ ∇(v)))dΩ
        l(v) = ∫((u_prev * v) / Δt)dΩ
        return a, l
    end

    t = 0.0
    pts = Gridap.Geometry.get_node_coordinates(grid)
    phi_0 = map(x -> analytical_levelset(t, x), pts)
    u0    = map(x -> u0_func(x), pts)

    geom   = EvolvingDiscreteGeometry(vec(collect(phi_0)), grid)
    u_grid = CartesianMeshField(vec(collect(u0)), info)

    nTime = 200
    timeHorizon = 0.5
    Δt = timeHorizon / nTime

    save_dumbbell_frame(0, 0.0, geom, u_grid, info)

    for step in 1:nTime
        t += Δt
        phi_new = map(x -> analytical_levelset(t, x), pts)
        set_levelset!(geom, vec(collect(phi_new)))
        get_active_indices(geom, :current)
        cutgeo = geom.cache.cut

        Ω     = Triangulation(cutgeo, PHYSICAL)
        Ω_act = Triangulation(cutgeo, ACTIVE)

        order = 1
        reffe = ReferenceFE(lagrangian, Float64, order)
        Vstd  = TestFESpace(Ω_act, reffe, conformity=:H1)

        strategy   = AggregateAllCutCells()
        aggregates = aggregate(strategy, cutgeo)

        V = AgFEMSpace(Vstd, aggregates)
        U = TrialFESpace(V)

        degree = 2 * order
        dΩ = Measure(Ω, degree)

        setup_transfer(geom, V)
        u_prev = grid_to_mesh(geom, u_grid)

        a, l = heat_varforms(u_prev, Δt, dΩ)
        op = AffineFEOperator(a, l, U, V)
        u_h = solve(op)

        u_grid = mesh_to_grid(geom, u_h)

        if step in [50, 100, 150, 200]
            save_dumbbell_frame(step, t, geom, u_grid, info)
        end
    end
end


# ============================================================
# PART 3: VISUAL TRANSFER  (9 individual panels)
# ============================================================
function generate_transfer_figures()
    n = 60
    domain = (-1, 1, -1, 1)
    grid = CartesianDiscreteModel(domain, (n, n))

    function analytical_field(x)
        return sin(π * x[1]) * sin(π * x[2])
    end

    function get_reference_field(geom)
        info = grid_info(geom.grid)
        pts = Gridap.Geometry.get_node_coordinates(geom.grid)
        u_data = vec(collect(map(analytical_field, pts)))
        return CartesianMeshField(u_data, info)
    end

    circle = EvolvingDomains.Geometric.Circle(VectorValue(0.0, 0.0), 0.6)
    pts = Gridap.Geometry.get_node_coordinates(grid)
    phi_values = vec(collect(map(x -> circle(x), pts)))
    geom = EvolvingDiscreteGeometry(phi_values, grid)

    active_indices = get_active_indices(geom)
    u_ref = get_reference_field(geom)

    # --- Exp 1: Grid → Mesh → Grid ---
    get_active_indices(geom, :current)
    cut_geo = geom.cache.cut
    Ω_act = Triangulation(cut_geo, ACTIVE)
    reffe = ReferenceFE(lagrangian, Float64, 1)
    V = TestFESpace(Ω_act, reffe)
    setup_transfer(geom, V)

    u_mesh = grid_to_mesh(geom, u_ref)
    u_back = mesh_to_grid(geom, u_mesh)

    err_transfer_data = zeros(length(u_ref.data))
    max_err_transfer = 0.0
    for idx in active_indices
        err = abs(u_back.data[idx] - u_ref.data[idx])
        err_transfer_data[idx] = err
        max_err_transfer = max(max_err_transfer, err)
    end
    u_err_transfer = CartesianMeshField(err_transfer_data, u_ref.grid)

    # --- Exp 2: Extension (Bulk → Void) ---
    masked_data = [phi <= 0 ? val : 0.0
                   for (phi, val) in zip(geom.levelset, u_ref.data)]
    u_masked = CartesianMeshField(vec(collect(masked_data)), u_ref.grid)
    u_ext = EvolvingDomains.Transfer.extend(geom, u_masked)
    u_err_ext_data = abs.(u_ext.data .- u_ref.data)
    u_err_ext = CartesianMeshField(u_err_ext_data, u_ref.grid)
    max_err_ext = maximum(u_err_ext_data)

    # --- Exp 3: Mesh → Grid → Mesh ---
    u_fem_init = interpolate(analytical_field, V)
    u_grid_intermediate = mesh_to_grid(geom, u_fem_init)
    u_mesh_final = grid_to_mesh(geom, u_grid_intermediate)
    u_diff_fem = interpolate(u_fem_init - u_mesh_final, V)
    u_err_mesh_data_on_grid = mesh_to_grid(geom, u_diff_fem).data
    u_err_mesh_grid = CartesianMeshField(abs.(u_err_mesh_data_on_grid), u_ref.grid)

    max_err_mesh = 0.0
    for idx in active_indices
        val = u_err_mesh_grid.data[idx]
        max_err_mesh = max(max_err_mesh, val)
    end

    u_fem_init_grid   = mesh_to_grid(geom, u_fem_init)
    u_mesh_final_grid = mesh_to_grid(geom, u_mesh_final)

    # --- Plotting ---
    experiments = [
        ("exp1_input",  u_ref,            "Exp 1: Grid → Mesh → Grid — Input",  nothing),
        ("exp1_result", u_back,           "Exp 1: Grid → Mesh → Grid — Result", nothing),
        ("exp1_error",  u_err_transfer,   "Exp 1: Grid → Mesh → Grid — Error",  max_err_transfer),
        ("exp2_input",  u_masked,         "Exp 2: Extension — Input",           nothing),
        ("exp2_result", u_ext,            "Exp 2: Extension — Result",            nothing),
        ("exp2_error",  u_err_ext,        "Exp 2: Extension — Error",             max_err_ext),
        ("exp3_input",  u_fem_init_grid,  "Exp 3: Mesh → Grid → Mesh — Input",  nothing),
        ("exp3_result", u_mesh_final_grid,"Exp 3: Mesh → Grid → Mesh — Result", nothing),
        ("exp3_error",  u_err_mesh_grid,  "Exp 3: Mesh → Grid → Mesh — Error",  max_err_mesh),
    ]

    # Fine grid shared by all transfer panels
    fine_n = 300
    bounds = (-1.0, 1.0, -1.0, 1.0)

    for (name, u_field, title, max_err) in experiments
        fig = Figure(size=(700, 600))
        ax = Axis(fig[1, 1], aspect=1, title=title,
                  xlabel="x", ylabel="y",
                  xgridvisible=false, ygridvisible=false)

        xs_fine, ys_fine, data_fine = interpolate_to_fine_grid(u_field, fine_n, fine_n, bounds)

        # mask with analytical circle (skip for extension result so the exterior
        # extension is visible)
        if !occursin("exp2", name)
            for j in 1:fine_n, i in 1:fine_n
                x = xs_fine[i]
                y = ys_fine[j]
                if sqrt(x^2 + y^2) > 0.6
                    data_fine[i, j] = NaN
                end
            end
        end

        if occursin("error", name)
            cr = (max_err !== nothing && max_err > 0) ? (0, max_err) : (0, 1)
            hm = heatmap!(ax, xs_fine, ys_fine, data_fine,
                          colormap=cgrad([:white, :gray50]), colorrange=cr)
        else
            hm = heatmap!(ax, xs_fine, ys_fine, data_fine,
                          colormap=cgrad([:white, :gray50]), colorrange=(-1, 1))
        end

        # smooth interface on the fine grid
        phi_fine = [ sqrt(x^2 + y^2) - 0.6  for x in xs_fine, y in ys_fine ]
        contour!(ax, xs_fine, ys_fine, phi_fine,
                 levels=[0.0], color=:black, linewidth=2)

        Colorbar(fig[1, 2], hm)

        save(joinpath(output_dir, "transfer_$(name).pdf"), fig)
    end
end


# ============================================================
# MAIN
# ============================================================
println("Generating publication figures …")
println("Output directory: $output_dir")

println("\n--- Conservative Transport ---")
generate_conservative_figures()

println("\n--- Dumbbell Temperature ---")
generate_dumbbell_figures()

println("\n--- Visual Transfer ---")
generate_transfer_figures()

println("\n✓ Done. Files saved to: $output_dir")
