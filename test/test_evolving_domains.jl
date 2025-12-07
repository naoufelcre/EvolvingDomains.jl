# =============================================================================
# Test: EvolvingDomains with LevelSetMethods.jl Backend
# =============================================================================
# This test demonstrates the full workflow:
# 1. Create a background Cartesian model
# 2. Define an evolving level set (circle advected by rotation velocity)
# 3. Advance the geometry
# 4. Use the resulting DiscreteGeometry with cut()
# 5. Visualize with EmbeddedViz (multiple GIF outputs)
#
# Includes performance timing for all operations.

using Gridap
using GridapEmbedded
using LinearAlgebra
using Printf
using Statistics
using Dates
using JSON

# Add the local modules
push!(LOAD_PATH, joinpath(@__DIR__, "..", "src"))
push!(LOAD_PATH, joinpath(@__DIR__, "..", "..", "EmbeddedViz", "src"))

using EvolvingDomains
using EmbeddedViz

# Try to load LevelSetMethods (optional)
lsm_available = try
    using LevelSetMethods
    true
catch
    false
end

println("=" ^ 70)
println("       CIRCLE EVOLUTION TEST - PERFORMANCE REPORT")
println("=" ^ 70)
println("LevelSetMethods.jl available: ", lsm_available)
println()

# Performance tracking
timings = Dict{String, Float64}()
step_times = Float64[]
cut_times = Float64[]

# Initialize visualization
setup_viz!()

if lsm_available
    # =========================================================================
    # 1. SETUP
    # =========================================================================
    println("--- 1. Setup ---")

    t_setup = @elapsed begin
        # Background model
        domain = (-1.0, 1.0, -1.0, 1.0)
        partition = (50, 50)
        model = CartesianDiscreteModel(domain, partition)
    end
    timings["Model Creation"] = t_setup
    println("  Background model: $(partition[1])×$(partition[2]) cells")
    println("  Domain: $domain")
    @printf("  Time: %.3f ms\n", t_setup * 1000)

    # Initial level set (circle at center)
    R = 0.2
    center = (0.5, 0.5)
    ϕ0(x) = -(R - norm(x .- center))

    # Velocity (rigid rotation around origin)
    u(x) = (-x[2], x[1])

    # Create evolver
    t_evolver = @elapsed begin
        evolver = LevelSetMethodsEvolver(;
            bg_model = model,
            initial_ls = ϕ0,
            velocity = u,
            integrator = :RK3,
            reinit_freq = nothing,
            bc = :Neumann
        )
    end
    timings["Evolver Creation"] = t_evolver
    @printf("  Evolver creation: %.3f ms\n", t_evolver * 1000)

    # Create evolving geometry
    t_eg = @elapsed begin
        eg = EvolvingDiscreteGeometry(evolver, model)
    end
    timings["EvolvingGeometry Creation"] = t_eg
    @printf("  EvolvingGeometry creation: %.3f ms\n", t_eg * 1000)

    # =========================================================================
    # 2. INITIAL CUT
    # =========================================================================
    println("\n--- 2. Initial Cut ---")

    t_cut_init = @elapsed begin
        geo0 = current_geometry(eg)
        cut0 = cut(model, geo0)
        trian_in0 = Triangulation(cut0, PHYSICAL_IN)
    end
    timings["Initial Cut"] = t_cut_init
    println("  Interior cells: ", num_cells(trian_in0))
    @printf("  Cut time: %.3f ms\n", t_cut_init * 1000)

    # =========================================================================
    # 3. EVOLUTION LOOP
    # =========================================================================
    println("\n--- 3. Evolution Loop ---")
    Δt = 0.1
    nsteps = 100

    # Create grid info for visualization
    grid = GridInfo((-1.0, -1.0), (1.0, 1.0), (50, 50))

    # Create multi-type recorder
    rec = Recorder()

    # Record initial state
    ϕ = current_levelset(eg)
    t = EvolvingDomains.current_time(eg)

    t_record_init = @elapsed begin
        record!(rec, :levelset, ϕ, grid; t=t, label="initial")
        record!(rec, :levelset3d, ϕ, grid; t=t, label="initial")
        record!(rec, :grid, ϕ, grid; t=t, label="initial")
    end
    timings["Initial Recording"] = t_record_init

    # Evolution loop with timing
    total_advance_time = 0.0
    total_record_time = 0.0

    println("  Running $nsteps steps with Δt = $Δt...")

    for step in 1:nsteps
        # Time the advance
        t_adv = @elapsed advance!(eg, Δt)
        global total_advance_time += t_adv
        push!(step_times, t_adv)

        local t_current = EvolvingDomains.current_time(eg)
        local ϕ_current = current_levelset(eg)

        # Time the cut
        local cut_geo
        t_cut = @elapsed begin
            geo = current_geometry(eg)
            cut_geo = cut(model, geo)
        end
        push!(cut_times, t_cut)

        # Time the recording
        t_rec = @elapsed begin
            record!(rec, :levelset, ϕ_current, grid; t=t_current, label="step$step")
            record!(rec, :levelset3d, ϕ_current, grid; t=t_current, label="step$step")
            record!(rec, :grid, ϕ_current, grid; t=t_current, label="step$step")
            record!(rec, :cutmesh, cut_geo, grid; t=t_current, label="step$step")
        end
        global total_record_time += t_rec

        # Progress indicator every 20 steps
        if step % 20 == 0
            @printf("    Step %3d/%-3d: t = %.2f, advance = %.2f ms, cut = %.2f ms, record = %.2f ms\n",
                    step, nsteps, t_current, t_adv * 1000, t_cut * 1000, t_rec * 1000)
        end
    end

    timings["Total Advance"] = total_advance_time
    timings["Total Recording"] = total_record_time
    timings["Avg Advance/Step"] = total_advance_time / nsteps
    timings["Avg Record/Step"] = total_record_time / nsteps

    # =========================================================================
    # 4. SAVE ANIMATIONS
    # =========================================================================
    println("\n--- 4. Saving Animations ---")

    t_save = @elapsed save_all(rec, "circle_evolution_"; fps=5)
    timings["Save Animations"] = t_save
    @printf("  Save time: %.3f ms\n", t_save * 1000)

    # =========================================================================
    # 5. FINAL CUT
    # =========================================================================
    println("\n--- 5. Final Cut ---")

    t_cut_final = @elapsed begin
        geo1 = current_geometry(eg)
        cut1 = cut(model, geo1)
        trian_in1 = Triangulation(cut1, PHYSICAL_IN)
    end
    timings["Final Cut"] = t_cut_final
    println("  Interior cells: ", num_cells(trian_in1))
    @printf("  Cut time: %.3f ms\n", t_cut_final * 1000)

    # =========================================================================
    # PERFORMANCE REPORT
    # =========================================================================
    println("\n" * "=" ^ 70)
    println("                    PERFORMANCE SUMMARY")
    println("=" ^ 70)
    println()

    # Grid info
    n_nodes = (partition[1] + 1) * (partition[2] + 1)
    println("Grid: $(partition[1])×$(partition[2]) cells, $n_nodes nodes")
    println("Steps: $nsteps, Δt = $Δt, Final time = $(nsteps * Δt)")
    println()

    # Timing breakdown
    println("┌─────────────────────────────────┬────────────────┐")
    println("│ Operation                       │ Time (ms)      │")
    println("├─────────────────────────────────┼────────────────┤")

    timings["Total Cutting"] = sum(cut_times)
    timings["Avg Cut/Step"] = sum(cut_times) / nsteps

    ordered_keys = ["Model Creation", "Evolver Creation", "EvolvingGeometry Creation",
                    "Initial Cut", "Total Advance", "Total Cutting", "Total Recording",
                    "Save Animations", "Final Cut"]

    for key in ordered_keys
        if haskey(timings, key)
            @printf("│ %-31s │ %12.3f   │\n", key, timings[key] * 1000)
        end
    end

    println("├─────────────────────────────────┼────────────────┤")
    total_time = sum(values(timings)) - timings["Avg Advance/Step"] - timings["Avg Record/Step"]
    @printf("│ %-31s │ %12.3f   │\n", "TOTAL", total_time * 1000)
    println("└─────────────────────────────────┴────────────────┘")
    println()

    # Per-step statistics
    println("Per-Step Statistics:")
    @printf("  Average advance time:  %.3f ms/step\n", timings["Avg Advance/Step"] * 1000)
    @printf("  Average record time:   %.3f ms/step\n", timings["Avg Record/Step"] * 1000)
    @printf("  Min advance time:      %.3f ms\n", minimum(step_times) * 1000)
    @printf("  Max advance time:      %.3f ms\n", maximum(step_times) * 1000)
    @printf("  Std dev advance time:  %.3f ms\n", std(step_times) * 1000)
    @printf("  Average cut time:      %.3f ms/step\n", timings["Avg Cut/Step"] * 1000)
    @printf("  Min cut time:          %.3f ms\n", minimum(cut_times) * 1000)
    @printf("  Max cut time:          %.3f ms\n", maximum(cut_times) * 1000)
    @printf("  Std dev cut time:      %.3f ms\n", std(cut_times) * 1000)
    println()

    # Performance metrics
    time_per_node = timings["Avg Advance/Step"] / n_nodes * 1e6  # microseconds
    @printf("Performance: %.3f μs/node/step\n", time_per_node)
    println()

    println("Generated files:")
    println("  - circle_evolution_levelset.gif (2D contour)")
    println("  - circle_evolution_levelset3d.gif (3D surface)")
    println("  - circle_evolution_grid.gif (grid with IN/OUT/CUT cells)")
    println("  - circle_evolution_cutmesh.gif (physical domain triangulation)")

    # =========================================================================
    # SAVE METADATA
    # =========================================================================
    metadata = Dict(
        "test_name" => "circle_evolution",
        "timestamp" => Dates.format(now(), "yyyy-mm-dd_HH:MM:SS"),
        "grid" => Dict(
            "domain" => domain,
            "partition" => partition,
            "n_cells" => partition[1] * partition[2],
            "n_nodes" => n_nodes
        ),
        "simulation" => Dict(
            "nsteps" => nsteps,
            "dt" => Δt,
            "final_time" => nsteps * Δt,
            "integrator" => "RK3",
            "spatial_scheme" => "WENO5",
            "reinit_freq" => nothing
        ),
        "geometry" => Dict(
            "type" => "circle",
            "center" => [center...],
            "radius" => R
        ),
        "timings_ms" => Dict(k => v * 1000 for (k, v) in timings),
        "step_times_ms" => step_times .* 1000,
        "cut_times_ms" => cut_times .* 1000,
        "statistics" => Dict(
            "avg_advance_ms" => timings["Avg Advance/Step"] * 1000,
            "avg_record_ms" => timings["Avg Record/Step"] * 1000,
            "avg_cut_ms" => timings["Avg Cut/Step"] * 1000,
            "min_advance_ms" => minimum(step_times) * 1000,
            "max_advance_ms" => maximum(step_times) * 1000,
            "std_advance_ms" => std(step_times) * 1000,
            "min_cut_ms" => minimum(cut_times) * 1000,
            "max_cut_ms" => maximum(cut_times) * 1000,
            "std_cut_ms" => std(cut_times) * 1000,
            "us_per_node_per_step" => time_per_node,
            "total_time_ms" => total_time * 1000
        ),
        "results" => Dict(
            "initial_interior_cells" => num_cells(trian_in0),
            "final_interior_cells" => num_cells(trian_in1)
        )
    )

    metadata_file = "circle_perf_" * Dates.format(now(), "yyyymmdd_HHMMSS") * ".json"
    open(metadata_file, "w") do f
        JSON.print(f, metadata, 2)
    end
    println()
    println("Performance metadata saved to: $metadata_file")
    println()
    println("=" ^ 70)
    println("                    TEST PASSED ✓")
    println("=" ^ 70)

else
    println("\n--- Fallback Test (No LevelSetMethods) ---")

    # Test the abstract interface exists
    println("AbstractLevelSetEvolver defined: ", isdefined(EvolvingDomains, :AbstractLevelSetEvolver))
    println("EvolvingDiscreteGeometry defined: ", isdefined(EvolvingDomains, :EvolvingDiscreteGeometry))

    # Test utilities
    model = CartesianDiscreteModel((0,1,0,1), (5,5))
    coords = EvolvingDomains.get_node_coords_as_vector(model)
    println("Node coords extracted: ", length(coords), " nodes")

    desc = EvolvingDomains.cartesian_descriptor(model)
    println("Descriptor: origin=$(desc[1]), corner=$(desc[2]), partition=$(desc[3])")

    println("\n--- Fallback Test Passed! ---")
    println("Install LevelSetMethods.jl to run full test:")
    println("  ] add LevelSetMethods")
end
