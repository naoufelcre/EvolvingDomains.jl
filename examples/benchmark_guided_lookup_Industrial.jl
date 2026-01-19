using EvolvingDomains
using Gridap
using Gridap.Geometry
using Gridap.ReferenceFEs
using Gridap.FESpaces
using Gridap.Arrays
using Gridap.CellData
using Gridap.Fields
using LinearAlgebra
using Random
using Printf

# =============================================================================
# Benchmark: Guided Velocity Source (Scalar Refactor Verification)
# =============================================================================
# Verification of the fix for the "Piola Mapping Artifact" (Vector Element Scaling).
# We now use two scalar fields (u_x, u_y) to represent velocity.

function run_benchmark()
    println("═" ^ 60)
    println("  BENCHMARK: Guided Velocity Source (Scalar Components)")
    println("═" ^ 60)

    # 1. Generate Quadtree Mesh (L=5, Gradient Sizing)
    println("\n[1] Generating Quadtree Mesh (L=5)...")
    L = 5
    h_fine = 1.0 / (2^L)
    
    # Sizing function for some variation
    function sizing(x)
        d = norm(x .- 0.5)
        return max(h_fine, 0.4 * d)
    end
    
    qmesh = generate_fine_mesh(1.0, L)
    bottom_up_coarsening!(qmesh, L, sizing)
    balance!(qmesh)
    elements = pave_mesh(qmesh)
    model, leaf_map = quadtree_to_discrete_model(elements)
    
    println("    Mesh Stats: $(num_cells(model)) cells, $(num_nodes(model)) nodes.")

    # 2. Define Analytical Velocity Field (Taylor-Green Vortex like)
    println("\n[2] Defining Velocity Field (Scalar Components)...")
    
    # Vectors for Truth
    function analytical_velocity(x)
        u = sin(2*π*x[1]) * cos(2*π*x[2])
        v = -cos(2*π*x[1]) * sin(2*π*x[2])
        return VectorValue(u, v)
    end
    
    # Scalar components for interpolation
    function analytical_u(x)
        return sin(2*π*x[1]) * cos(2*π*x[2])
    end
    
    function analytical_v(x)
        return -cos(2*π*x[1]) * sin(2*π*x[2])
    end
    
    # Project onto Scalar FE Space (Float64)
    # This is the KEY FIX: Using Scalar Elements avoids Piola mapping issues.
    reffe = ReferenceFE(lagrangian, Float64, 1)
    V_scal = FESpace(model, reffe)
    
    println("    Interpolating u_x...")
    u_h_x = interpolate(analytical_u, V_scal)
    println("    Interpolating u_y...")
    u_h_y = interpolate(analytical_v, V_scal)
    
    # 3. Setup Guided Velocity Source
    # Passing Tuple of Scalar Fields
    println("\n[3] Building GuidedVelocitySource...")
    vel_source = GuidedVelocitySource(u_h_x, u_h_y, qmesh.root, leaf_map)
    
    # 4. Verification Loop
    println("\n[4] Running Verification (N=1000 Samples)...")
    
    rng = MersenneTwister(1234)
    n_samples = 1000
    
    max_error = 0.0
    sum_sq_error = 0.0
    
    # Distribution buckets
    err_buckets = zeros(Int, 5) # <1e-6, <1e-4, <1e-2, <0.1, >=0.1
    
    for i in 1:n_samples
        # Random Point in [0,1]
        x = Point(rand(rng), rand(rng))
        
        # A. Analytical Truth
        v_true = analytical_velocity(x)
        
        # B. Guided Lookup
        # Returns Tuple (u, v)
        v_guided_tuple = sample_velocity(vel_source, (x[1], x[2]), 0.0)
        v_guided = VectorValue(v_guided_tuple[1], v_guided_tuple[2])
        
        # Metric
        diff = iszero(norm(v_guided)) ? norm(v_true) : norm(v_guided - v_true)
        
        push_bucket!(err_buckets, diff)
        
        max_error = max(max_error, diff)
        sum_sq_error += diff^2
        
        # Early warning for gross failure
        if diff > 0.1 && i <= 10
             println("    [High Error] At $x: True=$v_true, Guided=$v_guided, Diff=$diff")
        end
    end
    
    l2_error = sqrt(sum_sq_error / n_samples)
    
    println("\n" * "-"^40)
    println("  RESULTS")
    println("-"^40)
    @printf("  Max Error: %.2e\n", max_error)
    @printf("  L2 Error:  %.2e\n", l2_error)
    println("  Error Distribution:")
    println("    < 1e-6: $(err_buckets[1])")
    println("    < 1e-4: $(err_buckets[2])")
    println("    < 1e-2: $(err_buckets[3])")
    println("    < 0.1:  $(err_buckets[4])")
    println("    >= 0.1: $(err_buckets[5])")
    
    if max_error > 0.1
        println("\n  [FAIL] High Error Persists.")
    elseif max_error < 1e-2
        println("\n  [PASS] SCALAR REFACTOR SUCCESSFUL! (Max Error < 1e-2)")
    else
        println("\n  [WARN] Error is moderate. (1e-2 < Err < 0.1)")
    end
    println("═" ^ 60)
end

function push_bucket!(buckets, err)
    if err < 1e-6
        buckets[1] += 1
    elseif err < 1e-4
        buckets[2] += 1
    elseif err < 1e-2
        buckets[3] += 1
    elseif err < 0.1
        buckets[4] += 1
    else
        buckets[5] += 1
    end
end

if abspath(PROGRAM_FILE) == @__FILE__
    run_benchmark()
end
