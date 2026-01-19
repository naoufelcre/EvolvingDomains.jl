
# Stress Test: Uniform Grid Robustness
# =============================================================================
# Goal: Determine if the Uniform Grid EVER fails with "AssertionError: all_aggregated".
# If it never fails, then the Quadtree implementation is introducing the topological artifacts.

using Gridap
using Gridap.Geometry
using GridapEmbedded
using EvolvingDomains
using QuadtreeAgFEM
using LinearAlgebra
using GridapEmbedded.Interfaces: CUT

# === Geometry Definition (Unperturbed) ===
# We use the EXACT original geometry that caused the Quadtree L=6 crash
const R_geo = 0.5
const L_geo = 0.5 * R_geo
const center1 = Point(0.0, 0.0)
const center2 = center1 + VectorValue(-L_geo, L_geo)
function crescent_sdf(x)
    p = Point(x...)
    d1 = norm(p - center1) - R_geo
    d2 = norm(p - center2) - R_geo
    return max(d1, -d2)
end

function test_uniform_robustness()
    levels = 4:9 # Test a range of levels around the problematic L=6
    
    println("\n=== Uniform Grid Robustness Test ===")
    println("Geometry: Crescent Moon (Unperturbed)")
    println("Levels: $levels")
    
    for L in levels
        println("\nTesting Level L=$L...")
        
        # 1. Generate Uniform Model (Cartesian)
        n = 2^L
        domain = (-0.75, 0.75, -0.75, 0.75)
        partition = (n, n)
        bg_model = CartesianDiscreteModel(domain, partition)
        
        # 2. Cut Geometry
        # We need a sufficiently fine LSM for the geometry
        N_lsm = ceil(Int, 1.5 * n)
        lsm_model = CartesianDiscreteModel(domain, (N_lsm, N_lsm))
        geo = EvolvingDiscreteGeometry(bg_model, lsm_model, crescent_sdf)
        cut_geo = current_cut(geo)
        
        # 3. Aggregation (Standard)
        # We use the standard AggregateCutCellsByThreshold from GridapEmbedded
        # indirectly via RobustAggregation (or we can call standard directly)
        
        try
            # We use RobustAggregation because it wraps the standard logic but improves iteration.
            # However, the user asked if "Uniform grid failures" occur.
            # Standard GridapEmbedded usually asserts.
            
            strategy = RobustAggregation(0.5, 100)
            aggs = aggregate(strategy, cut_geo)
            
            # Check for orphans (RobustAgFEM doesn't crash, it warns)
            # But wait, RobustAgFEM is designed for Unstructured/Quadtree?
            # It works on any model.
            
            # Let's count orphans
            orphans = count(x -> x == 0, aggs)
            
            # Filter real orphans (CUT cells)
            # using GridapEmbedded.Interfaces: CUT # MOVED TO TOP
            cell_to_inoutcut = GridapEmbedded.AgFEM.compute_bgcell_to_inoutcut(bg_model, geo.ls_geo) # Assuming ls_geo fix
            
            real_orphans = 0
            for (c, root) in enumerate(aggs)
                if root == 0 && cell_to_inoutcut[c] == CUT
                    real_orphans += 1
                end
            end
            
            if real_orphans > 0
                println("  [FAIL] Level $L has $real_orphans orphans!")
            else
                println("  [PASS] Level $L: 0 orphans. Success.")
            end
            
        catch e
            println("  [CRASH] Level $L crashed with error: $e")
            # showstack(e)
        end
    end
end

# Fix for ls_geo access (from debugging session)
function EvolvingDomains.current_geometry(geom::EvolvingDiscreteGeometry)
    ϕ = current_levelset(geom)
    return EvolvingDomains._build_discrete_geometry(ϕ, geom.bg_model, geom.lsm_model)
end

# Helper to access the internal CSG.Geometry if needed, 
# but we can get classification from `cut_geo`.
# Actually `compute_bgcell_to_inoutcut` needs `(cut, geo)`. Here `cut` is `cut_geo`. `geo` is `cut_geo.geo`.

function test_robustness_fixed()
     levels = 4:9 
    
    println("\n=== Uniform Grid Robustness Test ===")
    println("Geometry: Crescent Moon (Unperturbed)")
    
    for L in levels
        print("Testing Level L=$L ... ")
        
        n = 2^L
        domain = (-0.75, 0.75, -0.75, 0.75)
        partition = (n, n)
        bg_model = CartesianDiscreteModel(domain, partition)
        
        N_lsm = ceil(Int, 1.5 * n)
        lsm_model = CartesianDiscreteModel(domain, (N_lsm, N_lsm))
        geo = EvolvingDiscreteGeometry(bg_model, lsm_model, crescent_sdf)
        cut_geo = current_cut(geo)
        
        strategy = RobustAggregation(0.5, 100)
        aggs = aggregate(strategy, cut_geo)
        
        # using GridapEmbedded.Interfaces: CUT # MOVED TO TOP
        cell_to_inoutcut = GridapEmbedded.AgFEM.compute_bgcell_to_inoutcut(bg_model, cut_geo.geo)
        
        real_orphans = 0
        for (c, root) in enumerate(aggs)
            if root == 0 && cell_to_inoutcut[c] == CUT
                real_orphans += 1
            end
        end
        
        if real_orphans > 0
            println("[FAIL] $real_orphans orphans.")
        else
            println("[PASS] 0 orphans.")
        end
    end
end

test_robustness_fixed()
