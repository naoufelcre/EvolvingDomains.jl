
# Benchmark Script: Poisson Problem on "Crescent Moon" (PERTURBED)
# Testing if L=6 singularity is resolved by breaking grid alignment.
# =============================================================================

push!(LOAD_PATH, joinpath(@__DIR__, "../src"))

using Gridap
using Gridap.Geometry
using GridapEmbedded
using EvolvingDomains
using QuadtreeAgFEM # Explicitly needed now!
using LinearAlgebra
using Printf
using GridapEmbedded.Interfaces: CUT, IN, OUT

# =============================================================================
# 1. Geometry with PERTURBATION
# =============================================================================

# We shift the origin by a tiny amount
const EPS_SHIFT = 1.0e-5
const OFFSET = VectorValue(EPS_SHIFT, EPS_SHIFT)

const R_geo = 0.5
const L_geo = 0.5 * R_geo
const center1 = Point(0.0, 0.0) + OFFSET
const center2 = center1 + VectorValue(-L_geo, L_geo)

println("[Test] Running with Origin Offset: $OFFSET")

function crescent_sdf(x)
    p = Point(x...)
    d1 = norm(p - center1) - R_geo
    d2 = norm(p - center2) - R_geo
    return max(d1, -d2)
end

# Analytical Solution: u(x) = x - y
u_func(x) = x[1] - x[2]
∇u_func(x) = VectorValue(1.0, -1.0)
f_func(x) = 0.0

# =============================================================================
# 2. Mesh Generation (Copied)
# =============================================================================

function generate_quadtree_model(L_refine::Int)
    h_fine = 1.0 / (2^L_refine)
    
    # Mapping [0,1] -> [-0.75, 0.75] (Length 1.5)
    origin = Point(-0.75, -0.75)
    width = 1.5
    map_to_physical(x_unit) = origin + width * VectorValue(x_unit[1], x_unit[2])
    
    function sizing_in_unit(x_unit)
        x_phys = map_to_physical(x_unit)
        dist = abs(crescent_sdf(x_phys))
        return max(h_fine, 0.5*dist)
    end

    qmesh = generate_fine_mesh(1.0, L_refine)
    bottom_up_coarsening!(qmesh, L_refine, sizing_in_unit)
    balance!(qmesh)
    elements = pave_mesh(qmesh)
    unit_model, _ = quadtree_to_discrete_model(elements)
    
    grid_unit = get_grid(unit_model)
    nodes = Gridap.Geometry.get_node_coordinates(grid_unit)
    new_nodes = [map_to_physical(n) for n in nodes]
    
    grid_phys = UnstructuredGrid(new_nodes, 
                                 get_cell_node_ids(grid_unit), 
                                 get_reffes(grid_unit), 
                                 get_cell_type(grid_unit), 
                                 Gridap.Geometry.NonOriented())
    return UnstructuredDiscreteModel(grid_phys)
end

# =============================================================================
# 3. Solver Workflow (Copied)
# =============================================================================

function solve_single_iteration(bg_model, L_refine)
    N_lsm = ceil(Int, 1.5 * (2^(L_refine+1)))
    domain_lsm = (-0.75, 0.75, -0.75, 0.75)
    lsm_model = CartesianDiscreteModel(domain_lsm, (N_lsm, N_lsm))
    
    try
        geo = EvolvingDiscreteGeometry(bg_model, lsm_model, crescent_sdf)
        cut_geo = current_cut(geo)
        
        # This is where it crashes usually (AggregateCutCellsByThreshold)
        # But here we use RobustAggregation
        strategy = RobustAggregation(0.5, 100)
        aggs = aggregate(strategy, cut_geo)
        
        # --- Visualization of Aggregates ---
        # RobustAgFEM returns the `cell_to_root` map for ALL background cells.
        # We visualize the full background mesh to match dimensions.
        
        trian_bg = Triangulation(bg_model)
        
        # --- Compute Classification (IN/OUT/CUT) ---
        # We need to visualize why these cells are orphans.
        # classification: IN, OUT, CUT
        cell_to_inoutcut = GridapEmbedded.AgFEM.compute_bgcell_to_inoutcut(bg_model, cut_geo.geo)
        
        # Map to Integer for Paraview: IN=1, OUT=-1, CUT=0
        class_field = map(cell_to_inoutcut) do c
            if c == IN; return 1; end
            if c == OUT; return -1; end
            if c == CUT; return 0; end
            return -99 # Should not happen
        end
        
        # Create output directory
        !isdir("output") && mkdir("output")
        
        # Write VTK with aggregation data
        # 'aggs' is a Vector{Int} mapping cell_id -> root_id
        writevtk(trian_bg, "output/debug_aggregation_L$L_refine", 
                 cellfields=["RootID" => aggs, 
                             "IsOrphan" => map(x -> x==0 ? 1.0 : 0.0, aggs),
                             "IsRoot" => map((r, c) -> r==c ? 1.0 : 0.0, aggs, 1:length(aggs)),
                             "Classification" => class_field])
        
        println("    [Viz] Debug VTK written to output/debug_aggregation_L$L_refine.vtu")
        # -----------------------------------
        
        println("    [Success] Aggregation complete for L=$L_refine")
        return true
    catch e
        println("    [Fail] Solver failed for L=$L_refine: ", e)
        return false
    end
end

# =============================================================================
# 4. Main Control
# =============================================================================

function run_test()
    # Only test the problematic level L=6
    L = 6
    println("Testing Level $L with perturbation...")
    
    model_q = generate_quadtree_model(L)
    success = solve_single_iteration(model_q, L)
    
    if success
        println("\n>>> TEST PASSED: Perturbation resolved the failure. <<<")
    else
        println("\n>>> TEST FAILED: Perturbation did not help. <<<")
    end
end

run_test()
