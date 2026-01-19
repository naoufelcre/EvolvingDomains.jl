
# Debug Script: Inspect Orphans
# =============================================================================
push!(LOAD_PATH, joinpath(@__DIR__, "../src"))
using Gridap
using Gridap.Geometry
using GridapEmbedded
using EvolvingDomains
using QuadtreeAgFEM
using LinearAlgebra
using GridapEmbedded.Interfaces: CUT, IN, OUT

# Copied Geometry Setup (Same as Benchmark)
const EPS_SHIFT = 1.0e-5
const OFFSET = VectorValue(EPS_SHIFT, EPS_SHIFT)
const R_geo = 0.5
const L_geo = 0.5 * R_geo
const center1 = Point(0.0, 0.0) + OFFSET
const center2 = center1 + VectorValue(-L_geo, L_geo)
function crescent_sdf(x)
    p = Point(x...)
    d1 = norm(p - center1) - R_geo
    d2 = norm(p - center2) - R_geo
    return max(d1, -d2)
end

function inspect_orphans()
    L = 6 # The problematic level
    println("Inspecting Orphans at Level $L...")
    
    # 1. Generate Model
    println("  Generating Model...")
    h_fine = 1.0 / (2^L)
    origin = Point(-0.75, -0.75)
    width = 1.5
    map_to_physical(x_unit) = origin + width * VectorValue(x_unit[1], x_unit[2])
    function sizing(x_unit)
        x_phys = map_to_physical(x_unit)
        return max(h_fine, 0.5*abs(crescent_sdf(x_phys)))
    end
    qmesh = generate_fine_mesh(1.0, L)
    bottom_up_coarsening!(qmesh, L, sizing)
    balance!(qmesh)
    elements = pave_mesh(qmesh)
    unit_model, _ = quadtree_to_discrete_model(elements)
    grid_unit = get_grid(unit_model)
    nodes = Gridap.Geometry.get_node_coordinates(grid_unit)
    new_nodes = [map_to_physical(n) for n in nodes]
    grid_phys = UnstructuredGrid(new_nodes, get_cell_node_ids(grid_unit), get_reffes(grid_unit), get_cell_type(grid_unit), Gridap.Geometry.NonOriented())
    bg_model = UnstructuredDiscreteModel(grid_phys)
    
    # 2. Cut Geometry
    N_lsm = ceil(Int, 1.5 * (2^(L+1)))
    lsm_model = CartesianDiscreteModel((-0.75, 0.75, -0.75, 0.75), (N_lsm, N_lsm))
    geo = EvolvingDiscreteGeometry(bg_model, lsm_model, crescent_sdf)
    cut_geo = current_cut(geo)
    
    # 3. Run Aggregation to get Orphans
    strategy = RobustAggregation(0.5, 100)
    aggs = aggregate(strategy, cut_geo)
    
    # 4. Analyze Orphans
    # aggs[i] == 0 means orphan
    orphans = findall(x -> x == 0, aggs)
    println("\nFound $(length(orphans)) orphans.")
    
    if isempty(orphans)
        println("No orphans found. Exiting.")
        return
    end
    
    # Get properties
    # FIX: `geo` does not have `ls_geo`. It has `ls_grid` or we access via `current_geometry(geo)` or `geo.ls_grid`?
    # Inspecting EvolvingDiscreteGeometry struct:
    # struct EvolvingDiscreteGeometry{Gm, Lm, F}
    #    bg_model::Gm
    #    ls_grid::Lm ...
    
    # We need the geometry used for cutting. `current_cut(geo)` returns an EmbeddedDiscretization which has `ls`?
    # Actually `AgFEM.compute_bgcell_to_inoutcut` needs `(cut, geo)`.
    # Based on `RobustAgFEM.jl`: `GridapEmbedded.AgFEM.compute_bgcell_to_inoutcut(cut,geo)`
    # Here `geo` is the `CSG.Geometry`.
    
    # In `EvolvingGeometry.jl`, `current_cut` returns `cut.embedded_disk`.
    # `EvolvingDiscreteGeometry` stores the LS geometry in `geo.geom`.
    # Let's check the file content first. 
    # Attempting to access `geo.ls` or similar.
    # In the `solve_single_iteration`, we saw `geo = EvolvingDiscreteGeometry(...)`.
    # `current_cut(geo)` returns the cut object.
    
    # Wait, `AgFEM.aggregate` took `cut_geo` and `cut_geo.geo`?
    # `cut_geo` is the `EmbeddedDiscretization`.
    # `cut_geo.geo` is the `CSG.Geometry`.
    
    # Let's get the classification from the cut object directly if possible, or recompute.
    cell_to_inoutcut = GridapEmbedded.AgFEM.compute_bgcell_to_inoutcut(cut_geo, cut_geo.geo)
    
    # Filter orphans: We only care about ORPHANS that are CUT.
    # aggs[c] == 0 means "No Root".
    # For OUT cells, this is normal.
    # For CUT cells, this is an error.
    
    real_orphans = Int[]
    for c in orphans
        if cell_to_inoutcut[c] == CUT
            push!(real_orphans, c)
        end
    end
    
    println("\nFound $(length(real_orphans)) CUT orphans (Expected ~12).")
    
    orphans = real_orphans # Replace for loop below
    println("Cell ID | Classification | Is CUT? | Neighbors (Class)")
    
    topo = Gridap.Geometry.get_grid_topology(bg_model)
    D = Gridap.Geometry.num_cell_dims(bg_model)
    cell_to_faces = Gridap.Geometry.get_faces(topo,D,D-1)
    face_to_cells = Gridap.Geometry.get_faces(topo,D-1,D)
    
    c1 = Gridap.Arrays.array_cache(cell_to_faces)
    c2 = Gridap.Arrays.array_cache(face_to_cells)
    
    for c in orphans
        class_code = cell_to_inoutcut[c]
        class_str = (class_code == CUT ? "CUT" : (class_code == IN ? "IN" : "OUT"))
        
        # Get neighbors
        faces = Gridap.Arrays.getindex!(c1, cell_to_faces, c)
        neigh_str = ""
        for f in faces
            neighs = Gridap.Arrays.getindex!(c2, face_to_cells, f)
            for n in neighs
                if n != c
                    n_class = cell_to_inoutcut[n]
                    n_sym = (n_class == CUT ? "C" : (n_class == IN ? "I" : "O"))
                    neigh_str *= " $n($n_sym)"
                end
            end
        end
        
        println("$c | $class_str | $(class_code == CUT) | $neigh_str")
    end
end

inspect_orphans()
