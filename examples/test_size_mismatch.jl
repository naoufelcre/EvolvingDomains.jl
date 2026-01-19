
# Compare Quadtree vs Uniform Resolution
# =============================================================================

using Gridap
using Gridap.Geometry
using QuadtreeAgFEM
using LinearAlgebra
using Statistics

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

function get_mesh_stats(model)
    trian = Triangulation(model)
    meas = get_cell_measure(trian)
    vols = collect(meas)
    
    # Approx diameter ~ sqrt(area)
    diams = sqrt.(vols)
    
    return length(vols), minimum(diams), mean(diams), maximum(diams)
end

function compare_resolution()
    L = 6
    println("Comparing Resolution at Level L=$L")
    
    # 1. Quadtree
    h_fine = 1.0 / (2^L)
    origin = Point(-0.75, -0.75)
    width = 1.5
    map_to_physical(x_unit) = origin + width * VectorValue(x_unit[1], x_unit[2])
    function sizing(x_unit)
        x_phys = map_to_physical(x_unit)
        # Sizing from benchmark
        val = max(h_fine, 0.5*abs(crescent_sdf(x_phys)))
        return val
    end
    
    println("\n[Quadtree Generation]")
    qmesh = generate_fine_mesh(1.0, L)
    bottom_up_coarsening!(qmesh, L, sizing)
    balance!(qmesh)
    elements = pave_mesh(qmesh)
    unit_model, _ = quadtree_to_discrete_model(elements)
    
    # Map to physical for fair comparison
    grid_unit = get_grid(unit_model)
    nodes = Gridap.Geometry.get_node_coordinates(grid_unit)
    new_nodes = [map_to_physical(n) for n in nodes]
    grid_phys = UnstructuredGrid(new_nodes, get_cell_node_ids(grid_unit), get_reffes(grid_unit), get_cell_type(grid_unit), Gridap.Geometry.NonOriented())
    qt_model = UnstructuredDiscreteModel(grid_phys)
    
    n_qt, min_d_qt, avg_d_qt, max_d_qt = get_mesh_stats(qt_model)
    
    # 2. Uniform
    println("\n[Uniform Generation]")
    n = 2^L
    domain = (-0.75, 0.75, -0.75, 0.75)
    partition = (n, n)
    uni_model = CartesianDiscreteModel(domain, partition)
    
    n_uni, min_d_uni, avg_d_uni, max_d_uni = get_mesh_stats(uni_model)
    
    # 3. Report
    println("\n=== Resolution Comparison (Physical Domain width 1.5) ===")
    println("Target h_fine (L=6) ~ $(1.5 * h_fine)")
    
    println("Quadtree:")
    println("  Cells: $n_qt")
    println("  Min h: $min_d_qt")
    println("  Avg h: $avg_d_qt")
    
    println("Uniform:")
    println("  Cells: $n_uni")
    println("  Min h: $min_d_uni")
    println("  Avg h: $avg_d_uni")
    
    ratio = n_uni / n_qt
    println("\nRatio (Uniform/Quadtree Cells): $ratio")
end

compare_resolution()
