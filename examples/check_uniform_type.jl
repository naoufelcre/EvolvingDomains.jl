
# Check Uniform Grid Element Type
using Gridap
using GridapEmbedded
using Gridap.ReferenceFEs

function check_types()
    println("=== Checking Element Types ===")
    
    # 1. Cartesian Model (Uniform)
    domain = (0,1,0,1)
    partition = (4,4)
    model = CartesianDiscreteModel(domain, partition)
    
    grid = get_grid(model)
    reffes = get_reffes(grid)
    println("Uniform Grid Reffes: ", reffes)
    
    # Check if it's QUAD or TRI
    is_quad = any(r -> get_polytope(r) == QUAD, reffes)
    is_tri = any(r -> get_polytope(r) == TRI, reffes)
    
    println("  Has Quads? $is_quad")
    println("  Has Tris? $is_tri")
    
    # 2. Check Input to Coarsening
    # The input to `bottom_up_coarsening!` is a `QuadMesh`.
    # A `QuadMesh` is a tree structure where leaves are squares.
    # It is effectively a grid of squares (Quads) before paving.
    # The "Paving" step converts it to a Gridap Model.
    
    println("\n=== Input to Coarsening ===")
    println("The input to coarsening is `QuadMesh`, which is a tree of squares.")
    println("Each node represents a spatial box (Square/Quad).")
end

check_types()
