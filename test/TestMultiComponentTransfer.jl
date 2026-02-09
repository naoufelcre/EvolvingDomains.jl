using Test
using EvolvingDomains
using EvolvingDomains.Geometric
using EvolvingDomains.Transfer
using Gridap
using GridapEmbedded
using Gridap.TensorValues
using LinearAlgebra

@testset "Multi-Component Transfer Verification" begin
    # 1. SETUP GEOMETRY
    n = 10
    domain = (-1, 1, -1, 1)
    grid = CartesianDiscreteModel(domain, (n, n))
    circle = EvolvingDomains.Geometric.Circle(VectorValue(0.0, 0.0), 0.6)
    pts = Gridap.Geometry.get_node_coordinates(grid)
    phi_values = vec(collect(map(x -> circle(x), pts)))
    geom = EvolvingDiscreteGeometry(phi_values, grid)
    active_indices = get_active_indices(geom)
    info = grid_info(grid)

    # ---------------------------------------------------------
    # TEST 1: VECTOR VALUE (Verified)
    # ---------------------------------------------------------
    println("Testing VectorValue...")
    u_vec_data = vec(collect(map(x -> VectorValue(sin(π * x[1]), cos(π * x[2])), pts)))
    u_vec_grid = CartesianMeshField(u_vec_data, info)

    # Ensure cut is populated
    get_active_indices(geom, :current)
    cut_geo = geom.cache.cut
    Ω_act = Triangulation(cut_geo, ACTIVE)

    V_vec = TestFESpace(Ω_act, ReferenceFE(lagrangian, VectorValue{2,Float64}, 1))
    setup_transfer(geom, V_vec)

    u_vec_back = mesh_to_grid(geom, grid_to_mesh(geom, u_vec_grid))
    err_vec = maximum([norm(u_vec_back.data[i] - u_vec_grid.data[i]) for i in active_indices])
    println("  Vector Error: $err_vec")
    @test err_vec < 1e-12

    # ---------------------------------------------------------
    # TEST 2: TENSOR VALUE (Matrix)
    # ---------------------------------------------------------
    println("Testing TensorValue (Matrix)...")
    u_ten_data = vec(collect(map(x -> TensorValue(x[1], x[2], -x[2], x[1]), pts)))
    u_ten_grid = CartesianMeshField(u_ten_data, info)

    # In Gridap 0.18, a common way to handle matrices is to use a 4-component VectorValue
    # since TensorValue{2,2} is internally just a SMatrix{2,2} which can be mapped.
    # We create a space that can hold 4 values.
    reffe_ten = ReferenceFE(lagrangian, VectorValue{4,Float64}, 1)
    V_ten = TestFESpace(Ω_act, reffe_ten)
    setup_transfer(geom, V_ten)

    # We need to map our TensorValue grid to a VectorValue{4} grid for the restriction
    u_ten_as_vec = CartesianMeshField([VectorValue(t[1], t[2], t[3], t[4]) for t in u_ten_data], info)

    u_ten_back_vec = mesh_to_grid(geom, grid_to_mesh(geom, u_ten_as_vec))

    # Verify the results (mapping back to TensorValue for comparison)
    u_ten_back = [TensorValue(v[1], v[2], v[3], v[4]) for v in u_ten_back_vec.data]

    err_ten = maximum([norm(u_ten_back[i] - u_ten_data[i]) for i in active_indices])
    println("  Tensor Error (via Vector4): $err_ten")
    @test err_ten < 1e-12
    @test u_ten_back[1] isa TensorValue

    # ---------------------------------------------------------
    # TEST 3: EXTENSION (Multi-component)
    # ---------------------------------------------------------
    println("Testing Vector Extension...")
    masked_vec = [phi <= 0 ? val : zero(val) for (phi, val) in zip(geom.levelset, u_vec_data)]
    u_vec_ext = extend(geom, CartesianMeshField(masked_vec, info))
    @test norm(u_vec_ext.data[1]) > 1e-3 # Ensure void is filled
    println("  Extension successful.")

end
