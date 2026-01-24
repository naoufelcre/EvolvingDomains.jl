using EvolvingDomains
using Gridap
using LinearAlgebra
using StaticArrays
using Test

@testset "Extension Physics & Accuracy" begin
    # 1. Setup
    domain = (0.0, 1.0, 0.0, 1.0)
    n = 40
    model = CartesianDiscreteModel(domain, (n, n))

    # 2. Circle in the middle
    radius = 0.2
    center = VectorValue(0.5, 0.5)
    geom = EvolvingDiscreteGeometry(model, x -> norm(VectorValue(x...) - center) - radius)

    # 3. Test Field: Constant velocity [1.0, 0.0] inside the circle
    # Grid has n cells -> n+1 nodes
    info = grid_info(geom.model)
    nx, ny = info.dims
    u_grid = fill(VectorValue(0.0, 0.0), nx, ny)
    
    # Define u_grid on the whole interior
    for j in 1:ny, i in 1:nx
        x_tuple = info.origin .+ (i-1, j-1) .* info.spacing
        x = VectorValue(x_tuple...)
        dist = norm(x - center)
        if dist <= radius + 1e-8
            u_grid[i,j] = VectorValue(1.0, 0.0)
        end
    end
    
    # 4. Extend
    u_ext = extend(geom, u_grid)

    # 5. Verify:
    # Point at [0.5, 0.5] is interior, should be [1.0, 0.0]
    # Point at [0.5, 0.5] is interior, should be [1.0, 0.0]
    # Use Tuples for index math to allow proper broadcasting
    pt_int = (0.5, 0.5)
    idx_int = floor.(Int, (pt_int .- info.origin) ./ info.spacing) .+ 1
    @test u_ext[idx_int[1], idx_int[2]] ≈ VectorValue(1.0, 0.0)

    # Outside point [0.8, 0.5] should project to interface -> see 1.0
    pt_out = (0.8, 0.5)
    idx_out = floor.(Int, (pt_out .- info.origin) ./ info.spacing) .+ 1
    # Compare norms explicitly to avoid issues with VectorValue isapprox
    val_out = u_ext[idx_out[1], idx_out[2]]
    target = VectorValue(1.0, 0.0)
    @test norm(val_out - target) < 0.1
end
