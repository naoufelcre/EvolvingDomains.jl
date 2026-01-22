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
    u_grid = fill(VectorValue(0.0, 0.0), n, n)
    info = grid_info(geom)
    
    # Define u_grid on the whole interior
    for j in 1:n, i in 1:n
        x = info.origin .+ (i-1, j-1) .* info.spacing
        dist = norm(x - center)
        if dist <= radius + 1e-8
            u_grid[i,j] = VectorValue(1.0, 0.0)
        end
    end
    
    # 4. Extend
    u_ext = extend(geom, u_grid)

    # 5. Verify:
    # Point at [0.5, 0.5] is interior, should be [1.0, 0.0]
    idx_int = floor.(Int, (VectorValue(0.5, 0.5) .- info.origin) ./ info.spacing) .+ 1
    @test u_ext[idx_int[1], idx_int[2]] ≈ VectorValue(1.0, 0.0)

    # Outside point [0.8, 0.5] should project to interface -> see 1.0
    idx_out = floor.(Int, (VectorValue(0.8, 0.5) .- info.origin) ./ info.spacing) .+ 1
    @test u_ext[idx_out[1], idx_out[2]] ≈ VectorValue(1.0, 0.0) atol=0.1
end
