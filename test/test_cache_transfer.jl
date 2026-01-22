using Test
using Gridap
using EvolvingDomains
import LevelSetMethods as LSM
import TransferOperator: prolong, restrict
import EvolvingDomains: get_transfer_op

@testset "Transfer Cache Memoization" begin
    # 1. Setup minimal geometry
    model = CartesianDiscreteModel((0,1,0,1), (4,4))
    geom = EvolvingDiscreteGeometry(model, x -> norm(x) - 0.5)
    
    # Initial state: step = 0, cache.transfer_step = -1
    @test geom.step == 0
    @test geom.cache.transfer_step == -1
    
    # 2. First prolong call should build the cache
    V = FESpace(model, ReferenceFE(lagrangian, Float64, 1))
    u_mesh = FEFunction(V, ones(num_free_dofs(V)))
    u_grid = prolong(geom, u_mesh)
    
    initial_op = geom.cache.transfer_op
    @test initial_op !== nothing
    @test geom.cache.transfer_step == 0
    
    # 3. Second call in same step should REUSE the cache
    u_grid_2 = prolong(geom, u_mesh)
    @test geom.cache.transfer_op === initial_op
    @test geom.cache.transfer_step == 0
    
    # 4. Advance! should invalidate (logically)
    advance!(geom, 0.1)
    @test geom.step == 1
    @test geom.cache.transfer_step == 0 # Still 0, but geom.step is 1
    
    # 5. prolong after advance should REBUILD
    u_grid_3 = prolong(geom, u_mesh)
    @test geom.cache.transfer_step == 1
    
    # 6. Manual Override (Expert Bridge)
    # Even if we are at step 1, if we manually write a new op and set step to 1,
    # it should be picked up.
    custom_op = setup_transfer(model, model)
    geom.cache.transfer_op = custom_op
    geom.cache.transfer_step = 1
    
    @test get_transfer_op(geom) === custom_op
end
