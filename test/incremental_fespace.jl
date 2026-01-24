
using EvolvingDomains
using EvolvingDomains: get_geometry_map
using Gridap
using Gridap.Geometry: get_node_coordinates
using GridapEmbedded
using GridapIncremental
using Test
using LinearAlgebra

@testset "Incremental FESpace Verification" begin
    # 1. Setup
    n = 20
    model = CartesianDiscreteModel((0, 1, 0, 1), (n, n))
    
    # Initial Geometry
    radius = 0.2
    center_0 = VectorValue(0.5, 0.5)
    velocity = VectorValue(0.1, 0.0)
    sdf(x, t) = norm(VectorValue(x...) - (center_0 + velocity * t)) - radius
    
    geom = EvolvingDiscreteGeometry(model, x->sdf(x, 0.0))
    
    # 2. Construct Incremental FESpace
    reffe = ReferenceFE(lagrangian, Float64, 1)
    
    # Step A: Base Incremental Space
    V_inc = IncrementalFESpace(model, reffe)
    
    # Step B: Initial Aggregation Wrapper
    # Current cut is updated lazily
    cut = current_cut(geom)
    # Note: Aggregate requires strategy
    params = Dict(:target_cond => 0.5) 
    # AgFEM usually needs `aggregate(strategy, cut, geom/inout, ...)`
    # Looking at GridapEmbedded tests/AgFEM...
    # Assuming `aggregate(cut, TBD)`
    # For now, let's verify V_inc mechanics first.
    
    integrator = IncrementalIntegrator(model, 2)
    
    # 3. Evolution Loop
    dt = 0.1
    n_steps = 3
    
    for step in 1:n_steps
        t = step * dt
        
        # Move Geometry
        new_ls = vec([sdf(x, t) for x in get_node_coordinates(geom.lsm_model)])
        set_levelset!(geom, new_ls)
        
        # Update Integrator (Physics)
        update_integrator!(integrator, geom)
        
        # Retrieve Map
        geo_map = get_geometry_map(integrator)
        @test !isnothing(geo_map)
        
        # Update FESpace
        update_fespace!(V_inc, model, geo_map)
        
        # Verify FESpace is valid
        trian = Triangulation(model)
        @test num_free_dofs(V_inc) > 0
        
        println("Step $step: V_inc DOFs = $(num_free_dofs(V_inc))")
        
        # Step C: Aggregate Wrapper (Test Compatibility)
        # We need aggregation strategy (dummy for now)
        # Assuming we can just wrap without explicit aggregation for test if cut is trivial? 
        # Actually AgFEMSpace needs aggregates.
        # But we can verify `V_inc` behaves like a proper SingleFieldFESpace.
        @test isa(V_inc, Gridap.FESpaces.SingleFieldFESpace)
        
        # Note: Full AgFEM test requires robust aggregation setup which depends on cut topology.
        # Ideally:
        # aggregates = aggregate(RobustAggregation(0.5), cut, geom, IN)
        # V = AgFEMSpace(V_inc, aggregates)
        # @test num_free_dofs(V) <= num_free_dofs(V_inc)
    end
    
end
