
using EvolvingDomains
using Gridap
using Gridap.Geometry: get_node_coordinates
using GridapEmbedded
using Test
using LinearAlgebra
using SparseArrays

@testset "Incremental Assembly Verification" begin
    # 1. Setup
    n = 10
    model = CartesianDiscreteModel((0, 1, 0, 1), (n, n))
    
    # Static Circle (to avoid time stepping complexity)
    radius = 0.3
    center = VectorValue(0.5, 0.5)
    sdf(x) = norm(VectorValue(x...) - center) - radius
    
    geom = EvolvingDiscreteGeometry(model, sdf)
    
    # FESpace
    reffe = ReferenceFE(lagrangian, Float64, 1)
    V = FESpace(model, reffe)
    U = V
    
    # 2. Standard Assembly (Reference)
    cut = current_cut(geom)
    PHYSICAL = GridapEmbedded.Interfaces.PHYSICAL
    trian_cut = Triangulation(cut, PHYSICAL)
    dΩ_std = Measure(trian_cut, 2)
    
    # Mass Matrix form
    a(u,v) = ∫(u*v)dΩ_std
    M_std = assemble_matrix(a, U, V)
    
    # 3. Incremental Assembly
    integrator = IncrementalIntegrator(model, 2)
    update_integrator!(integrator, geom) # Populate cache
    
    dΩ_inc = measure_Ω(integrator, geom)
    
    a_inc(u,v) = ∫(u*v)dΩ_inc
    M_inc = assemble_matrix(a_inc, U, V)
    
    # 4. Compare
    diff = M_std - M_inc
    err = norm(diff, Inf)
    println("Assembly Error (Max Norm): ", err)
    
    # Note: Machine epsilon differences expected due to summation order?
    # Or strict equality?
    # Given consistent quadrature, should be very close.
    @test err < 1e-12
    
    # 5. Moving Circle Test (Evolution)
    # Move circle slightly -> Update -> Verify again
    new_sdf(x) = norm(x - (center + VectorValue(0.05, 0.0))) - radius
    new_ls = vec([new_sdf(x) for x in get_node_coordinates(geom.lsm_model)])
    set_levelset!(geom, new_ls)
    
    # Update Standard
    cut_2 = current_cut(geom)
    trian_cut_2 = Triangulation(cut_2, PHYSICAL)
    dΩ_std_2 = Measure(trian_cut_2, 2)
    M_std_2 = assemble_matrix(u->∫(u)*dΩ_std_2, U, V) # Linear form? No Bilinear
    # Gridap API quirk: a(u,v)
    a_2(u,v) = ∫(u*v)dΩ_std_2
    M_std_2 = assemble_matrix(a_2, U, V)
    
    # Update Incremental
    update_integrator!(integrator, geom)
    dΩ_inc_2 = measure_Ω(integrator, geom)
    a_inc_2(u,v) = ∫(u*v)dΩ_inc_2
    M_inc_2 = assemble_matrix(a_inc_2, U, V)
    
    err_2 = norm(M_std_2 - M_inc_2, Inf)
    println("Step 2 Error (Max Norm): ", err_2)
    @test err_2 < 1e-12
    
end
