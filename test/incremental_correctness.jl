
using EvolvingDomains
using Gridap
using Gridap.Geometry: get_node_coordinates
using GridapEmbedded
using Test
using LinearAlgebra

@testset "Incremental Integration Correctness (Symmetric API)" begin
    # 1. Setup
    n = 20
    model = CartesianDiscreteModel((0, 1, 0, 1), (n, n))
    
    # Moving Circle
    radius = 0.2
    center_0 = VectorValue(0.5, 0.5)
    velocity = VectorValue(0.1, 0.0)
    
    # SDF
    sdf(x, t) = norm(VectorValue(x...) - (center_0 + velocity * t)) - radius
    
    geom = EvolvingDiscreteGeometry(model, x->sdf(x, 0.0))
    
    # Integrator
    degree = 2
    integrator = IncrementalIntegrator(model, degree)
    
    # 2. Loop
    dt = 0.1
    n_steps = 5
    
    expected_area = π * radius^2
    
    for step in 1:n_steps
        t = step * dt
        
        # Move Geometry
        new_ls = vec([sdf(x, t) for x in get_node_coordinates(geom.lsm_model)])
        set_levelset!(geom, new_ls)
        
        # Update Integrator
        update_integrator!(integrator, geom)
        
        # Define Measure (User API)
        dΩ = measure_Ω(integrator, geom)
        
        # Integrate Area (Scalar) - NOT SUPPORTED by current IncrementalIntegrator (Matrix Cache only)
        # contrib = integrate(x->1.0, dΩ)
        # area = sum(contrib)
        # println("Step $step: Area = $area (Expected: $expected_area)")
        
        # Verify Matrix Assembly (Intended Use Case)
        V = FESpace(model, ReferenceFE(lagrangian, Float64, 1))
        U = TrialFESpace(V)
        
        a(u,v) = ∫(u*v)*dΩ
        A = assemble_matrix(a, U, V)
        
        # We can check the mass matrix sum corresponds roughly to area?
        # sum(A) ≈ Area
        total_mass = sum(A)
        println("Step $step: Mass Matrix Sum = $total_mass (Expected Area: $expected_area)")
        @test isapprox(total_mass, expected_area, atol=0.05)
        
        # Verify Boundary Measure works (length)
        dΓ = measure_Γ(integrator, geom)
        perimeter = sum(integrate(x->1.0, dΓ))
        expected_perm = 2 * π * radius
        println("Step $step: Perimeter = $perimeter (Expected: $expected_perm)")
        @test isapprox(perimeter, expected_perm, atol=0.1) 
    end
end
