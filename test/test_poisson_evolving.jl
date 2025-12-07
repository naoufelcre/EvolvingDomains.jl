# =============================================================================
# Test: Poisson Problem on Evolving Geometry
# =============================================================================
# This test solves a Poisson problem on a domain defined by an evolving level set.
#
# Problem: -Δu = f in Ω(t)
#           u = g on ∂Ω(t)
#
# We use a manufactured solution to verify the accuracy.

using Gridap
using GridapEmbedded
using LinearAlgebra
using Printf
using Test

# Add the local modules
push!(LOAD_PATH, joinpath(@__DIR__, "..", "src"))
using EvolvingDomains

# Try to load LevelSetMethods
lsm_available = try
    using LevelSetMethods
    true
catch
    false
end

if lsm_available
    println("=" ^ 70)
    println("       POISSON ON EVOLVING DOMAIN TEST")
    println("=" ^ 70)

    # =========================================================================
    # 1. SETUP
    # =========================================================================

    # Background model
    domain = (-1.0, 1.0, -1.0, 1.0)
    partition = (40, 40)
    model = CartesianDiscreteModel(domain, partition)

    # Initial level set (circle at (0.5, 0.0))
    R = 0.2
    center = (0.5, 0.0)
    ϕ0(x) = norm(x .- center) - R # Negative inside, positive outside

    # Velocity (rotation)
    u_vel(x) = (-x[2], x[1])

    # Create evolver
    evolver = LevelSetMethodsEvolver(;
        bg_model = model,
        initial_ls = ϕ0,
        velocity = u_vel,
        integrator = :RK3,
        reinit_freq = nothing,
        bc = :Neumann
    )

    # Create evolving geometry
    eg = EvolvingDiscreteGeometry(evolver, model)

    # =========================================================================
    # 2. MANUFACTURED SOLUTION
    # =========================================================================

    # u(x, t) = sin(πx) * sin(πy) * cos(t)
    u_exact(x, t) = sin(π*x[1]) * sin(π*x[2]) * cos(t)

    # f(x, t) = -Δu = -(-2π² * u) = 2π² * u
    f_source(x, t) = 2 * π^2 * u_exact(x, t)

    # g(x, t) = u(x, t) on boundary
    g_bc(x, t) = u_exact(x, t)

    # =========================================================================
    # 3. TIME LOOP
    # =========================================================================

    Δt = 0.1
    nsteps = 100

    println("Running $nsteps steps with Δt = $Δt...")

    for step in 1:nsteps
        # Advance geometry
        advance!(eg, Δt)
        t = EvolvingDomains.current_time(eg)

        # Get current geometry and cut
        geo = current_geometry(eg)
        cut_geo = cut(model, geo)

        # Setup integration meshes
        # Interior Ω
        trian_in = Triangulation(cut_geo, PHYSICAL_IN)
        dΩ = Measure(trian_in, 2)

        # Boundary Γ (interface)
        trian_Γ = EmbeddedBoundary(cut_geo)
        dΓ = Measure(trian_Γ, 2)

        # Boundary Γ_bg (background boundary, if any active)
        # For this problem, the circle is fully inside, so we only care about the interface.
        # But strictly speaking we might need to handle background boundaries if the domain touches them.
        # Here we keep the circle small enough.

        # FE Spaces
        # Standard H1 space on background
        reffe = ReferenceFE(lagrangian, Float64, 1)
        V = TestFESpace(model, reffe, conformity=:H1)
        U = TrialFESpace(V)

        # Weak form (Nitsche's method for Dirichlet BC on Γ)
        γ = 10.0 # Penalty parameter
        h = 2.0 / 40.0 # Approximate mesh size

        # Functions at current time
        u_ex(x) = u_exact(x, t)
        f(x) = f_source(x, t)
        g(x) = g_bc(x, t)

        # Background triangulation for regularization
        trian_bg = Triangulation(model)
        dΩ_bg = Measure(trian_bg, 2)

        # Nitsche terms
        # a(u, v) = ∫∇u⋅∇v dΩ - ∫(∇u⋅n)v dΓ - ∫(∇v⋅n)u dΓ + ∫(γ/h)uv dΓ
        # l(v) = ∫fv dΩ - ∫(∇v⋅n)g dΓ + ∫(γ/h)gv dΓ

        n_Γ = get_normal_vector(trian_Γ)

        # Define scaled functions to avoid Number * Function error
        penalty_g(x) = (γ/h) * g(x)

        # Add small regularization 1e-10 on background to fix singular matrix (outside DOFs)
        a(u, v) = ∫( ∇(u)⋅∇(v) )dΩ +
                  ∫( (γ/h)*(u*v) - (n_Γ⋅∇(u))*v - (n_Γ⋅∇(v))*u )dΓ +
                  ∫( 1.0e-10*(u*v) )dΩ_bg

        l(v) = ∫( f*v )dΩ +
               ∫( penalty_g*v - (n_Γ⋅∇(v))*g )dΓ

        # Solve
        op = AffineFEOperator(a, l, U, V)
        uh = solve(op)

        # Compute error
        e = u_ex - uh
        l2_err = sqrt(sum(∫( e*e )dΩ))

        @printf("Step %d (t=%.2f): L2 Error = %.2e\n", step, t, l2_err)

        # Assert error is small (expecting ~O(h^2) or at least reasonable)
        @test l2_err < 1e-2

        # Visualization
        # Save solution and error to VTK
        # Note: We need to interpolate the exact solution to the FE space to visualize it easily on the same mesh
        u_interp = interpolate_everywhere(u_ex, U)

        # Write to VTK
        # We save the cut mesh (trian_in) with the solution and error fields
        writevtk(trian_in, "poisson_evolving_step_$(step)", cellfields=["uh"=>uh, "u_exact"=>u_interp, "error"=>e])

    end

    println("=" ^ 70)
    println("       TEST PASSED ✓")
    println("=" ^ 70)

else
    println("LevelSetMethods.jl not available. Skipping Poisson test.")
end
