module TestDumbellParabolic

using Test
using EvolvingDomains
using EvolvingDomains.Geometric
using EvolvingDomains.Kinematic
using EvolvingDomains.Transfer
using Gridap
using Gridap.Geometry
using Gridap.TensorValues
using GridapEmbedded

# This test Implements a numerical test inspired by the paper
# https://arxiv.org/pdf/2504.14116
# Analysis of a finite element method for PDEs in evolving domains with topological changes
# MA Olshanskii, A Reusken

@testset "Dumbbell Parabolic Benchmark (AgFEM)" begin

    # --- 1. Setup ---
    domain = (-2, 2, -2.0, 2.0)
    n = 75
    partition = (n, n)
    grid = CartesianDiscreteModel(domain, partition)
    info = grid_info(grid)

    # Analytical Level Set
    # φ(t, x) < 0 is inside
    function analytical_levelset(t, x)
        x1, x2 = x[1], x[2]
        return (t - 0.25) - ((x1^2 - 0.3 * x1^4) - x2^2)
    end

    # Initial Condition
    function u0_func(x)
        x1, x2 = x[1], x[2]
        return analytical_levelset(0, x)^2 + (x1 + x2 + 1)
    end

    function heat_varforms(u_prev, Δt, dΩ)
        a(u, v) = ∫((u * v) / Δt + (∇(u) ⋅ ∇(v)))dΩ
        l(v) = ∫((u_prev * v) / Δt)dΩ
        return a, l
    end

    # --- 2. Initial State (t=0) ---
    t = 0.0
    pts = Gridap.Geometry.get_node_coordinates(grid)

    # Initialize Geometry and initial temperature dsitribution
    phi_0 = map(x -> analytical_levelset(t, x), pts)
    u0 = map(x -> u0_func(x), pts)

    geom = EvolvingDiscreteGeometry(vec(collect(phi_0)), grid)
    u_grid = CartesianMeshField(vec(collect(u0)), info)

    # --- 3. Time Loop ---
    nTime = 200
    timeHorizon = 0.5
    Δt = timeHorizon / nTime

    plot_interval = 1
    plot(geom; label="Dumbbell step 0 / $nTime   t = 0.0")
    started = time_ns()

    for step in 1:nTime
        t += Δt
        # Update Geometry
        phi_new = map(x -> analytical_levelset(t, x), pts)
        set_levelset!(geom, vec(collect(phi_new)))
        get_active_indices(geom, :current)
        cutgeo = geom.cache.cut

        # Cutting
        Ω = Triangulation(cutgeo, PHYSICAL)
        Ω_act = Triangulation(cutgeo, ACTIVE)

        #FESpace management
        order = 1
        reffe = ReferenceFE(lagrangian, Float64, order)
        Vstd = TestFESpace(Ω_act, reffe, conformity=:H1)

        strategy = AggregateAllCutCells()
        aggregates = aggregate(strategy, cutgeo)

        V = AgFEMSpace(Vstd, aggregates)
        U = TrialFESpace(V)

        #Measure management
        degree = 2 * order
        dΩ = Measure(Ω, degree)

        # Setup transfer operator for the CURRENT geometry and AgFEM space
        transfer_op = setup_transfer(geom, V)

        # Restrict u_grid (from previous step) to the new mesh
        # This effectively interpolates u_n onto Ω_{n+1}
        u_prev = grid_to_mesh(geom, u_grid)

        #Solve Heat Equation
        a, l = heat_varforms(u_prev, Δt, dΩ)
        op = AffineFEOperator(a, l, U, V)
        u_h = solve(op)

        # Store solution back to grid
        u_grid = mesh_to_grid(geom, u_h)

        if step % plot_interval == 0 || step == nTime
            mean_ms = (time_ns() - started) / (1e6 * step)
            plot(geom; label="Dumbbell step $step / $nTime   t = $(round(t, digits=4))   mean = $(round(mean_ms, digits=2)) ms/iteration")
        end
    end

    @test true
end

end # module
