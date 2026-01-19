
# Benchmark Script: Poisson Problem on "Crescent Moon"
# Comparing Quadtree-Adapted Mesh vs. Uniform Cartesian Mesh
# Solution: Localized Gaussian Peak on Interface
# =============================================================================

push!(LOAD_PATH, joinpath(@__DIR__, "../src"))

using Gridap
using Gridap.Geometry
using GridapEmbedded
using EvolvingDomains
using LinearAlgebra
using Printf
using Dates

# =============================================================================
# 1. Geometry & Manufactured Solution
# =============================================================================

const R_geo = 0.5 + 1.0e-5  # Perturbed to avoid grid alignment singularities
const L_geo = 0.5 * R_geo
const center1 = Point(0.0, 0.0)
const center2 = center1 + VectorValue(-L_geo, L_geo)

function crescent_sdf(x)
    p = Point(x...)
    d1 = norm(p - center1) - R_geo
    d2 = norm(p - center2) - R_geo
    return max(d1, -d2) 
end

# Gaussian Peak parameters
const x0 = Point(0.5, 0.0) # On the interface of Circle 1
const α = 100.0            # Sharp decay

# Analytical Solution: u(x) = exp(-α * |x - x0|^2)
u_func(x) = exp(-α * norm(x - x0)^2)

# Gradient: ∇u = -2α(x-x0) * u
∇u_func(x) = -2 * α * (x - x0) * u_func(x)

# Laplacian: Δu = ∇⋅∇u
# = -2α [ ∇⋅((x-x0)u) ]
# = -2α [ (∇⋅(x-x0))u + (x-x0)⋅∇u ]
# = -2α [ 2u + (x-x0)⋅(-2α(x-x0)u) ]
# = -2α [ 2u - 2α|x-x0|^2 u ]
# = -4α u + 4α^2 |x-x0|^2 u
# = 4α u ( α|x-x0|^2 - 1 )

# Poisson Equation: -Δu = f
f_func(x) = - (4 * α * u_func(x) * (α * norm(x - x0)^2 - 1))

# =============================================================================
# 2. Mesh Generation
# =============================================================================

function generate_quadtree_model(L_refine::Int)
    h_fine = 1.0 / (2^L_refine)
    
    # Mapping [0,1] -> [-0.75, 0.75]
    epsilon = 1.0e-7  # Perturbation to avoid grid alignment singularities (L=6 fix)
    origin = Point(-0.75 + epsilon, -0.75 + epsilon)
    width = 1.5
    map_to_physical(x_unit) = origin + width * VectorValue(x_unit[1], x_unit[2])
    
    function sizing_in_unit(x_unit)
        x_phys = map_to_physical(x_unit)
        
        # New Logistic Sizing Function:
        # - h_min = h_fine (finest level)
        # - h_max = 0.5 (very coarse in bulk)
        # - steepness = 20.0 (smooth but firm transition)
        # - transition = 0.3 (keep fine resolution until 0.3 distance away)
        sizer = logistic_sizing_function(crescent_sdf, h_fine, 0.5, 20.0, 0.3)
        
        return sizer(x_phys)
    end

    qmesh = generate_fine_mesh(1.0, L_refine) 
    bottom_up_coarsening!(qmesh, L_refine, sizing_in_unit)
    balance!(qmesh)
    elements = pave_mesh(qmesh)
    # valid_elements = [e for e in elements if !any(isnan, e)] # Removing this as QuadElement is not iterable

    unit_model, _ = quadtree_to_discrete_model(elements)
    
    grid_unit = get_grid(unit_model)
    nodes = Gridap.Geometry.get_node_coordinates(grid_unit)
    new_nodes = [map_to_physical(n) for n in nodes]
    
    grid_phys = UnstructuredGrid(new_nodes, 
                                 get_cell_node_ids(grid_unit), 
                                 get_reffes(grid_unit), 
                                 get_cell_type(grid_unit), 
                                 Gridap.Geometry.NonOriented())
    bg_model = UnstructuredDiscreteModel(grid_phys)
    
    return bg_model
end

function generate_uniform_model(L_refine::Int)
    N = ceil(Int, 1.5 * (2^L_refine))
    domain = (-0.75, 0.75, -0.75, 0.75)
    partition = (N, N)
    model = CartesianDiscreteModel(domain, partition)
    return simplexify(model)
end

# =============================================================================
# 3. Solver
# =============================================================================

function solve_single_iteration(bg_model, L_refine; vtk_prefix=nothing)
    N_lsm = ceil(Int, 1.5 * (2^(L_refine+1))) 
    domain_lsm = (-0.75, 0.75, -0.75, 0.75)
    lsm_model = CartesianDiscreteModel(domain_lsm, (N_lsm, N_lsm))
    
    t_geo_init = @elapsed begin
        geo = EvolvingDiscreteGeometry(bg_model, lsm_model, crescent_sdf)
    end
    
    t_update = @elapsed begin
        cut_geo = current_cut(geo)
        Ω_act = Triangulation(cut_geo, ACTIVE)
    end
    
    l2_err = NaN
    h1_err = NaN
    area = NaN
    n_dofs = 0

    t_solve = @elapsed begin
        try
            order = 1
            refFe_scalar = ReferenceFE(lagrangian, Float64, order)
            
            V_std = TestFESpace(Ω_act, refFe_scalar)
            
            strategy = AggregateCutCellsByThreshold(0.5)
            aggs = aggregate(strategy, cut_geo)
            V = AgFEMSpace(V_std, aggs)
            U = TrialFESpace(V)

            Ω = Triangulation(cut_geo, PHYSICAL)
            Γ = EmbeddedBoundary(cut_geo)
            
            degree = 2*order + 2
            dΩ = Measure(Ω, degree)
            dΓ = Measure(Γ, degree)
            
            area = sum(∫(1)dΩ)
            
            n_Γ = get_normal_vector(Γ)
            γd = 10.0
            h = 1.5 / (2^L_refine) 
            
            a(u,v) = ∫( ∇(v)⋅∇(u) )dΩ + ∫( (γd/h)*v*u - v*(n_Γ⋅∇(u)) - (n_Γ⋅∇(v))*u )dΓ
            l(v) = ∫( f_func*v )dΩ + ∫( (γd/h)*v*u_func - (n_Γ⋅∇(v))*u_func )dΓ
            
            op = AffineFEOperator(a, l, U, V)
            uh = solve(op)
            
            n_dofs = num_free_dofs(V)

            # D. Error (Manual calculation to avoid AD issues)
            # L2: ∫ (u - uh)^2
            l2_err = sqrt(sum( ∫( (u_func - uh)*(u_func - uh) )dΩ ))
            
            # H1: ∫ (u - uh)^2 + (∇u - ∇uh)^2
            # We use the analytical gradient ∇u_func
            diff_grad(x) = ∇u_func(x) - ∇(uh)(x)
            # But ∇(uh) is a Gridap object. 
            # Easier: Use integrated terms
            # ∇(e)⋅∇(e) = |∇u|^2 - 2∇u⋅∇uh + |∇uh|^2
            # Or define a function for gradient difference?
            # Simplest: ∫( ∇(uh)⋅∇(uh) - 2*∇(uh)⋅∇u_func + ∇u_func⋅∇u_func )dΩ
            # Note: ∇u_func returns a VectorValue matching Gridap's expectations?
            
            h1_semiconorm_sq = sum( ∫( (∇(uh) - ∇u_func)⋅(∇(uh) - ∇u_func) )dΩ )
            h1_err = sqrt(l2_err^2 + h1_semiconorm_sq)
            
            e = u_func - uh
            
            if !isnothing(vtk_prefix)
                output_dir = dirname(vtk_prefix)
                !isdir(output_dir) && mkdir(output_dir)
                writevtk(Ω, vtk_prefix * "_solution", cellfields=["uh"=>uh, "u_exact"=>u_func, "error"=>e])
                writevtk(Ω_act, vtk_prefix * "_active_mesh")
            end
            
        catch e
            println("    [Error] Solver failed: ", e)
            isa(e, InterruptException) && rethrow(e)
        end
    end
    
    n_cells = num_cells(bg_model)
    
    return (
        time_update = t_update,
        time_solve = t_solve,
        n_cells = n_cells,
        n_dofs = n_dofs,
        l2_err = l2_err,
        h1_err = h1_err,
        area = area
    )
end

function run_benchmark()
    println("Starting Benchmark: Quadtree vs Uniform (Localized Gaussian)")
    println("==========================================================")
    
    levels = [5, 7, 8, 9]  # Skipping L=6 (known aggregation bug)
    
    println("\n[Warmup] Pre-compiling functions with L=4...")
    solve_single_iteration(generate_quadtree_model(4), 4)
    solve_single_iteration(generate_uniform_model(4), 4)
    println("[Warmup] Complete. Starting Measurements.\n")
    
    @printf("%-5s | %-10s | %-8s | %-8s | %-10s | %-10s | %-10s | %-10s | %-10s\n", 
            "L", "Type", "Cells", "DoFs", "Gen(s)", "Upd(s)", "Sol(s)", "L2 Error", "Area")
    println("-"^100)
    
    for L in levels
        # Quadtree
        GC.gc()
        t_gen_q = @elapsed model_q = generate_quadtree_model(L)
        m_q = solve_single_iteration(model_q, L; vtk_prefix="output/quadtree_local_L$L")
        
        @printf("%-5d | %-10s | %-8d | %-8d | %-10.4f | %-10.4f | %-10.4f | %-10.2e | %-10.4f\n",
                L, "Quadtree", m_q.n_cells, m_q.n_dofs, t_gen_q, m_q.time_update, m_q.time_solve, m_q.l2_err, m_q.area)

        # Uniform
        GC.gc()
        t_gen_u = @elapsed model_u = generate_uniform_model(L)
        m_u = solve_single_iteration(model_u, L; vtk_prefix="output/uniform_local_L$L")
        
        @printf("%-5d | %-10s | %-8d | %-8d | %-10.4f | %-10.4f | %-10.4f | %-10.2e | %-10.4f\n",
                L, "Uniform", m_u.n_cells, m_u.n_dofs, t_gen_u, m_u.time_update, m_u.time_solve, m_u.l2_err, m_u.area)
        
        println("-"^100)
    end
end

run_benchmark()
