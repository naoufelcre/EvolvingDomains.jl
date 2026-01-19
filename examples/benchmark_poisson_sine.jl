
# Benchmark Script: Poisson Problem on "Crescent Moon"
# Comparing Quadtree-Adapted Mesh vs. Uniform Cartesian Mesh
# Solution: u(x,y) = sin(πx) * sin(πy) (Non-polynomial)
# =============================================================================

# Add package source to path
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

# Crescent Moon Geometry
const R_geo = 0.5
const L_geo = 0.5 * R_geo
const center1 = Point(0.0, 0.0)
const center2 = center1 + VectorValue(-L_geo, L_geo)

function crescent_sdf(x)
    p = Point(x...)
    d1 = norm(p - center1) - R_geo
    d2 = norm(p - center2) - R_geo
    return max(d1, -d2) 
end

# Analytical Solution: u(x) = sin(πx)sin(πy)
u_func(x) = sin(π*x[1]) * sin(π*x[2])
∇u_func(x) = VectorValue(
    π * cos(π*x[1]) * sin(π*x[2]),
    π * sin(π*x[1]) * cos(π*x[2])
)
# Laplacian defined as Δu = ∂xx + ∂yy
# ∂xx = -π^2 sin(πx)sin(πy)
# ∂yy = -π^2 sin(πx)sin(πy)
# Δu = -2π^2 u
# Equation: -Δu = f  => f = 2π^2 u
f_func(x) = 2 * π^2 * u_func(x)

# =============================================================================
# 2. Mesh Generation Functions (Same as before)
# =============================================================================

function generate_quadtree_model(L_refine::Int)
    h_fine = 1.0 / (2^L_refine)
    
    # Mapping [0,1] -> [-0.75, 0.75] (Length 1.5)
    origin = Point(-0.75, -0.75)
    width = 1.5
    map_to_physical(x_unit) = origin + width * VectorValue(x_unit[1], x_unit[2])
    
    # Sizing function in UNIT coordinates
    function sizing_in_unit(x_unit)
        x_phys = map_to_physical(x_unit)
        dist = abs(crescent_sdf(x_phys))
        return max(h_fine, 0.5*dist) 
    end

    # Generate fine base
    qmesh = generate_fine_mesh(1.0, L_refine) 
    
    # Coarsen
    bottom_up_coarsening!(qmesh, L_refine, sizing_in_unit)
    
    # Balance
    balance!(qmesh)
    
    # Pave 
    elements = pave_mesh(qmesh)
    
    # Convert to Gridap Model (on Unit Square)
    unit_model, _ = quadtree_to_discrete_model(elements)
    
    # Map Model to Physical Domain
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
# 3. Solver Workflow
# =============================================================================

function solve_single_iteration(bg_model, L_refine; vtk_prefix=nothing)
    # A. Setup Evolving Geometry
    N_lsm = ceil(Int, 1.5 * (2^(L_refine+1))) 
    domain_lsm = (-0.75, 0.75, -0.75, 0.75)
    lsm_model = CartesianDiscreteModel(domain_lsm, (N_lsm, N_lsm))
    
    t_geo_init = @elapsed begin
        geo = EvolvingDiscreteGeometry(bg_model, lsm_model, crescent_sdf)
    end
    
    # B. Update (Cut)
    t_update = @elapsed begin
        cut_geo = current_cut(geo)
        Ω_act = Triangulation(cut_geo, ACTIVE)
    end
    
    # C. Physics Solve (Poisson)
    l2_err = NaN
    h1_err = NaN
    area = NaN
    n_dofs = 0

    t_solve = @elapsed begin
        try
            # Spaces
            order = 1
            refFe_scalar = ReferenceFE(lagrangian, Float64, order)
            
            V_std = TestFESpace(Ω_act, refFe_scalar)
            
            strategy = AggregateCutCellsByThreshold(0.5)
            aggs = aggregate(strategy, cut_geo)
            V = AgFEMSpace(V_std, aggs)
            U = TrialFESpace(V)

            # Terms
            Ω = Triangulation(cut_geo, PHYSICAL)
            Γ = EmbeddedBoundary(cut_geo)
            
            degree = 2*order # Integration degree must handle sin(x)^2 reasonably well
            dΩ = Measure(Ω, degree+2) # Bump degree for Sine
            dΓ = Measure(Γ, degree+2)
            
            # Compute Area
            area = sum(∫(1)dΩ)
            
            n_Γ = get_normal_vector(Γ)
            γd = 10.0
            h = 1.5 / (2^L_refine) 
            
            # Weak form (Nitsche)
            a(u,v) = ∫( ∇(v)⋅∇(u) )dΩ + ∫( (γd/h)*v*u - v*(n_Γ⋅∇(u)) - (n_Γ⋅∇(v))*u )dΓ
            l(v) = ∫( f_func*v )dΩ + ∫( (γd/h)*v*u_func - (n_Γ⋅∇(v))*u_func )dΓ
            
            op = AffineFEOperator(a, l, U, V)
            uh = solve(op)
            
            # D. Error
            e = u_func - uh
            l2_err = sqrt(sum( ∫( e*e )dΩ ))
            h1_err = sqrt(sum( ∫( e*e + ∇(e)⋅∇(e) )dΩ ))
            
            # Visualization
            if !isnothing(vtk_prefix)
                output_dir = dirname(vtk_prefix)
                !isdir(output_dir) && mkdir(output_dir)
                writevtk(Ω, vtk_prefix * "_solution", cellfields=["uh"=>uh, "u_exact"=>u_func, "error"=>e])
                writevtk(Ω_act, vtk_prefix * "_active_mesh")
            end
            
            n_dofs = num_free_dofs(V)
            
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

# =============================================================================
# 4. Main Benchmark Control
# =============================================================================

function run_benchmark()
    println("Starting Benchmark: Quadtree vs Uniform (Sine)")
    println("==============================================")
    
    levels = 5:9
    
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
        m_q = solve_single_iteration(model_q, L; vtk_prefix="output/quadtree_sine_L$L")
        
        @printf("%-5d | %-10s | %-8d | %-8d | %-10.4f | %-10.4f | %-10.4f | %-10.2e | %-10.4f\n",
                L, "Quadtree", m_q.n_cells, m_q.n_dofs, t_gen_q, m_q.time_update, m_q.time_solve, m_q.l2_err, m_q.area)

        # Uniform
        GC.gc()
        t_gen_u = @elapsed model_u = generate_uniform_model(L)
        m_u = solve_single_iteration(model_u, L; vtk_prefix="output/uniform_sine_L$L")
        
        @printf("%-5d | %-10s | %-8d | %-8d | %-10.4f | %-10.4f | %-10.4f | %-10.2e | %-10.4f\n",
                L, "Uniform", m_u.n_cells, m_u.n_dofs, t_gen_u, m_u.time_update, m_u.time_solve, m_u.l2_err, m_u.area)
        
        println("-"^100)
    end
end

run_benchmark()
