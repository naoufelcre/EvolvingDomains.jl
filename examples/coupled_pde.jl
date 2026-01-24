# =============================================================================
# coupled_pde.jl — Skeleton for Dynamic FEM-Driven Evolution
# =============================================================================

using EvolvingDomains
import TransferOperator: prolong
using Gridap
using GridapEmbedded
using GridapIncremental  # [NEW] For IncrementalFESpace
using LevelSetMethods
using CairoMakie
using LinearAlgebra: norm

# 1. Setup Standard Model (Single mesh for physics and LSM)
domain = (0.0, 1.0, 0.0, 1.0)
n_cells = (64, 64)
model = CartesianDiscreteModel(domain, n_cells)

# 2. Initial Geometry
center = VectorValue(0.5, 0.5)
radius = 0.125
geom = EvolvingDiscreteGeometry(model,
    x -> norm(VectorValue(x...) - center) - radius;
    reinit_freq = 5
)

# [OPTIMIZATION 1] Incremental FESpace
# Instead of rebuilding the FESpace from scratch every step (which involves 
# expensive DOF re-numbering and allocation), we use an IncrementalFESpace.
# This space maintains a map of DOFs and only updates the "new" cells coming 
# from the moving interface, reusing the bulk DOFs.
reffe = ReferenceFE(lagrangian, Float64, 1)
V_inc = IncrementalFESpace(model, reffe) 

# [OPTIMIZATION 2] Cached Measure (IncrementalIntegrator)
# Standard `Measure(Ω, degree)` recomputes integration points for the *entire* domain.
# IncrementalIntegrator caches the quadrature for cells that remain purely "Inside"
# between steps, and only recomputes for cells near the interface (CUT) or 
# those that changed status.
integrator = IncrementalIntegrator(model, 4)

# 3. Dynamic Simulation Loop
dt = 0.01
t_end = 0.75
n_steps = round(Int, t_end / dt)

fig = Figure(size=(700, 600))
ax = Axis(fig[1, 1], aspect=1, limits=(0, 1, 0, 1),
          xlabel="x", ylabel="y", title="Coupled PDE Evolution")

output_path = joinpath(@__DIR__, "coupled_pde.gif")

function step_simulation(step)
    t = step * dt
    
    # A. Update Domain Information (Cut logic)
    # returns an EmbeddedDiscretization (cached inside geom)
    cut = current_cut(geom)
    
    # [OPTIMIZATION] Update Physics-Aware Machinery
    # 1. Update the Integrator (computes the "GeometryMap" of cell changes)
    update_integrator!(integrator, geom)
    
    # 2. Retrieve the change map (which cells are NEW, which are REUSED)
    geo_map = get_geometry_map(integrator)
    
    # 3. Smart Update of FESpace using the change map
    # This is O(Δinterface), vastly faster than O(TOTAL) for large meshes.
    update_fespace!(V_inc, model, geo_map)

    # B. Define Physics Domain and Measure
    # Use the cached measure instead of creating a fresh one.
    dΩ = measure_Ω(integrator, geom)
    
    # C. TODO: Solve your dynamic PDE on (Ω, dΩ)
    # Example: 
    # V = RobustAggFESpace(model, reffe, cut, aggregate(RobustAggregation(0.5), cut, geo, IN))
    # OR directly use V_inc if suitable (Incremental space handles raw cut DOFs)
    # uh = solve_physics(V_inc, dΩ, ...)
    
    # Placeholder: Simple analytic motion
    uh = x -> VectorValue(0.1, 0.0) 
    
    # D. Derive Advection Velocity
    # Velocity can be any function of the FE solution: v = f(uh)
    v_phys = x -> uh(x)
    
    # E. Transfer Workflow (Physics Domain -> LSM Grid)
    # prolong: pointwise sampling on grid nodes (Optimized internally via TransferOperator cache)
    # extend: closest-point extension for stability (Optimized internally via ExtensionOperator cache)
    # set_velocity!: writes to LSM buffer
    v_grid = prolong(geom, v_phys) 
    v_ext  = extend(geom, v_grid)  
    set_velocity!(geom, v_ext)     
    
    # F. Advance Geometry
    step > 0 && advance!(geom, dt)
    
    # G. Visualize
    empty!(ax)
    plot_levelset!(ax, geom)
    ax.title = "Step $step (t=$(round(t, digits=2)))"
end

println("Starting simulation...")
println("  - Using IncrementalFESpace (Smart DOF Reuse)")
println("  - Using IncrementalIntegrator (Cached Quadrature)")
record(step_simulation, fig, output_path, 0:n_steps; framerate=20)
println("✅ Animation saved to: $output_path")
