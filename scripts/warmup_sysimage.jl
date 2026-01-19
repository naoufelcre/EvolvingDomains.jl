# =============================================================================
# warmup_sysimage.jl — Exercise code paths to trigger compilation
# =============================================================================
#
# This script runs during sysimage build to precompile common operations.
# Add any workflows you frequently use to minimize runtime JIT.
#
# =============================================================================

@info "Warming up EvolvingDomains..."

using EvolvingDomains
using Gridap
using GridapEmbedded
using StaticArrays

# -----------------------------------------------------------------------------
# 1. Basic Geometry Creation & Evolution
# -----------------------------------------------------------------------------
@info "  Compiling geometry operations..."

model = CartesianDiscreteModel((0,1,0,1), (16,16))

# Circle SDF
sdf = x -> sqrt((x[1]-0.5)^2 + (x[2]-0.5)^2) - 0.25

# Create geometry with velocity
vel = TimeDependentVelocity((x, t) -> (-x[2] + 0.5, x[1] - 0.5))
geom = EvolvingDiscreteGeometry(model, sdf; velocity=vel, reinit_freq=10)

# Advance a few steps
for _ in 1:3
    advance!(geom, 0.01)
end

# Access cut geometry (triggers GridapEmbedded compilation)
cut = current_cut(geom)
ϕ = current_levelset(geom)

# -----------------------------------------------------------------------------
# 2. Quadtree Mesh Generation
# -----------------------------------------------------------------------------
@info "  Compiling quadtree operations..."

mesh = generate_fine_mesh(1.0, 5)
sizing = sizing_function_from_distance(sdf, 0.3, 0.02)
bottom_up_coarsening!(mesh, 5, sizing)
balance!(mesh)
elements = pave_mesh(mesh)

# Convert to Gridap model
qt_model = quadtree_to_discrete_model(elements)

# -----------------------------------------------------------------------------
# 3. FE Spaces and Assembly (major compilation target)
# -----------------------------------------------------------------------------
@info "  Compiling FE assembly..."

# Simple Poisson on the cut domain
Ω_act = Triangulation(cut, ACTIVE)
Ω_phys = Triangulation(cut, PHYSICAL_IN)
Γ = EmbeddedBoundary(cut)

order = 1
reffe = ReferenceFE(lagrangian, Float64, order)

# Aggregated FE Space
strategy = AggregateCutCellsByThreshold(1.0)
aggregates = aggregate(strategy, cut)
Vstd = FESpace(Ω_act, reffe; conformity=:H1)
V = AgFEMSpace(Vstd, aggregates)
U = TrialFESpace(V)

# Measures
degree = 2*order
dΩ = Measure(Ω_phys, degree)
dΓ = Measure(Γ, degree)

# Weak form (Poisson with Nitsche BC)
n_Γ = get_normal_vector(Γ)
h = 1/16
γ = 10.0 / h

a(u,v) = ∫( ∇(u)⋅∇(v) )dΩ + 
         ∫( γ*(u*v) - (∇(u)⋅n_Γ)*v - u*(∇(v)⋅n_Γ) )dΓ

f(x) = 1.0
g(x) = 0.0

l(v) = ∫( f*v )dΩ + ∫( γ*g*v - g*(∇(v)⋅n_Γ) )dΓ

# Assemble and solve
op = AffineFEOperator(a, l, U, V)
uh = solve(op)

# -----------------------------------------------------------------------------
# 4. Velocity Sources
# -----------------------------------------------------------------------------
@info "  Compiling velocity sources..."

# Static velocity
v_static = StaticFunctionVelocity(x -> (1.0, 0.0))
sample_velocity(v_static, (0.5, 0.5), 0.0)

# Time-dependent velocity
v_time = TimeDependentVelocity((x, t) -> (sin(t), cos(t)))
sample_velocity(v_time, (0.5, 0.5), 1.0)

# FE velocity source (if we had a vector FE function)
# This path is covered by Stokes solves in practice

# -----------------------------------------------------------------------------
# Done
# -----------------------------------------------------------------------------
@info "✅ Warmup complete — all major code paths compiled"
