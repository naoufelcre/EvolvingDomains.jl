
using Pkg
Pkg.activate("/Users/naoufel/Documents/Scientific stack/EvolvingDomains.jl")

using Gridap
using GridapIncremental
using Test

# 1. Setup Static Background
model = CartesianDiscreteModel((0,1,0,1), (10,10))

# 2. Emulate "View" Triangulation 1 (Left Half)
cells_1 = collect(1:50) # First 50 cells
trian_1 = Triangulation(model, cells_1) 
# Note: FESpace is usually built on the Model, then restricted?
# Or built on the Triangulation?
# Standard CutFEM: FESpace(model) -> restricted by integration measure?
# AgFEM: FESpace(model) -> aggregated DOFs.

# Let's try standard FESpace on model (CASE A: Global Space, shifting active domain)
reffe = ReferenceFE(lagrangian, Float64, 1)
V_global = FESpace(model, reffe) 
V_inc = IncrementalFESpace(model, reffe)

# CASE B: Space on a View?
# Gridap doesn't easily support FESpace on a View Triangulation directly without a Model context usually.
# But `RestrictedDiscreteModel` exists.

# Let's focus on CachedMeasure first, as that is the primary user request.
# CachedMeasure wraps a Measure.
# A Measure is typically `Measure(trian, degree)`.
degree = 2
dΩ_1 = Measure(trian_1, degree)
cm = CachedMeasure(dΩ_1)

# Integrate constant function
vol_1 = sum(integrate(x->1.0, cm))
println("Step 1 Volume: ", vol_1) # Should be 0.5

# 3. Emulate "View" Triangulation 2 (Left Half + 1 Cell)
cells_2 = collect(1:51)
trian_2 = Triangulation(model, cells_2)

# Update CachedMeasure?
# update_measure!(cm, new_model, geo_map)
# But `trian_2` is not a "new model". 
# The expected API for update_measure! takes a DiscreteModel and a GeometryMap (cell map).
# If we are just changing the *Integration Domain*, likely standard CachedMeasure doesn't support "changing triangulation" out of the box
# unless we treat the Triangulation change as a Model change.

# Let's see if we can trick it or if we need to overload it.
# If we pass `model` (background) to update_measure!, it assumes the measure is over `model`.
# But our measure is over `trian_1`.

# HYPOTHESIS: CachedMeasure needs to be defined on the Background Mesh, 
# and we provide a "Dirty Mask" to restrict integration to the "Active" part + "Dirty" part.
# But for CutFEM, the "Active" part isn't just a list of cells, the GEOMETRY inside the cells changes.

println("Done.")
