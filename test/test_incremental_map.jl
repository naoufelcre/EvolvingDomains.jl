
using Pkg
Pkg.activate("/Users/naoufel/Documents/Scientific stack/EvolvingDomains.jl")

using Gridap
using GridapIncremental
using Test
using GridapIncremental.Utils: GeometryMap

# 1. Setup
model = CartesianDiscreteModel((0,1,0,1), (10,10))
reffe = ReferenceFE(lagrangian, Float64, 1) # Not needed for Measure?
trian = Triangulation(model)
dΩ = Measure(trian, 2)
cm = CachedMeasure(dΩ)

# 2. Initial Integration
# Integrate 1.0 everywhere.
vol_1 = sum(integrate(x->1.0, cm))
println("Vol 1: ", vol_1) # 1.0

# 3. Emulate Physics Change
# Let's say cells [1, 2] became "Dirty" (CUT or Status Change).
# All other cells are "Clean" (IN -> IN).
n_cells = num_cells(model)
cell_map = collect(1:n_cells) # Default identity
cell_map[1] = 0
cell_map[2] = 0

new_cells = [1, 2]
reused = collect(3:n_cells)
geo_map = GeometryMap(cell_map, new_cells, reused)

# 4. Update
update_measure!(cm, model, geo_map)

# 5. Check internals
println("Dirty Cells: ", cm.dirty_cells)
@test 1 in cm.dirty_cells
@test 2 in cm.dirty_cells
@test !(3 in cm.dirty_cells)

# 6. Verify Re-Integration works
# If we define a function that changed value in checking cells?
# Well, the measure logic usually assumes the domain didn't change if mapped.
# But for CachedMeasure, "Mapped" means "Copy Value".
# So if we mapped i->0, we forced it to be empty, so `integrate` will recompute it.
# Let's verify `integrate` recomputes dirty cells.
vol_2 = sum(integrate(x->1.0, cm))
println("Vol 2: ", vol_2) # Should still be 1.0

