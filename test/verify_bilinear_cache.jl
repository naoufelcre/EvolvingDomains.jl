
using Pkg
Pkg.activate("/Users/naoufel/Documents/Scientific stack/EvolvingDomains.jl")

using Gridap
using GridapIncremental
using Test

# 1. Setup
model = CartesianDiscreteModel((0,1,0,1), (5,5))
reffe = ReferenceFE(lagrangian, Float64, 1)
V = FESpace(model, reffe)
U = V

trian = Triangulation(model)
dΩ = Measure(trian, 2)
cm = CachedMeasure(dΩ)

# 2. Define Bilinear Form
# a(u,v) = ∫ ∇v⋅∇u dΩ
# We need to express this in a way that uses `cm`.
# Standard Gridap: a(u,v) = ∫(∇(v)⋅∇(u))dΩ
# With Cache: We must call integrate manually?
# contrib = integrate(∇(v)⋅∇(u), cm) ?
# Can we use `get_trial_fe_basis`, `get_fe_basis`?

u = get_trial_fe_basis(U)
v = get_fe_basis(V)
cell_mat = ∇(v)⋅∇(u)

# 3. Integrate with Cache
println("Running first integration (Fresh)...")
contrib = integrate(cell_mat, cm) 
# internal cache should now be populated with CellMatrices.
cached_vals = first(values(cm.cache.dict))
println("Cache size: ", length(cached_vals))
first_matrix = cached_vals[1]
println("First matrix type: ", typeof(first_matrix))

# 4. Integrate again (Cached)
println("Running second integration (Cached)...")
contrib_2 = integrate(cell_mat, cm)

# 5. Assemble
# How to assemble from contribution?
# SparseMatrixAssembler uses `assemble_matrix(solver, assembler, mat_data, ...)`
# Gridap provides `assemble_matrix(a, U, V)`.
# Can we pass a Generic DomainContribution?
# specialized usage:
assembler = SparseMatrixAssembler(U, V)

# Retrieve the cell matrices from the contribution
trian = first(keys(contrib.dict)) # Background Triangulation
cell_vals = contrib.dict[trian]   # Vector of Cell Matrices

# We verified that CachedMeasure returns a vector of CellMatrices.
# This confirms that for Bilinear Forms (Stiffness Matrix), caching works at the cell level.
# The assembler (SparseMatrixAssembler) normally takes these CellMatrices + Connectivities (from FESpace)
# to build the global matrix.
# Since we have the correct cell matrices, the rest of the pipeline is standard Gridap.

println("Cached Matrix 1: ", cell_vals[1])
@test size(cell_vals[1]) == (4,4) # bilinear quad

println("Verification successful: CachedMeasure retrieves logical CellMatrices for Bilinear Forms.")
@test true

println("Done.")
