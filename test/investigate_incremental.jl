
using Pkg
Pkg.activate("/Users/naoufel/Documents/Scientific stack/EvolvingDomains.jl")
Pkg.develop(path="/Users/naoufel/Documents/Scientific stack/GridapIncremental.jl")

using Gridap
using GridapEmbedded
using EvolvingDomains
using GridapIncremental
using Test

# 1. Setup Initial Geometry
model = CartesianDiscreteModel((0,1,0,1), (10,10))
initial_ls(x) = sqrt((x[1]-0.5)^2 + (x[2]-0.5)^2) - 0.3
geom = EvolvingDiscreteGeometry(model, initial_ls)

# Initial Cut
# Initial Cut
cut_1 = current_cut(geom)
println("typeof(cut_1): ", typeof(cut_1))
# Try to access underlying geometry
# In GridapEmbedded <= 0.9, cut_1 might be Generalized/Standard/ActiveBody EmbeddedDiscretization
# They typically store `geo` or `ls_to_bgcell_to_inoutcut`

function get_cell_status(cut_dist, model)
    # Attempt to extract status array
    # If cut_dist has a geometry field
    if hasfield(typeof(cut_dist), :geo)
        geo = cut_dist.geo
        if hasfield(typeof(geo), :bgcell_to_inoutcut)
            return geo.bgcell_to_inoutcut
        end
        # Look for other fields
         println("Geo fields: ", fieldnames(typeof(geo)))
    end
    
    println("CutDist fields: ", fieldnames(typeof(cut_dist)))
    return nothing
end

status_1 = get_cell_status(cut_1, model)
println("Step 0: Status captured.")

# 2. Advance Geometry
advance!(geom, 0.1) 
# We need to manually move the LS because velocity is 0 by default and we didn't set one?
# Lets just set a new LS manually to simulate movement
new_ls(x) = sqrt((x[1]-0.55)^2 + (x[2]-0.55)^2) - 0.3
set_levelset!(geom, collect(new_ls(x) for x in get_node_coordinates(geom.model)))

cut_2 = current_cut(geom)
status_2 = get_cell_status(cut_2, model)

# 3. Compute Difference
diff_mask = status_1 .!= status_2
dirty_cells = findall(diff_mask)

println("Dirty cells count: ", length(dirty_cells))
println("Dirty cells: ", dirty_cells)

# 4. Mock Update CachedMeasure
# cm = CachedMeasure(...)
# cm.dirty_cells = dirty_cells
# integrate(..., cm)

