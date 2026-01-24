# =============================================================================
# GridInfo.jl — External Solver Integration
# =============================================================================

using Gridap.Geometry: get_cartesian_descriptor, CartesianDiscreteModel

"""
    CartesianGridInfo

Metadata about the background Cartesian grid for external solver integration.

# Fields
- `origin`: Grid origin (x₀, y₀)
- `spacing`: Cell sizes (Δx, Δy)
- `dims`: Node counts (nx, ny)
- `cells`: Cell counts (nx-1, ny-1)
"""
struct CartesianGridInfo
    origin::NTuple{2,Float64}
    spacing::NTuple{2,Float64}
    dims::NTuple{2,Int}
    cells::NTuple{2,Int}
end

"""
    grid_info(model::CartesianDiscreteModel) -> CartesianGridInfo

Extract grid metadata from a Gridap CartesianDiscreteModel.
"""
function grid_info(model::CartesianDiscreteModel)
    desc = get_cartesian_descriptor(model)
    origin = Tuple(desc.origin)
    spacing = Tuple(desc.sizes)  # desc.sizes is already cell sizes
    partition = Tuple(desc.partition)
    dims = partition .+ 1
    return CartesianGridInfo(origin, spacing, dims, partition)
end



