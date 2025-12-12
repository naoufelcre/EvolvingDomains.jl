# =============================================================================
# Grid Information for External Solvers
# =============================================================================

"""
    CartesianGridInfo{N}

Immutable struct containing Cartesian grid metadata for external solvers.

# Fields
- `origin`: Grid origin coordinates (NTuple{N,Float64})
- `spacing`: Cell sizes Δx, Δy, ... (NTuple{N,Float64}) 
- `dims`: Node counts nx, ny, ... (NTuple{N,Int})
- `cells`: Cell counts (dims .- 1) (NTuple{N,Int})

# Example
```julia
info = CartesianGridInfo(model)
Δx, Δy = info.spacing
nx, ny = info.dims
```
"""
struct CartesianGridInfo{N}
    origin::NTuple{N,Float64}
    spacing::NTuple{N,Float64}
    dims::NTuple{N,Int}
    cells::NTuple{N,Int}
end

"""
    CartesianGridInfo(model::CartesianDiscreteModel)

Construct grid info from a Gridap CartesianDiscreteModel.
"""
function CartesianGridInfo(model::CartesianDiscreteModel)
    (origin, corner, partition) = cartesian_descriptor(model)
    dims = partition .+ 1
    spacing = (corner .- origin) ./ partition
    return CartesianGridInfo(origin, spacing, dims, partition)
end

# Convenience methods
Base.ndims(::CartesianGridInfo{N}) where {N} = N
