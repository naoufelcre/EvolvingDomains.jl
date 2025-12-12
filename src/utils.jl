# =============================================================================
# Utility Functions
# =============================================================================

"""
    get_node_coords_as_vector(model::CartesianDiscreteModel) -> Vector{VectorValue}

Extract node coordinates from a CartesianDiscreteModel as a flat vector.
"""
function get_node_coords_as_vector(model::CartesianDiscreteModel)
    trian = Triangulation(model)
    return Gridap.Geometry.get_node_coordinates(trian)
end

"""
    cartesian_descriptor(model::CartesianDiscreteModel) -> (origin, corner, partition)

Extract the origin, corner, and partition from a CartesianDiscreteModel.
"""
function cartesian_descriptor(model::CartesianDiscreteModel)
    # Access internal CartesianDescriptor
    # This is version-dependent, so we provide a safe fallback
    try
        desc = model.grid_topology.grid.desc
        origin = desc.origin
        sizes = desc.sizes
        partition = desc.partition
        corner = origin .+ sizes
        return (Tuple(origin), Tuple(corner), Tuple(partition))
    catch
        # Fallback: infer from node coordinates
        coords = get_node_coords_as_vector(model)
        xs = [c[1] for c in coords]
        ys = [c[2] for c in coords]
        x_unique = sort(unique(xs))
        y_unique = sort(unique(ys))
        origin = (x_unique[1], y_unique[1])
        corner = (x_unique[end], y_unique[end])
        partition = (length(x_unique) - 1, length(y_unique) - 1)
        return (origin, corner, partition)
    end
end

# =============================================================================
# Domain Mask Utilities (for external solvers)
# =============================================================================

"""
    domain_mask(ϕ::Vector{Float64}) -> BitVector

Returns a mask where `true` indicates inside the domain (ϕ < 0).

# Example
```julia
ϕ = current_levelset(eg)
mask = domain_mask(ϕ)
inside_nodes = findall(mask)
```
"""
domain_mask(ϕ::AbstractVector{<:Real}) = ϕ .< 0

"""
    narrow_band_mask(ϕ::Vector{Float64}, bandwidth::Real) -> BitVector

Returns a mask for nodes within `bandwidth` of the interface (|ϕ| < bandwidth).
Useful for narrow-band level set methods.

# Arguments
- `ϕ`: Level set values at grid nodes
- `bandwidth`: Width of the narrow band (typically k × Δx for stencil width k)

# Example
```julia
ϕ = current_levelset(eg)
info = grid_info(eg)
γ = 6 * info.spacing[1]  # 6 cells for WENO5
band = narrow_band_mask(ϕ, γ)
```
"""
narrow_band_mask(ϕ::AbstractVector{<:Real}, bandwidth::Real) = abs.(ϕ) .< bandwidth
