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
