module Transfer

using ..Geometric: EvolvingDiscreteGeometry, current_levelset, grid_info, CartesianGridInfo
using ..Geometric: GeometryCache, CartesianMeshField, get_interpolator, get_active_indices
import TransferOperator
using Gridap: CartesianDiscreteModel, FESpace, FEFunction, get_triangulation

include("ExtensionOperator.jl")
include("GridMeshTransfer.jl")

export get_transfer_op, prolong, restrict, get_extension_op, extend
export update_transfer_cache!, update_extension_cache!
export GridMeshTransfer, setup_transfer
export grid_to_mesh, mesh_to_grid

# --- Transfer Logic ---

"""
    setup_transfer(geom::EvolvingDiscreteGeometry, target_space::FESpace)

Initialize the GridMeshTransfer operator for a given FE space and current geometry.
"""
function setup_transfer(geom::EvolvingDiscreteGeometry, target_space::FESpace)
    info = grid_info(geom.grid)
    active_indices = get_active_indices(geom, :current)
    op = GridMeshTransfer(info, target_space, active_indices)
    update_transfer_cache!(geom, op)
    return op
end

function get_transfer_op(geom::EvolvingDiscreteGeometry)
    if isnothing(geom.cache.transfer_op)
        error("Transfer operator not initialized. Call setup_transfer(geom, target_space) first.")
    end
    return geom.cache.transfer_op
end

"""
    grid_to_mesh(geom::EvolvingDiscreteGeometry, u_grid)
    mesh_to_grid(geom::EvolvingDiscreteGeometry, u_mesh)

High-level aliases for spatial transfer between Eulerian background grid 
and Lagrangian/Cut-cell mesh. 
`grid_to_mesh` maps from Grid to Mesh (Restriction).
`mesh_to_grid` maps from Mesh to Grid (Prolongation).
"""
function grid_to_mesh(geom::EvolvingDiscreteGeometry, u_grid)
    op = get_transfer_op(geom)
    return TransferOperator.restrict(op, u_grid)
end

function mesh_to_grid(geom::EvolvingDiscreteGeometry, u_mesh)
    op = get_transfer_op(geom)
    return TransferOperator.prolong(op, u_mesh)
end

# --- Interface compliance with TransferOperator ---

# We extend TransferOperator methods so that TransferOperator.restrict(geom, ...) works.
# These simply delegate to the grid_to_mesh/mesh_to_grid functions above.

function TransferOperator.restrict(geom::EvolvingDiscreteGeometry, u_grid)
    return grid_to_mesh(geom, u_grid)
end

function TransferOperator.prolong(geom::EvolvingDiscreteGeometry, u_mesh)
    return mesh_to_grid(geom, u_mesh)
end

function get_extension_op(geom::EvolvingDiscreteGeometry)

    if isnothing(geom.cache.extension_op)

        # ϕ_vals = current_levelset(geom) # Not implemented yet
        ϕ_vals = geom.levelset

        info = grid_info(geom.grid)

        new_op = ClosestPointExtension(info, ϕ_vals)

        geom.cache.extension_op = new_op

    end

    return geom.cache.extension_op

end

function extend(geom::EvolvingDiscreteGeometry, u_grid)
    op = get_extension_op(geom)
    return extend_field(op, u_grid)
end

function update_transfer_cache!(geom::EvolvingDiscreteGeometry, op)
    geom.cache.transfer_op = op
    return geom
end

function update_extension_cache!(geom::EvolvingDiscreteGeometry, op)
    geom.cache.extension_op = op
    return geom
end

end # module
