using Gridap: FESpace, FEFunction, CellField, interpolate, Point, get_triangulation, get_free_dof_values
using ..Geometric: CartesianGridInfo, CartesianMeshField, get_interpolator
import TransferOperator
import Gridap

"""
    GridMeshTransfer

Operator to transfer fields between a Cartesian background grid and a 
finite element space (typically on a CutFEM mesh).
"""
struct GridMeshTransfer
    grid_info::CartesianGridInfo
    target_space::FESpace
    active_indices::Vector{Int}
end

"""
    restrict(op::GridMeshTransfer, u_grid::CartesianMeshField)

Map data from Eulerian Grid to FE Mesh (FESpace).
Leverages the grid's interpolator to provide a continuous CellField for Gridap.
"""
function TransferOperator.restrict(op::GridMeshTransfer, u_grid::CartesianMeshField)
    itp = get_interpolator(u_grid)
    # Wrap the interpolator in a CellField that Gridap can sample at quadrature points
    trian = get_triangulation(op.target_space)
    cf = CellField(x -> itp(x[1], x[2]), trian)
    
    # Project/Interpolate the CellField into the FE space
    return interpolate(cf, op.target_space)
end

"""
    prolong(op::GridMeshTransfer, u_mesh::FEFunction)

Map data from FE Mesh (FEFunction) to Eulerian Grid.
Optimized to use direct index mapping (Glue) if available, falling back to 
batch evaluation otherwise.
"""
function TransferOperator.prolong(op::GridMeshTransfer, u_mesh::FEFunction)
    nx, ny = op.grid_info.dims
    
    # 1. Attempt High-Efficiency Path (Direct Index Mapping)
    # This avoids all geometric searches and point evaluations.
    trian = get_triangulation(u_mesh)
    
    # BodyFittedTriangulation from a GridPortion contains node_to_parent_node
    # which maps: active_node_idx -> background_node_idx
    if hasproperty(trian, :grid) && hasproperty(trian.grid, :node_to_parent_node)
        mapping = trian.grid.node_to_parent_node
        
        # Apply constraints if in AgFEMSpace to get actual nodal values
        V = u_mesh.fe_space
        u_vals = if hasproperty(V, :space)
            Vstd = V.space
            u_std = interpolate(u_mesh, Vstd)
            get_free_dof_values(u_std)
        else
            get_free_dof_values(u_mesh)
        end
        
        # Determine if we need to reinterpret (for VectorValue fields)
        # We compare the number of DOFs to the number of nodes.
        n_nodes = Gridap.Geometry.num_nodes(trian)
        n_dofs = length(u_vals)
        
        T_dof = eltype(u_vals)
        
        if n_dofs == n_nodes
            # Scalar case
            new_data = zeros(T_dof, nx * ny)
            new_data[mapping] .= u_vals
            return CartesianMeshField(new_data, op.grid_info)
        elseif n_dofs > n_nodes && n_dofs % n_nodes == 0
            # Vector case: DOFs are interleaved [u1, v1, u2, v2...]
            # We determine the dimension from the ratio
            dim = n_dofs ÷ n_nodes
            
            # Since we only support 2D currently, we can assume VectorValue{dim, T_dof}
            # To be robust, we create the type dynamically.
            T_vec = Gridap.TensorValues.VectorValue{dim, T_dof}
            u_vals_view = reinterpret(T_vec, u_vals)
            
            new_data = zeros(T_vec, nx * ny)
            new_data[mapping] .= u_vals_view
            return CartesianMeshField(new_data, op.grid_info)
        else
            # Unexpected DOF layout (e.g. Higher order or complex constraints)
            # Fallback to Batch Evaluation
            return _prolong_batch(op, u_mesh)
        end
    end

    return _prolong_batch(op, u_mesh)
end

function _prolong_batch(op::GridMeshTransfer, u_mesh::FEFunction)
    nx, ny = op.grid_info.dims
    origin = op.grid_info.origin
    dx, dy = op.grid_info.spacing
    
    if isempty(op.active_indices)
         error("Cannot prolong field: No active indices found in geometry.")
    end

    # Pre-calculate points to avoid overhead in evaluation
    points = Vector{Point{2,Float64}}(undef, length(op.active_indices))
    for (k, idx) in enumerate(op.active_indices)
        i = (idx - 1) % nx + 1
        j = (idx - 1) ÷ nx + 1
        points[k] = Point(origin[1] + (i-1)*dx, origin[2] + (j-1)*dy)
    end

    # Batch evaluate on the FEFunction
    vals = u_mesh(points)
    
    T = eltype(vals)
    new_data = zeros(T, nx * ny)
    new_data[op.active_indices] .= vals
    
    return CartesianMeshField(new_data, op.grid_info)
end
