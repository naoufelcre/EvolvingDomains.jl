using Gridap: FESpace, FEFunction, CellField, interpolate, Point, get_triangulation
using ..Geometric: CartesianGridInfo, CartesianMeshField, get_interpolator
import TransferOperator

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
Evaluates the FEFunction at grid nodes corresponding to active_indices.
"""
function TransferOperator.prolong(op::GridMeshTransfer, u_mesh::FEFunction)
    nx, ny = op.grid_info.dims
    origin = op.grid_info.origin
    dx, dy = op.grid_info.spacing
    
    # Determine result type by evaluating at a representative point.
    # We MUST use an active point, otherwise Gridap throws an error 
    # if the point is outside the active triangulation.
    if isempty(op.active_indices)
         error("Cannot prolong field: No active indices found in geometry.")
    end
    
    first_idx = op.active_indices[1]
    i_0 = (first_idx - 1) % nx + 1
    j_0 = (first_idx - 1) ÷ nx + 1
    x_sample = Point(origin[1] + (i_0-1)*dx, origin[2] + (j_0-1)*dy)

    T = typeof(u_mesh(x_sample))
    new_data = zeros(T, nx * ny)
    
    # Evaluate FEFunction at the physical coordinates of the active nodes
    for idx in op.active_indices
        # Convert linear index to i, j (1-based)
        i = (idx - 1) % nx + 1
        j = (idx - 1) ÷ nx + 1
        
        # Physical coordinates of the grid node
        x_phys = Point(origin[1] + (i-1)*dx, origin[2] + (j-1)*dy)
        
        # Evaluate the FE function at this point
        new_data[idx] = u_mesh(x_phys)
    end
    
    return CartesianMeshField(new_data, op.grid_info)
end
