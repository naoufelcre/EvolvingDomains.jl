# =============================================================================
# GridMeshTransfer.jl — Built-in Transfer Operator for Grid ↔ Mesh
# =============================================================================

import TransferOperator
using TransferOperator: AbstractTransferOperator
using Gridap.Geometry: get_node_coordinates
using Gridap.Arrays: evaluate
using Gridap.TensorValues: VectorValue

export GridMeshTransfer, setup_transfer

"""
    GridMeshTransfer{M} <: AbstractTransferOperator

Default transfer operator between a Cartesian grid and an arbitrary mesh.
Uses batch evaluation for prolong, bilinear interpolation for restrict.
"""
struct GridMeshTransfer{M} <: AbstractTransferOperator
    origin::NTuple{2,Float64}
    spacing::NTuple{2,Float64}
    dims::NTuple{2,Int}
    mesh::M
end

"""
    setup_transfer(grid_info, mesh) -> GridMeshTransfer

Construct a transfer operator from grid info and a mesh reference.
"""
function setup_transfer(info::CartesianGridInfo, mesh)
    return GridMeshTransfer(info.origin, info.spacing, info.dims, mesh)
end

# Convenience: from CartesianDiscreteModel
function setup_transfer(model::CartesianDiscreteModel, mesh)
    return setup_transfer(grid_info(model), mesh)
end

"""
    prolong(op::GridMeshTransfer, mesh_field) -> Array{T,2}

Transfer mesh field to grid via batch evaluation.
Void points are filled with `NaN`.
"""
function TransferOperator.prolong(op::GridMeshTransfer, mesh_field)
    nx, ny = op.dims

    # 1. Generate all grid node coordinates
    coords = Vector{VectorValue{2,Float64}}(undef, nx * ny)
    for j in 1:ny, i in 1:nx
        x = op.origin[1] + (i-1) * op.spacing[1]
        y = op.origin[2] + (j-1) * op.spacing[2]
        coords[i + (j-1)*nx] = VectorValue(x, y)
    end

    # 2. Manual Evaluation Loop
    # Detect type from first evaluation
    val0 = mesh_field(coords[1])
    T = typeof(val0)
    result = Matrix{T}(undef, nx, ny)

    for k in 1:length(coords)
        idx_i = (k - 1) % nx + 1
        idx_j = div(k - 1, nx) + 1
        result[idx_i, idx_j] = mesh_field(coords[k])
    end

    return result
end

"""
    restrict(op::GridMeshTransfer, grid_field) -> Function

Create interpolator from grid field for mesh evaluation.
"""
function TransferOperator.restrict(op::GridMeshTransfer, grid_field::AbstractArray{T,2}) where T
    return x -> begin
        return bilinear_interpolate_scalar(
            vec(grid_field), op.origin, op.spacing, op.dims, (x[1], x[2])
        )
    end
end

"""
    locate(op::GridMeshTransfer, x) -> Union{Int, Nothing}

Find grid cell containing point x.
"""
function TransferOperator.locate(op::GridMeshTransfer, x)
    i = floor(Int, (x[1] - op.origin[1]) / op.spacing[1]) + 1
    j = floor(Int, (x[2] - op.origin[2]) / op.spacing[2]) + 1
    nx, ny = op.dims
    (1 <= i <= nx && 1 <= j <= ny) || return nothing

    return i + (j-1) * nx
end
