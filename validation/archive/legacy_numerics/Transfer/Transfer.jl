module Transfer

using ..Geometric: CartesianGridInfo, get_interpolator
import TransferOperator
using Gridap: FESpace

include("ExtensionOperator.jl")
include("GridMeshTransfer.jl")

export GridMeshTransfer, ClosestPointExtension, extend

"""Apply a closest-point extension operator to a flat grid buffer."""
extend(op::ClosestPointExtension, data::AbstractVector) = extend_field(op, data)

end # module
