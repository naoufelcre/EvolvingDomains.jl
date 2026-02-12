module CartesianField

using ..Geometric: CartesianGridInfo
using Interpolations

export CartesianMeshField
export meshsize
export get_interpolator

#For defining fields over the domains and to use it conveniently between the mesh and the grid representation

"""
    CartesianMeshField{T}

Wrapper around a flat data vector and grid info to support Cartesian indexing
and boundary handling (clamping).
"""
struct CartesianMeshField{T}
    data::Vector{T}
    grid::CartesianGridInfo
end

@inline function Base.getindex(f::CartesianMeshField, I::CartesianIndex{2})
    nx, ny = f.grid.dims
    i = clamp(I[1], 1, nx)
    j = clamp(I[2], 1, ny)
    @inbounds f.data[i+(j-1)*nx]
end

@inline function Base.setindex!(f::CartesianMeshField, v, I::CartesianIndex{2})
    nx, ny = f.grid.dims
    if 1 <= I[1] <= nx && 1 <= I[2] <= ny
        @inbounds f.data[I[1]+(I[2]-1)*nx] = v
    else
        error("Index out of bounds")
    end
end

@inline function meshsize(f::CartesianMeshField, dim::Int)
    @inbounds f.grid.spacing[dim]
end

"""
    get_interpolator(f::CartesianMeshField)

Returns an object that supports evaluation at arbitrary coordinates `itp(x, y)`.
Uses linear interpolation with flat extrapolation (clamping) at boundaries.
"""
function get_interpolator(f::CartesianMeshField)
    # Create a 2D view of the flat data without copying
    nx, ny = f.grid.dims
    data_2d = reshape(f.data, nx, ny)

    # Define physical ranges
    # Note: Interpolations.jl expects ranges to match dimensions
    x0, y0 = f.grid.origin
    dx, dy = f.grid.spacing

    # Create ranges (using standard StepRangeLen for efficiency)
    xaxis = range(x0, step=dx, length=nx)
    yaxis = range(y0, step=dy, length=ny)

    # 1. Interpolate on the integer grid (BSpline Linear)
    itp = interpolate(data_2d, BSpline(Linear()))

    # 2. Scale to physical coordinates
    sitp = scale(itp, xaxis, yaxis)

    # 3. Extrapolate with Flat (clamping) to match getindex behavior
    eitp = extrapolate(sitp, Flat())

    return eitp
end

end # module
