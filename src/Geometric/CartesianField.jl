module CartesianField

using ..Geometric: CartesianGridInfo

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
    ClampedBilinear{T}

Callable bilinear interpolant over a *snapshot* of a [`CartesianMeshField`](@ref)'s data.

For a point `(x, y)` let `sx = clamp((x - x0)/dx, 0, nx-1)` and
`sy = clamp((y - y0)/dy, 0, ny-1)` be the fractional grid indices, and set
`i = floor(sx)`, `j = floor(sy)`, `fx = sx - i`, `fy = sy - j`. The value is the
tensor-product blend of the four surrounding nodes with weights
`(1-fx)(1-fy)`, `fx(1-fy)`, `(1-fx)fy`, `fx·fy`. Nodes are read from the flat
vector with x-fast indexing `i + j*nx + 1` (0-based `i`, `j`).

Boundary rules: the clamp is `Flat()` extrapolation, so coordinates at or outside
the grid return the nearest edge/corner value and `±Inf` clamps to an edge; `NaN`
coordinates return `NaN`. The data is copied at construction, so later in-place
writes to `field.data` are not observed.
"""
struct ClampedBilinear{T}
    data::Vector{T}
    nx::Int
    ny::Int
    x0::Float64
    y0::Float64
    dx::Float64
    dy::Float64
end

@inline function (itp::ClampedBilinear{T})(x::Real, y::Real) where {T}
    # `Flat()` extrapolation propagates NaN through the index arithmetic; match it.
    if isnan(x) || isnan(y)
        return itp.data[1] * NaN
    end

    # Fractional grid indices (0-based). `clamp` also implements the `Flat()`
    # boundary rule: it turns ±Inf into an edge, and the boundary value is then
    # read exactly because the clamped fraction is 0 or 1.
    sx = clamp((x - itp.x0) / itp.dx, 0.0, itp.nx - 1)
    sy = clamp((y - itp.y0) / itp.dy, 0.0, itp.ny - 1)

    # Trust boundary: after the clamp 0 <= sx <= nx-1 (same for y), so `floor`
    # lands in [0, n-1] and the upper neighbour is clamped into the same range.
    # All four flat reads below are therefore in 1:length(data); `@inbounds` is safe.
    i0 = floor(Int, sx); i1 = min(i0 + 1, itp.nx - 1)
    j0 = floor(Int, sy); j1 = min(j0 + 1, itp.ny - 1)
    fx = sx - i0; fy = sy - j0

    nx = itp.nx
    k00 = i0 + j0 * nx + 1
    k01 = i0 + j1 * nx + 1
    k10 = i1 + j0 * nx + 1
    k11 = i1 + j1 * nx + 1
    @inbounds begin
        v00 = itp.data[k00]
        v01 = itp.data[k01]
        v10 = itp.data[k10]
        v11 = itp.data[k11]
    end

    w00 = (1 - fx) * (1 - fy)
    w10 = fx * (1 - fy)
    w01 = (1 - fx) * fy
    w11 = fx * fy
    return v00 * w00 + v10 * w10 + v01 * w01 + v11 * w11
end

"""
    get_interpolator(f::CartesianMeshField)

Returns a callable [`ClampedBilinear`](@ref) supporting evaluation at arbitrary
coordinates `itp(x, y)`. Uses linear interpolation with flat extrapolation
(clamping) at boundaries.
"""
function get_interpolator(f::CartesianMeshField)
    nx, ny = f.grid.dims
    (nx >= 1 && ny >= 1) || throw(ArgumentError(
        "CartesianMeshField grid dimensions must be positive, got ($nx, $ny)."))
    # Interpolations.jl's BSpline(Linear()) rejected singleton dimensions.
    (nx >= 2 && ny >= 2) || throw(ArgumentError(
        "CartesianMeshField size ($nx, $ny) is inconsistent with linear interpolation; " *
        "each dimension needs at least 2 nodes."))

    expected = Base.checked_mul(nx, ny)
    length(f.data) == expected || throw(DimensionMismatch(
        "CartesianMeshField has $(length(f.data)) values; grid requires $expected."))

    x0, y0 = f.grid.origin
    dx, dy = f.grid.spacing
    isfinite(x0) || throw(ArgumentError(
        "CartesianMeshField x origin must be finite, got $x0."))
    isfinite(y0) || throw(ArgumentError(
        "CartesianMeshField y origin must be finite, got $y0."))
    # Match Interpolations.jl `scale`, which rejects non-positive/non-finite steps.
    (isfinite(dx) && dx > 0) || throw(ArgumentError(
        "CartesianMeshField x spacing must be finite and positive, got $dx."))
    (isfinite(dy) && dy > 0) || throw(ArgumentError(
        "CartesianMeshField y spacing must be finite and positive, got $dy."))
    # The far edge must be representable, otherwise the index map is not finite.
    isfinite(x0 + (nx - 1) * dx) || throw(ArgumentError(
        "CartesianMeshField x upper bound x0 + (nx-1)*dx is not representable."))
    isfinite(y0 + (ny - 1) * dy) || throw(ArgumentError(
        "CartesianMeshField y upper bound y0 + (ny-1)*dy is not representable."))

    # Snapshot the data, matching the old interpolator (Interpolations.jl copied
    # the coefficient array): later in-place writes to `f.data` are not observed.
    return ClampedBilinear{eltype(f.data)}(copy(f.data), nx, ny, x0, y0, dx, dy)
end

end # module
