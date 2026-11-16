module GridInterpolation

using ..Geometric: _GridMetrics, _checkmetrics

# Internal flat-grid interpolation helpers. The kernels pass their value arrays and a
# _GridMetrics explicitly; there is no separate user-facing field container.

"""Internal callable bilinear interpolant over copied flat 2D grid data; a queried point outside the grid returns the value of the nearest boundary node."""
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
    if isnan(x) || isnan(y)
        return itp.data[1] * NaN
    end
    sx = clamp((x - itp.x0) / itp.dx, 0.0, itp.nx - 1)
    sy = clamp((y - itp.y0) / itp.dy, 0.0, itp.ny - 1)
    i0 = floor(Int, sx); i1 = min(i0 + 1, itp.nx - 1)
    j0 = floor(Int, sy); j1 = min(j0 + 1, itp.ny - 1)
    fx = sx - i0; fy = sy - j0
    nx = itp.nx
    @inbounds begin
        v00 = itp.data[i0 + j0 * nx + 1]
        v01 = itp.data[i0 + j1 * nx + 1]
        v10 = itp.data[i1 + j0 * nx + 1]
        v11 = itp.data[i1 + j1 * nx + 1]
    end
    return (1-fx)*(1-fy)*v00 + fx*(1-fy)*v10 + (1-fx)*fy*v01 + fx*fy*v11
end

"""Internal callable trilinear interpolant over copied flat 3D grid data; a queried point outside the grid returns the value of the nearest boundary node."""
struct ClampedTrilinear{T}
    data::Vector{T}
    dims::NTuple{3,Int}
    origin::NTuple{3,Float64}
    spacing::NTuple{3,Float64}
end

@inline function (itp::ClampedTrilinear)(x::Real, y::Real, z::Real)
    any(isnan, (x, y, z)) && return itp.data[1] * NaN
    dims = itp.dims
    s = ntuple(d -> clamp(((x, y, z)[d] - itp.origin[d]) / itp.spacing[d], 0.0,
                          dims[d] - 1), 3)
    lo = ntuple(d -> floor(Int, s[d]), 3)
    hi = ntuple(d -> min(lo[d] + 1, dims[d] - 1), 3)
    t = ntuple(d -> s[d] - lo[d], 3)
    value = zero(itp.data[1])
    @inbounds for a in 0:1, b in 0:1, c in 0:1
        i = (a == 0 ? lo[1] : hi[1])
        j = (b == 0 ? lo[2] : hi[2])
        k = (c == 0 ? lo[3] : hi[3])
        w = (a == 0 ? 1-t[1] : t[1]) * (b == 0 ? 1-t[2] : t[2]) *
            (c == 0 ? 1-t[3] : t[3])
        value += itp.data[i + j*dims[1] + k*dims[1]*dims[2] + 1] * w
    end
    return value
end

"""    _getInterpolator(data, info)
Return a clamped multilinear interpolator over a copy of flat grid data in first-axis-fastest order; a queried point outside the grid returns the value of the nearest boundary node."""
function _getInterpolator(data::AbstractVector{T}, info::_GridMetrics{N}) where {T,N}
    N in (2, 3) || throw(ArgumentError("flat-grid interpolation supports only 2D or 3D."))
    _checkmetrics(info)
    length(data) == prod(info.dims) || throw(DimensionMismatch(
        "grid data has $(length(data)) values; grid requires $(prod(info.dims))."))
    values = collect(data)
    if N == 2
        return ClampedBilinear{T}(values, info.dims[1], info.dims[2],
            info.origin[1], info.origin[2], info.spacing[1], info.spacing[2])
    end
    return ClampedTrilinear{T}(values, info.dims, info.origin, info.spacing)
end

end # module
