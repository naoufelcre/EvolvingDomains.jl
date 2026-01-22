module InterpolationUtils

using StaticArrays

export bilinear_interpolate_scalar
export smootherstep

"""
    smootherstep(x)

C² smooth transition function:
- 1.0 for x <= 0
- 0.0 for x >= 1
- Smooth blend in between
"""
@inline function smootherstep(x::Real)
    x ≤ 0 && return 1.0
    x ≥ 1 && return 0.0
    return 1.0 - x^3 * (x * (6x - 15) + 10)
end

"""
    bilinear_interpolate_scalar(values, origin, spacing, dims, x)

Bilinear interpolation of a scalar field defined on a Cartesian grid.

# Arguments
- `values`: Scalar values stored in a flat vector (column-major order)
- `origin`: Tuple (x0, y0) of the grid origin
- `spacing`: Tuple (dx, dy) of the grid spacing
- `dims`: Tuple (nx, ny) of the grid dimensions (number of nodes)
- `x`: Tuple or Vector (px, py) of the query point

Returns extrapolated values (nearest cell) if outside the grid (clamped indices).
"""
@inline function bilinear_interpolate_scalar(values::AbstractVector, origin, spacing, dims, x)
    ox, oy = origin
    dx, dy = spacing
    nx, ny = dims

    # Grid-relative coordinates
    fx = (x[1] - ox) / dx
    fy = (x[2] - oy) / dy

    i = clamp(floor(Int, fx) + 1, 1, nx - 1)
    j = clamp(floor(Int, fy) + 1, 1, ny - 1)

    tx = (x[1] - (ox + (i-1)*dx)) / dx
    ty = (x[2] - (oy + (j-1)*dy)) / dy

    tx = clamp(tx, 0.0, 1.0)
    ty = clamp(ty, 0.0, 1.0)

    idx00 = i + (j-1)*nx
    idx10 = (i+1) + (j-1)*nx
    idx01 = i + j*nx
    idx11 = (i+1) + j*nx

    # Bilinear interpolation
    val = (1-tx)*(1-ty)*values[idx00] +
          tx*(1-ty)*values[idx10] +
          (1-tx)*ty*values[idx01] +
          tx*ty*values[idx11]

    return val
end

end # module
