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
    
    # Lower corner indices (1-based)
    # Clamp to ensure we are within valid cell range [1, nx-1] x [1, ny-1]
    i = clamp(floor(Int, fx) + 1, 1, nx - 1)
    j = clamp(floor(Int, fy) + 1, 1, ny - 1)
    
    # Local coordinates in cell [0, 1]
    # We re-calculate based on clamped i,j to handle extrapolation consistently
    # or just clamp the local coords?
    # Logic from VelocitySource:
    # tx = clamp(fx - (i - 1), 0.0, 1.0)
    # Logic from EvolvingGeometry:
    # ξ = (px - (ox + (i-1)*dx)) / dx
    
    # Let's use the explicit local coord calculation which is safer for extrapolation
    tx = (x[1] - (ox + (i-1)*dx)) / dx
    ty = (x[2] - (oy + (j-1)*dy)) / dy
    
    # If point is outside, linear extrapolation is risky for Level Sets?
    # EvolvingGeometry logic: "If point is outside, we extrapolate from nearest cell"
    # VelocitySource logic: "clamp tx/ty to 0..1" -> This effectively does nearest-neighbor extrapolation relative to the boundary cell's field
    # Let's stick to Clamped Local Coords for safety (Nearest Neighbor of border value effectively if far out)
    # But wait, true bilinear extrapolation might be better for SDFs? 
    # Actually, clamped is safer to avoid wild values. `VelocitySource` does clamping. `prior EvolvingGeometry` did not explicitly clamp local coords but relied on clamped indices.
    # Let's use clamping for `tx` `ty` to 0.0-1.0 to stay within the cell values.
    tx = clamp(tx, 0.0, 1.0)
    ty = clamp(ty, 0.0, 1.0)

    # 1D Linear indices (Column-major)
    # idx(i, j) = i + (j-1)*nx
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
