# This file contains mandatory information for integration with external solvers.

"""
    _GridMetrics

Metadata about the background Cartesian grid for external solver integration.

# Fields
- `origin`: Grid origin
- `spacing`: Cell sizes
- `dims`: Node counts

Cell counts are `dims .- 1` and are derived at the few call sites that need them.
"""
struct _GridMetrics{N}
    origin::NTuple{N,Float64}
    spacing::NTuple{N,Float64}
    dims::NTuple{N,Int}
end

function _GridMetrics(origin::NTuple{N,<:Real}, spacing::NTuple{N,<:Real},
                           dims::NTuple{N,<:Integer}) where N
    N in (2, 3) || throw(ArgumentError("Cartesian grids support only 2D or 3D."))
    return _GridMetrics{N}(Float64.(origin), Float64.(spacing), Int.(dims))
end

"""
    _checkmetrics(info) -> info

Validate a `_GridMetrics` that a caller built by hand. Metrics produced by `_gridMetrics`
already satisfy every clause, so kernels call this only on caller-supplied metrics.
"""
function _checkmetrics(info::_GridMetrics)
    all(>=(2), info.dims) ||
        throw(ArgumentError("grid axes require at least two nodes."))
    all(isfinite, info.origin) || throw(ArgumentError("grid origins must be finite."))
    all(h -> isfinite(h) && h > 0, info.spacing) ||
        throw(ArgumentError("grid spacings must be finite and positive."))
    all(eachindex(info.dims)) do d
        n, origin, h = info.dims[d], info.origin[d], info.spacing[d]
        upper = origin + (n - 1) * h
        isfinite(upper) && origin + h > origin && upper > origin + (n - 2) * h
    end || throw(ArgumentError("grid nodes must be finite and distinctly representable."))
    return info
end

"""
    _gridMetrics(grid::Grid) -> _GridMetrics

Derive uniform-grid metrics from a package `Grid`. Throws when the grid is not 2D or 3D or
when an axis is not uniformly spaced. This is the only metric builder in the package: CIP,
WENO5, cut geometry and transport all measure a grid the same way.
"""
function _gridMetrics(grid::Grid)
    axes = grid.axes
    N = length(axes)
    N in (2, 3) || throw(ArgumentError("numerical operators support only 2D or 3D grids."))
    dims = ntuple(d -> length(axes[d]), N)
    origin = ntuple(d -> Float64(first(axes[d])), N)
    spacing = ntuple(d -> Float64((last(axes[d]) - first(axes[d])) / (dims[d] - 1)), N)
    for d in 1:N
        axis, h = axes[d], spacing[d]
        all(isfinite, axis) && isfinite(h) && h > 0 ||
            throw(ArgumentError("grid axes must contain finite, increasing coordinates."))
        # Tolerance scaled by the axis magnitude, not by the local coordinate: the rounding
        # error in a node built as `first + (i-1)h` accumulates from the axis magnitude, so
        # a tolerance that shrinks near a zero-crossing node rejects an axis that is exactly
        # uniform (e.g. range(-2.0, 2.0; length = 321)).
        x1 = Float64(axis[1])
        tolerance = _coordinateTolerance(x1, axis[end], h)
        for i in eachindex(axis)
            abs(Float64(axis[i]) - (x1 + (i - 1) * h)) <= tolerance ||
                throw(ArgumentError("numerical operators support only uniformly spaced grids."))
        end
    end
    return _checkmetrics(_GridMetrics(origin, spacing, dims))
end

# Mesh identity of a Grid or a Cut. Both carry `mesh_id`; nothing else is accepted.
@inline sameid(a, b) = getproperty(a, :mesh_id) === getproperty(b, :mesh_id)

@inline samemetrics(a::_GridMetrics, b::_GridMetrics) =
    a.origin == b.origin && a.spacing == b.spacing && a.dims == b.dims

function _cutGridInfo(cut, grid)
    sameid(cut, grid) || throw(ArgumentError(
        "Cut and Field Grid mesh identities differ."))
    return _gridMetrics(grid)
end

"""
    _checklevel(source, level_id, info, label) -> data

Check the dense level array belonging to the grid. `source` is a `Field` with a level ID,
or the array itself when `level_id === nothing`. Container exactness is not required here;
kernels that index a contiguous buffer assert `Array{Float64,N}` at their own call site.
"""
function _checklevel(source, level_id, info::_GridMetrics, label::AbstractString)
    data = if level_id === nothing
        source
    else
        arrays = getproperty(source, :data)
        haskey(arrays, level_id) ||
            throw(ArgumentError("$label has no values attached to level $(repr(level_id))."))
        arrays[level_id]
    end
    data isa AbstractArray{Float64} && ndims(data) == length(info.dims) || throw(ArgumentError(
        "$label must be a rank-$(length(info.dims)) Float64 array."))
    size(data) == info.dims || throw(DimensionMismatch(
        "$label has size $(size(data)); grid requires $(info.dims)."))
    all(isfinite, data) || throw(ArgumentError(
        "$label contains a non-finite value at Cartesian index $(findfirst(!isfinite, data))."))
    return data
end

function _cutFieldInfo(cut, grid, data)
    info = _cutGridInfo(cut, grid)
    _checklevel(data, nothing, info, "Level-set Field values")
    return info
end
