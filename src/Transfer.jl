# Multilevel transfer: interpolate the source level, evaluate at the target level's nodes.
#
# The whole algebra is one verb (`transfer!`) and one scheme protocol. A scheme is any
# callable `itp(x, values, grid) -> Float64`, so a user plugs in their own function or a
# closure over per-grid precomputation without touching the package.

"""
    linearInterpolator(x, values, grid)

A linear interpolator, so the library comes batteries included. Returns the interpolated
value of a nodal field at the point `x`.
"""
function linearInterpolator(x::NTuple{N,Float64}, values::AbstractArray{<:Real,N},
                            grid::Grid) where {N}
    length(grid.axes) == N || throw(DimensionMismatch(
        "point is $N-dimensional but the grid is $(length(grid.axes))-dimensional"))
    size(values) == _nodeDims(grid) || throw(DimensionMismatch(
        "values have size $(size(values)); grid requires $(_nodeDims(grid))"))
    # Per axis: lower node and the fraction towards the next one.
    lower = ntuple(N) do d
        axis = grid.axes[d]
        i = clamp(searchsortedlast(axis, x[d]), 1, length(axis) - 1)
        width = axis[i + 1] - axis[i]
        (i, width == 0 ? 0.0 : clamp((x[d] - axis[i]) / width, 0.0, 1.0))
    end
    value = 0.0
    for corner in 0:(2^N - 1)
        weight = 1.0
        index = ntuple(N) do d
            upper = (corner >> (d - 1)) & 1 == 1
            weight *= upper ? lower[d][2] : 1.0 - lower[d][2]
            lower[d][1] + (upper ? 1 : 0)
        end
        value += weight * values[index...]
    end
    return value
end

# True when one grid's axes are a subsampling of the other's on every axis, which is what
# refine!/coarsen! produce. Transfer itself would interpolate between any two grids over the
# shared domain; the check exists so that moving data between levels always means "nested
# transfer", never an unannounced resampling.
#
# The nesting relation is exactly `x_longer(1 + (i-1)k) = x_shorter(i)` for one integer stride
# `k` per axis: the divisibility of `n - 1`, then the node test, and nothing else. Coordinates
# are compared with the package tolerance rather than with `==`, because the same set of nodes
# written two ways (refine! by a factor of 3 versus 49 equally spaced points) differs by an ulp;
# a genuinely different mesh misses by orders of magnitude more.
function _nested(a::Grid, b::Grid)
    length(a.axes) == length(b.axes) || return false
    for k in eachindex(a.axes)
        x, y = a.axes[k], b.axes[k]
        longer, shorter = length(x) >= length(y) ? (x, y) : (y, x)
        (length(longer) - 1) % (length(shorter) - 1) == 0 || return false
        stride = (length(longer) - 1) ÷ (length(shorter) - 1)
        coarse = (Float64(last(longer)) - Float64(first(longer))) / (length(shorter) - 1)
        tol = _coordinateTolerance(first(longer), last(longer), coarse)
        all(i -> abs(Float64(longer[1 + stride * (i - 1)]) - Float64(shorter[i])) <= tol,
            eachindex(shorter)) || return false
    end
    return true
end

"""
    transfer!(field, from::Symbol, to::Symbol; itp=linearInterpolator)

Transfer the nodal values stored at level `from` to level `to` using an interpolator.
`itp(x, values, grid)` receives a tuple of physical coordinates, and fine-to-coarse transfer is
point sampling, so it is not conservative.
"""
function transfer!(field::Field{T,N}, from::Symbol, to::Symbol;
                   itp=linearInterpolator) where {T,N}
    from === to && throw(ArgumentError("transfer! needs two different levels"))
    d = getproperty(field, :discretization)
    values = getproperty(field, :data)
    source_grid, target_grid = _levelgrid(d, from), _levelgrid(d, to)
    _levelvalues(field, from)
    _levelvalues(field, to)
    _nested(source_grid, target_grid) || throw(ArgumentError(
        "levels $from and $to are not nested: neither axis set is a subsampling of the " *
        "other. Build them with refine!/coarsen!, or interpolate by hand."))
    target = values[to]
    size(target) == _nodeDims(target_grid) || throw(DimensionMismatch(
        "level $to holds values of size $(size(target)); its grid requires $(_nodeDims(target_grid))"))
    source = values[from]
    for index in CartesianIndices(target)
        x = ntuple(k -> Float64(target_grid.axes[k][index[k]]), N)
        target[index] = itp(x, source, source_grid)
    end
    return field
end
