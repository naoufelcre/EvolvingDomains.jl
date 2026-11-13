# Multilevel transfer: interpolate the source level, evaluate at the target level's nodes.
#
# The whole algebra is one verb (`transfer!`) and one scheme protocol. A scheme is any
# callable `itp(x, values, grid) -> Float64`, so a user plugs in their own function or a
# closure over per-grid precomputation without touching the package.

"""
    linearInterpolator(x, values, grid) -> Float64

Multilinear interpolation of `values` (defined on `grid`) at the physical point
`x::NTuple{N,Float64}`. Axes are located with `searchsortedlast`, so nonuniform grids work
and nothing is assumed about uniformity; coordinates outside the grid clamp to the boundary
node.

`values` must be finite. Interpolation is local but not forgiving: every corner of the
enclosing cell is multiplied by its weight, zero weight included (`0 * NaN` is `NaN`), so a
single `NaN` returns `NaN` for that whole cell — even at its finite nodes, where the
offending corner carries no weight — while cells not containing it are unaffected. The
caller owns finiteness, since nothing is scanned or thrown here, matching the rule that the
container stays permissive and each operator validates only when it is used.

This is the default `transfer!` scheme. It is deliberately separate from the transport
kernels' internal snapshot interpolator (`ClampedBilinear`/`ClampedTrilinear`), which is
uniform-grid only, copies its data at construction and exists for per-ray evaluation inside
`lagrangianTransport!`.
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
    transfer!(field, from, to; itp = linearInterpolator) -> field

Fill level `to` from level `from` of the same `Field` by evaluating `itp` at every node of
`to`. Both levels must be attached; the grids must be nested — one axis set a
**constant-stride** subsample of the other, which is what `refine!`/`coarsen!` produce. A
genuine but uneven subset of nodes is refused. Because a `Discretization` enforces one
physical domain, nested levels necessarily share the same box.

Direction is implied by the levels: coarser→finer is prolongation, finer→coarser is
point-sampling (which is **not** conservative — a conserved quantity needs an averaging
restriction, deliberately left to a future `restrict!`).

`itp` is any callable `(x, values, grid) -> Float64`; the default is `linearInterpolator`.

Transfer is interpolation between function spaces, so it cannot create information: a
coarse level's error is reproduced faithfully on the fine level (measured 12% volume error
and 0.59 cells of displacement after two upward transfers from a 9-node start), and no
amount of later refinement removes it. The directions are not symmetric —
coarse→fine→coarse reproduces the coarse values exactly when the two grids coincide node for
node, and only to the nesting tolerance otherwise, while fine→coarse→fine always loses
something — so a
zero round-trip drift in one order is not evidence of a lossless transfer.
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
