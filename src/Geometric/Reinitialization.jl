module Reinitialization

using ..Geometric: _gridMetrics, _GridMetrics, _cutFieldInfo
using ....EvolvingDomains: Cut, Field, CUT_CUT, _levelgrid, _levelvalues

# The `N`-dimensional index helpers every operator in this file uses.
@inline _strides(dims::NTuple{N,Int}) where {N} =
    ntuple(d -> d == 1 ? 1 : prod(dims[1:d-1]), N)

@inline _linear(index, strides::NTuple{N,Int}) where {N} =
    1 + sum((index[d] - 1) * strides[d] for d in 1:N)

@inline _coordinate(linear::Int, axis::Int, dims, strides) =
    (linear - 1) ÷ strides[axis] % dims[axis] + 1

"""
    _levelsetScale(phi, info, cut_nodes) -> Float64

One scale converting level-set units to length, taken as the median |∇φ| over the cut
nodes.

Deliberately a **single number** shared by every cut node, not a per-node value. The seed
below divides by it, and the crossing fraction only survives if both nodes of an edge
divide by the *same* thing — a per-node |∇φ| reintroduces exactly the inconsistency this
is built to remove. The median rather than the mean because a few cut nodes sit where the
central difference straddles a kink and reports a wild gradient.
"""
function _levelsetScale(phi::Vector{Float64}, info::_GridMetrics{N},
                        cut_nodes::Vector{Int}) where {N}
    dims, spacing = info.dims, info.spacing
    strides = _strides(dims)
    g = Float64[]
    sizehint!(g, length(cut_nodes))
    for idx in cut_nodes
        # Central difference of φ along each axis, clamped into the box at the edge nodes.
        gradient = ntuple(N) do d
            c = _coordinate(idx, d, dims, strides)
            lo, hi = max(c - 1, 1), min(c + 1, dims[d])
            (phi[idx + (hi - c) * strides[d]] - phi[idx + (lo - c) * strides[d]]) /
                ((hi - lo) * spacing[d])
        end
        m = hypot(gradient...)
        m > 0 && push!(g, m)
    end
    isempty(g) && throw(ArgumentError(
        "reinitialize!: the level set has zero gradient at every cut node, so it carries no " *
        "length scale to convert to."))
    return sort!(g)[(length(g) + 1) ÷ 2]
end

"""
    reinitialize!(levelset::Field, cut::Cut)

Reinitialize the Cut's level-set values by fast sweeping, seeding nodes of explicit CUT
cells so the represented zero crossings remain fixed.

**One pass, so it restores a usable level set, not an exact distance function.** Signs and
zero crossings are preserved, but a curved front is not a fixed point: the sweeping carries
first-order truncation error, of order one cell at 20–25 cells from the interface and
growing with distance (0.8h at 81² measured, and a disc of radius 0.3 reads −0.27 at its
centre on 17²). An axis-aligned planar front *is* a fixed point. That is harmless for
sign-based work and for cut geometry, and it is not harmless for anything that consumes φ as
a distance (moment-fitted quadratures, φ-valued lifting, extrapolation): iterate or solve a
transport equation for those. Used in a loop it can also amplify a small cut-volume bias
(measured ~10× more drift with reinitialization than without over 100 steps).
"""
function reinitialize!(levelset::Field, cut::Cut)
    level_id = cut.level_id
    d = levelset.discretization
    g = _levelgrid(d, level_id)
    phi = _levelvalues(levelset, level_id)
    info = _cutFieldInfo(cut, g, phi)
    dims = info.dims
    node_dim = length(dims)
    strides = _strides(dims)
    cell_dims = ntuple(k -> dims[k] - 1, node_dim)
    cell_strides = _strides(cell_dims)

    # Seed nodes: every vertex of every background cell the cut marks as straddling.
    cut_nodes = Int[]
    corners = CartesianIndices(ntuple(_ -> 0:1, node_dim))
    for cell_idx in eachindex(cut.bg_states)
        cut.bg_states[cell_idx] == CUT_CUT || continue
        base = ntuple(k -> (cell_idx - 1) ÷ cell_strides[k] % cell_dims[k] + 1, node_dim)
        for corner in corners
            push!(cut_nodes, _linear(ntuple(k -> base[k] + corner[k], node_dim), strides))
        end
    end
    unique!(cut_nodes)
    if isempty(cut_nodes)
        @warn "reinitialize!: level set is single-signed (no interface on this grid); " *
              "leaving it unchanged." maxlog = 1
        return levelset
    end

    original = vec(phi)
    t = fill(1e10, length(original))
    frozen = fill(false, length(original))
    scale = _levelsetScale(original, info, cut_nodes)
    for idx in cut_nodes
        t[idx] = abs(original[idx]) / scale
        frozen[idx] = true
    end
    _fastSweeping!(t, info, frozen)
    for i in eachindex(original)
        original[i] = sign(original[i]) * t[i]
    end
    return levelset
end

"""
    _fastSweeping!(t, info, frozen) -> t

Solve `|∇T| = 1` in place by fast sweeping: one Gauss–Seidel pass per monotone direction
(`2^N` of them), each node taking the backtracking neighbours of its own sweep direction,
which is the upwind term for that quadrant. `frozen` marks the seeded Dirichlet nodes.
"""
function _fastSweeping!(t::Vector{Float64}, info::_GridMetrics{N}, frozen) where {N}
    dims, spacing = info.dims, info.spacing
    strides = _strides(dims)
    outer = CartesianIndices(ntuple(_ -> 1:2, N))
    for direction in outer
        ranges = ntuple(d -> direction[d] == 1 ? (1:dims[d]) : (dims[d]:-1:1), N)
        for I in CartesianIndices(ranges)
            idx = _linear(I, strides)
            frozen[idx] && continue
            neighbours = ntuple(N) do d
                from = I[d] + (direction[d] == 1 ? -1 : 1)
                1 <= from <= dims[d] ? (t[idx + (from - I[d]) * strides[d]], spacing[d]) :
                    (1e10, spacing[d])
            end
            t[idx] = min(t[idx], _solveEikonal(sort(collect(neighbours); by=first)))
        end
    end
    return t
end

# `Iterators.product` with a single range yields scalars, so iterate a one-tuple of the value.
"""
    _solveEikonal(vals) -> Float64

Upwind eikonal update: `vals` holds `(neighbour_time, spacing)` pairs sorted by time.
Start from the 1-D update and add neighbours while the solution stays above the last one
included — the standard fast-sweeping "drop the far neighbours" reduction. One
implementation serves both dimensions.
"""
@inline function _solveEikonal(vals)
    for n in 2:length(vals)
        u = vals[n][1]
        solution = _solveEikonalSubset(vals, n)
        solution <= u && return _solveEikonalSubset(vals, n - 1)
    end
    return _solveEikonalSubset(vals, length(vals))
end

@inline function _solveEikonalSubset(vals, n)
    A = B = C = 0.0
    for q in 1:n
        u, h = vals[q]
        w = 1 / (h * h); A += w; B -= 2u * w; C += u * u * w
    end
    max((-B + sqrt(max(B * B - 4A * (C - 1), 0.0))) / (2A), vals[n][1])
end

end # module
