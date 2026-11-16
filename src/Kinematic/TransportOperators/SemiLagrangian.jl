module SemiLagrangian

using StaticArrays

using ...Geometric: _GridMetrics, _getInterpolator, _gridMetrics, _sameid, _samemetrics
using ...EvolvingDomains: Cut, Field, Grid

@inline function _quadraticInterpolationWeights(α, stencil::Symbol)
    if stencil == :left
        return 0.5α*(α-1), 1-α^2, 0.5α*(α+1)
    elseif stencil == :right
        return 0.5*(α-1)*(α-2), α*(2-α), 0.5α*(α-1)
    end
    throw(ArgumentError("Unknown stencil: $stencil"))
end

# The backward pull integrates over the whole neighbourhood of each ray destination by
# supersampling. Without it, diverging characteristics destroy mass in the pull phase and the
# forward push has to invent it back, which produces unphysical results.

"""
    TransportMap(previous_cut, current_cut, grid, velocity, dt; cache_weights=false)

Build a reusable conservative transport map between two cut snapshots, for geometry-aware
Lagrangian transport. It redistributes only from the nodes incident to in-phase or cut cells, and
conserves the unweighted nodal sum only for contributions that reach the target support.

Algorithm: supersampled backward pull with demand scaling and forward leakage redistribution.
Reference: Lentine, Grétarsson & Fedkiw (2011), doi:10.1016/j.jcp.2010.12.036
"""
struct TransportMap{N,Q}
    # All fields below are construction snapshots: rays, masks, offsets and the
    # packed stencils are read during lagrangianTransport! but are not designed to be mutated.
    # (source_mask is the exception: lagrangianTransport! detects a changed mask and falls back.)
    #
    # --- Backward (Pull Phase) ---
    # target_active_idx[k] -> bundle of 4 (2D) or 8 (3D) departure points
    active_indices::Vector{Int}
    backward_rays::Vector{SVector{Q,SVector{N,Float64}}}

    # --- Conservation ---
    # source_idx -> total_weight_pulled (Used to scale weights for conservation)
    demand_map::Vector{Float64}

    # --- Optional packed pull-weight cache ---
    # _computeConservativeWeights{2,3} is step-invariant:
    # it depends only on the cached backward rays, the grid and the source
    # support, never on the field or dt. With `cache_weights=true` its nonzero
    # outputs are stored flat in ray order, then x-fast stencil order:
    #
    #   wq = 2^-N * w_m  (stencil weight with the supersample share folded in)
    #   weight_offsets[k] <= p < weight_offsets[k+1] selects node k's entries.
    #
    # Dropping zero weights and keeping their original order makes the packed sum
    # reproduce the uncached loop bit for bit.
    weight_offsets::Vector{Int}      # empty exactly when no packed weights were stored
    weight_indices::Vector{Int32}    # linear source-node index per entry
    weight_values::Vector{Float64}   # wq = 0.25 * w_m per entry
    # source_mask snapshot; lagrangianTransport! ignores the cache when the live mask differs,
    # so editing the support never silently reuses stale weights.
    # Also empty when no cache was stored. Whether the cache exists is read from
    # `isempty(weight_offsets)`; storing a second flag allowed the two to disagree.
    cached_source_mask::BitVector

    # --- Forward (Push Phase / Leakage Correction) ---
    # source_idx -> arrival_point (Only for nodes not fully resolved by Pull phase)
    leakage_indices::Vector{Int}
    leakage_rays::Vector{SVector{N,Float64}}

    # --- Metadata ---
    grid_meta::_GridMetrics{N}
    # Source and target supports are the discrete geometry at t⁻ and t⁺.
    # The level set itself is deliberately not retained in the map.
    source_mask::BitVector
    target_mask::BitVector
    is_identity::Bool
    # Identity of the mesh the map was built for. Concrete: every builder path sets them,
    # so a `nothing` state cannot arise and is not representable.
    grid_id::UInt64
    level_id::Symbol
end

function TransportMap(
    previous_cut::Cut{P,N},
    current_cut::Cut{Q,N},
    grid::Grid{D},
    velocity::AbstractVector{<:SVector{N}},
    dt::Real;
    cache_weights::Bool=false,
) where {P,Q,N,D}
    meta, coords, active_previous, active_current =
        _cutTransportInputs(previous_cut, current_cut, grid)
    Base.require_one_based_indexing(velocity)
    expected = prod(meta.dims)
    length(velocity) == expected || throw(DimensionMismatch(
        "TransportMap velocity has $(length(velocity)) values; grid requires $expected."))
    all(v -> all(isfinite, v), velocity) ||
        throw(ArgumentError("TransportMap velocity contains non-finite values."))
    sampled = velocity isa Vector{SVector{N,Float64}} ? velocity :
        [SVector(ntuple(d -> Float64(v[d]), N)) for v in velocity]
    sampled === velocity || all(v -> all(isfinite, v), sampled) || throw(ArgumentError(
        "TransportMap velocity contains non-finite Float64 values."))
    interpolator = _getInterpolator(sampled, meta)
    return _buildTransportMap(meta, coords, active_previous, active_current,
        x -> interpolator(x...), dt; cache_weights,
        grid_id=previous_cut.mesh_id, level_id=previous_cut.level_id)
end

function _cutTransportInputs(previous_cut::Cut{P,N}, current_cut::Cut{Q,N},
    grid::Grid{D}) where {P,Q,N,D}
    N == D || throw(ArgumentError(
        "TransportMap cuts are $N-dimensional but the grid is $D-dimensional."))
    _sameid(previous_cut, current_cut) || throw(ArgumentError(
        "TransportMap previous and current cuts belong to different grids."))
    previous_cut.level_id == current_cut.level_id || throw(ArgumentError(
        "TransportMap previous and current cuts belong to different levels."))
    _sameid(current_cut, grid) || throw(ArgumentError(
        "TransportMap cuts and grid belong to different meshes."))
    meta = _gridMetrics(grid)
    # Linear node order: first axis fastest, as in the nested per-dimension comprehensions.
    coords = Vector{SVector{N,Float64}}(undef, prod(meta.dims))
    for (node_index, node) in enumerate(Iterators.product(grid.axes...))
        coords[node_index] = SVector(ntuple(dim -> Float64(node[dim]), N))
    end
    return meta, coords, _cutActiveNodes(previous_cut, meta.dims),
        _cutActiveNodes(current_cut, meta.dims)
end

# The nodal band a transport step may read: the corners of every in-phase or cut cell.
#
# shortcut: this band is one cell wide while the interpolation stencil spans four nodes, so a
# target at the band edge has stencil nodes outside the band. The product clip below (and not
# the per-axis clip the comment describes) is what keeps such a target from having no support
# at all. Widening the band does not fix it: the target set grows with it. The real remedy is
# the ghost-cell layer of Lentine et al., which is not implemented.
function _cutActiveNodes(cut::Cut{T,N}, dims::NTuple{N,Int}) where {T,N}
    cell_dims = ntuple(dim -> dims[dim] - 1, N)
    cell_index, node_index = LinearIndices(cell_dims), LinearIndices(dims)
    active = Int[]
    for cell in CartesianIndices(cell_dims)
        state = cut.bg_states[cell_index[cell]]
        (state == -1 || state == 0) || continue
        for corner in CartesianIndices(ntuple(_ -> 0:1, N))
            push!(active, node_index[cell + corner])
        end
    end
    sort!(unique!(active))
    return active
end

"""
    _buildTransportMap(meta, coords, active_previous, active_current, velocity_at, dt;
                       cache_weights, grid_id, level_id)

Trace the backward rays of every active target node with `2^N` supersample offsets, record the
conservation demand they create, optionally pack the pull stencils, and collect the forward
leakage rays for mass the pull phase leaves behind. The map follows Lentine, Grétarsson &
Fedkiw (2011), doi:10.1016/j.jcp.2010.12.036.
"""
function _buildTransportMap(meta::_GridMetrics{N}, coords, active_previous,
    active_current, velocity_at, dt::Real; cache_weights::Bool=false,
    grid_id::UInt64, level_id::Symbol) where {N}
    isfinite(dt) || throw(ArgumentError("TransportMap time step must be finite."))
    dt >= 0 || throw(ArgumentError("TransportMap time step must be non-negative."))

    n_nodes = foldl(Base.checked_mul, meta.dims)
    source_mask = falses(n_nodes)
    target_mask = falses(n_nodes)
    source_mask[active_previous] .= true
    target_mask[active_current] .= true
    any(source_mask) && !any(target_mask) && throw(ArgumentError(
        "TransportMap cannot conserve mass from a nonempty source mask to an empty target mask."))

    rays_per_node = 2^N
    # Supersample offsets: one per corner of the ±h/4 box, first axis fastest.
    offsets = ntuple(q -> SVector(ntuple(
            d -> ((q - 1) >> (d - 1)) & 1 == 1 ? 0.25 * meta.spacing[d] : -0.25 * meta.spacing[d],
            N)), rays_per_node)
    backward_rays = Vector{SVector{rays_per_node,SVector{N,Float64}}}(undef, length(active_current))
    stationary = true
    for (k, i) in enumerate(active_current)
        rays = MVector{rays_per_node,SVector{N,Float64}}(undef)
        for q in 1:rays_per_node
            departure = coords[i] + offsets[q]
            rays[q] = _traceRay(departure, velocity_at, -dt)
            stationary &= rays[q] == departure
        end
        backward_rays[k] = SVector(rays)
    end

    # 2. Conservation demand, plus the optional packed pull weights. The weights are only
    # materialized when requested (default off) and when the node index fits in Int32, so
    # the default map allocates no large cache arrays.
    weights_requested = cache_weights && n_nodes <= typemax(Int32)
    stencil = 2^(2N)
    share = 2.0^-N                      # 0.25 in 2D, 0.125 in 3D; exact in binary
    demand = zeros(Float64, n_nodes)
    weight_offsets = weights_requested ? Vector{Int}(undef, length(backward_rays) + 1) : Int[]
    weight_indices, weight_values = Int32[], Float64[]
    if weights_requested
        weight_offsets[1] = 1
        sizehint!(weight_indices, stencil * length(backward_rays))
        sizehint!(weight_values, stencil * length(backward_rays))
    end
    for k in eachindex(backward_rays)
        weights_requested && (weight_offsets[k] = length(weight_indices) + 1)
        for departure in backward_rays[k]
            indices, weights = _computeConservativeWeights(departure, meta, source_mask)
            # A stencil with no visible node contributes nothing. That is the correct open
            # boundary when the ray leaves the box (there is no source data outside), and a
            # map defect when the departure point is inside the domain.
            if sum(weights) == 0.0 && _insideDomain(departure, meta)
                node = Tuple(coords[active_current[k]])
                throw(ArgumentError(
                    "TransportMap: the backward ray from node $(active_current[k]) at " *
                    "$(node) to $(Tuple(departure)) reaches no active source node, so the " *
                    "map cannot be built for this step. This happens when the step is long " *
                    "enough that the band no longer covers the characteristic, or when the " *
                    "band is thinner than the stencil. Choosing the timestep is the " *
                    "caller's responsibility; see the CFL note in README.md."))
            end
            for m in 1:stencil
                s_idx = indices[m]
                wq = share * weights[m]
                source_mask[s_idx] && (demand[s_idx] += wq)
                if weights_requested && wq > 0.0
                    push!(weight_indices, Int32(s_idx))
                    push!(weight_values, wq)
                end
            end
        end
        weights_requested && (weight_offsets[k + 1] = length(weight_indices) + 1)
    end

    if stationary
        for i in eachindex(source_mask)
            source_mask[i] && (stationary &= _traceRay(coords[i], velocity_at, dt) == coords[i])
        end
    end
    stationary &= source_mask == target_mask

    # 3. Leakage rays (forward) for mass the pull phase does not collect.
    leak_idx, leak_rays = Int[], SVector{N,Float64}[]
    for i in 1:n_nodes
        source_mask[i] || continue
        if demand[i] < 1.0
            arrival = _traceRay(coords[i], velocity_at, dt)
            _, forward_weights = _computeConservativeWeights(arrival, meta, target_mask)
            # Same split as the pull side: a ray that leaves the box carries its leftover mass
            # out of the domain (open outflow), while an inside ray with no visible target is
            # a support hole the map cannot conserve through.
            if sum(forward_weights) == 0.0 && _insideDomain(arrival, meta)
                throw(ArgumentError(
                    "TransportMap: the forward ray from source node $(i) at " *
                    "$(Tuple(coords[i])) to $(Tuple(arrival)) reaches no active target node, " *
                    "so the leftover mass of this node has nowhere to go and the map cannot " *
                    "conserve it. Same two causes as for a pull stencil: a step long enough " *
                    "to leave the band, or a band thinner than the stencil."))
            end
            push!(leak_idx, i)
            push!(leak_rays, arrival)
        end
    end

    return TransportMap(active_current, backward_rays, demand,
        weight_offsets, weight_indices, weight_values,
        weights_requested ? copy(source_mask) : BitVector(),
        leak_idx, leak_rays, meta, source_mask, target_mask, stationary,
        grid_id, level_id)
end

"""
    lagrangianTransport!(target_data, source_data, map)
    lagrangianTransport!(target::Field, source::Field, map)

Using a transport map, do Lagrangian transport of the nodal values, mutating the target in
place. Source and target must not alias each other or the map, and entries outside the target
support are set to zero.

Algorithm: supersampled pull with per-source demand scaling and forward leakage push.
Reference: Lentine, Grétarsson & Fedkiw (2011), doi:10.1016/j.jcp.2010.12.036
"""
function lagrangianTransport!(target_data::Vector{Float64}, source_data::Vector{Float64},
    map::TransportMap{N,Q}) where {N,Q}
    Q == 2^N || throw(ArgumentError(
        "TransportMap{$N,$Q} is inconsistent: a $N-dimensional map traces 2^$N rays per node."))
    candidates = 2^(2N)
    share = 2.0^-N              # exact in binary: 0.25 in 2D, 0.125 in 3D
    expected = length(map.source_mask)
    length(target_data) == expected || throw(DimensionMismatch(
        "lagrangianTransport!: target has $(length(target_data)) values; map requires $expected."))
    length(source_data) == expected || throw(DimensionMismatch(
        "lagrangianTransport!: source has $(length(source_data)) values; map requires $expected."))
    Base.mightalias(target_data, source_data) && throw(ArgumentError(
        "Conservative lagrangianTransport! source and target must not alias."))
    (Base.mightalias(target_data, map.demand_map) ||
     Base.mightalias(source_data, map.demand_map) ||
     Base.mightalias(target_data, map.weight_values) ||
     Base.mightalias(source_data, map.weight_values)) && throw(ArgumentError(
        "Conservative lagrangianTransport! fields must not alias TransportMap storage."))
    for i in eachindex(source_data)
        map.source_mask[i] && !isfinite(source_data[i]) && throw(ArgumentError(
            "Conservative lagrangianTransport! source contains non-finite values on its support."))
    end

    if map.is_identity
        fill!(target_data, 0.0)
        target_data[map.active_indices] .= source_data[map.active_indices]
        return target_data
    end

    fill!(target_data, 0.0)

    # --- Pass 1: Backward Pull (with 2^N supersampling) ---
    # The packed stencils are used only when the map carries a valid cache and
    # the live support still matches its construction snapshot; otherwise the
    # stencils are recomputed from the cached rays and the live mask.
    use_cache = !isempty(map.weight_offsets) && map.source_mask == map.cached_source_mask
    if use_cache
        Base.Threads.@threads for k in eachindex(map.active_indices)
            target_idx = map.active_indices[k]
            value = 0.0
            @inbounds for p in map.weight_offsets[k]:(map.weight_offsets[k + 1] - 1)
                s_idx = Int(map.weight_indices[p])
                demand = map.demand_map[s_idx]
                # Scale by 1/demand on the source node if over-requested (mass
                # conservation); the packed weight already carries the 2^-N share.
                scale = demand > 1.0 ? (1.0 / demand) : 1.0
                value += map.weight_values[p] * scale * source_data[s_idx]
            end
            target_data[target_idx] = value
        end
    else
        Base.Threads.@threads for k in eachindex(map.active_indices)
            target_idx = map.active_indices[k]
            value = 0.0
            for departure in map.backward_rays[k]
                indices, weights = _computeConservativeWeights(departure, map.grid_meta, map.source_mask)
                for m in 1:candidates
                    s_idx = indices[m]
                    weight = weights[m]
                    if weight > 0.0
                        demand = map.demand_map[s_idx]
                        scale = demand > 1.0 ? (1.0 / demand) : 1.0
                        value += share * weight * scale * source_data[s_idx]
                    end
                end
            end
            target_data[target_idx] = value
        end
    end

    # --- Pass 2: Forward Push (Leakage Correction) ---
    for (k, s_idx) in enumerate(map.leakage_indices)
        value = source_data[s_idx]
        if value != 0.0
            leftover = (1.0 - map.demand_map[s_idx]) * value
            indices, weights = _computeConservativeWeights(
                map.leakage_rays[k], map.grid_meta, map.target_mask)
            for m in 1:candidates
                target_data[indices[m]] += weights[m] * leftover
            end
        end
    end

    return target_data
end

function lagrangianTransport!(target::Field{T,N}, source::Field{S,N},
    map::TransportMap{N,Q}) where {T,S,N,Q}
    td, sd = target.discretization, source.discretization
    (haskey(td.levels, map.level_id) && haskey(sd.levels, map.level_id) &&
     haskey(target.data, map.level_id) && haskey(source.data, map.level_id)) ||
        throw(ArgumentError("Conservative lagrangianTransport! Fields must contain map level $(map.level_id)."))
    tg, sg = td.levels[map.level_id], sd.levels[map.level_id]
    tg.mesh_id === map.grid_id || throw(ArgumentError(
        "Conservative lagrangianTransport! target belongs to a different grid identity."))
    sg.mesh_id === map.grid_id || throw(ArgumentError(
        "Conservative lagrangianTransport! source belongs to a different grid identity."))
    _samemetrics(_gridMetrics(sg), map.grid_meta) || throw(ArgumentError(
        "Conservative lagrangianTransport! grid coordinates changed after map construction."))
    target_data, source_data = target.data[map.level_id], source.data[map.level_id]
    target_data isa Array{Float64,N} || throw(ArgumentError(
        "Conservative lagrangianTransport! target Field values must be a $(N)D Float64 array."))
    source_data isa Array{Float64,N} || throw(ArgumentError(
        "Conservative lagrangianTransport! source Field values must be a $(N)D Float64 array."))
    dims = ntuple(d -> length(sg.axes[d]), N)
    size(source_data) == dims || throw(DimensionMismatch(
        "Conservative lagrangianTransport! source Field has size $(size(source_data)); grid requires $dims."))
    size(target_data) == dims || throw(DimensionMismatch(
        "Conservative lagrangianTransport! target Field has size $(size(target_data)); grid requires $dims."))
    lagrangianTransport!(vec(target_data), vec(source_data), map)
    return target
end

# Utilities

"""    _traceRay(x, velocity_at, dt)

Trace a characteristic with the explicit trapezoidal rule.
"""
function _traceRay(x::SVector{D,T}, velocity_at, dt) where {D,T}
    raw_v1 = velocity_at(x)
    length(raw_v1) == D || throw(ArgumentError(
        "Velocity dimension does not match the $D-dimensional transport grid."))
    v1 = SVector(ntuple(d -> Float64(raw_v1[d]), D)...)
    all(isfinite, v1) || throw(ArgumentError(
        "Non-finite velocity encountered while tracing from $x."))

    x_mid = x + v1 * dt
    raw_v2 = velocity_at(x_mid)
    length(raw_v2) == D || throw(ArgumentError(
        "Velocity dimension does not match the $D-dimensional transport grid."))
    v2 = SVector(ntuple(d -> Float64(raw_v2[d]), D)...)
    all(isfinite, v2) || throw(ArgumentError(
        "Non-finite velocity encountered while tracing through $x_mid."))

    x_new = x + 0.5 * (v1 + v2) * dt
    return x_new
end

# Is the point inside the closed box? Used to tell "the characteristic left the domain"
# (contributes nothing, mass leaves) from "the stencil found no support" (a map defect).
@inline function _insideDomain(x::SVector{N}, meta::_GridMetrics{N}) where {N}
    for d in 1:N
        lower = meta.origin[d]
        upper = lower + (meta.dims[d] - 1) * meta.spacing[d]
        (lower <= x[d] <= upper) || return false
    end
    return true
end

# One merged per-axis kernel: the average of the left- and right-anchored quadratic Lagrange
# stencils, `(α(α-1)/4, 1-(3α+α²)/4, (5α-α²)/4, α(α-1)/4)` on nodes `i-1 .. i+2`. It sums to
# 1, its mean offset is α, and it is exact for quadratics; the two outer lobes are negative.
@inline function _mergedQuadraticWeights(α)
    left = _quadraticInterpolationWeights(α, :left)
    right = _quadraticInterpolationWeights(α, :right)
    return (0.5 * left[1], 0.5 * (left[2] + right[1]), 0.5 * (left[3] + right[2]),
            0.5 * right[3])
end

# The `4^N` stencil of the piecewise-quadratic reconstruction: one body for 2D and 3D.
#
# Lentine, Grétarsson & Fedkiw (2011): a stencil point outside the domain wall or outside the
# allowed support is "not visible", takes weight 0, and the visible weights are rescaled to sum
# 1. The sample point itself is never relocated, so a ray outside the box stays invisible
# instead of folding back onto the boundary ring.
#
# The quadratic kernel's two outer lobes are negative and are clipped, which is what makes
# `Σ w = 1` a positive partition. The clip is applied to the weight *product*, so a corner whose
# two lobes are both negative survives with positive weight: the applied kernel is then not a
# tensor product and its mean offset is α + 0.023 (measured), not α.
#
# shortcut: clipping the lobes one axis at a time is the formulation the comment above
# describes and removes that 0.023 bias (measured: the box boundary layer then thins from five
# nodes to two), but it empties the stencil of every target at the band edge and the map can
# only refuse those steps. It needs the ghost layer first; see `_cutActiveNodes`.
@inline function _stencilWeights!(indices, weights, x::SVector{N,T}, grid::_GridMetrics{N},
                                  allowed) where {N,T}
    strides = ntuple(d -> d == 1 ? 1 : prod(grid.dims[1:d-1]), N)
    raw = ntuple(d -> (x[d] - grid.origin[d]) / grid.spacing[d] + 1, N)
    if !all(isfinite, raw)
        # A non-finite departure point has no support at all; giving it the first node would
        # silently weigh against the wrong corner instead of raising the caller's band error.
        fill!(weights, 0.0)
        fill!(indices, 1)
        return nothing
    end
    position = ntuple(d -> clamp(Float64(raw[d]), -1.0e9, 1.0e9), N)
    base = ntuple(d -> floor(Int, position[d]), N)
    kernel = ntuple(d -> _mergedQuadraticWeights(position[d] - base[d]), N)

    total = 0.0
    q = 1
    for corner in CartesianIndices(ntuple(_ -> 1:4, N))
        linear = 1
        visible = true
        lobe = 1.0
        for d in 1:N
            node = base[d] + corner[d] - 2
            visible &= 1 <= node <= grid.dims[d]
            linear += (clamp(node, 1, grid.dims[d]) - 1) * strides[d]
            lobe *= kernel[d][corner[d]]
        end
        indices[q] = linear
        weight = (visible && (allowed === nothing || allowed[linear])) ? max(0.0, lobe) : 0.0
        weights[q] = weight
        total += weight
        q += 1
    end

    if total > 0.0
        for k in eachindex(weights)
            weights[k] /= total
        end
    else
        fill!(weights, 0.0)
    end
    return nothing
end

"""    _computeConservativeWeights(x, grid, allowed=nothing)

Construct the conservative quadratic interpolation stencil at `x`, used by the pull and leakage
phases of a `TransportMap`.
"""
@inline function _computeConservativeWeights(x::SVector{N,T}, grid::_GridMetrics{N},
                                            allowed=nothing) where {N,T}
    indices = MVector{4^N,Int}(undef)
    weights = MVector{4^N,Float64}(undef)
    _stencilWeights!(indices, weights, x, grid, allowed)
    return SVector(indices), SVector(weights)
end

end # module
