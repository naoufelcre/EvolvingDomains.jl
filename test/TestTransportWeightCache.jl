# Focused test for the opt-in packed semi-Lagrangian pull-weight cache.
#
# `TransportMap(geom, velocity, dt; cache_weights=true)` caches, at construction,
# the nonzero conservative stencils of its four backward rays.  `advect!` then
# reads them instead of recomputing `compute_conservative_weights` for every
# active node on every field.  The default (`cache_weights=false`) is the old
# uncached behaviour and allocates no large cache arrays.
#
# This file checks the cached map, the default map and a literal re-implementation
# of the original uncached algorithm (`reference_advect`) all agree bit for bit.
#
# Runnable directly:  julia --project=. test/TestTransportWeightCache.jl
# (Not wired into runtests.jl on purpose — that file belongs to another owner.)

using Test
using EvolvingDomains
using EvolvingDomains.Geometric
using EvolvingDomains.Kinematic
using EvolvingDomains.Kinematic.SemiLagrangian
using Gridap
using Gridap.TensorValues

const SL = EvolvingDomains.Kinematic.SemiLagrangian

# --- Fixtures -----------------------------------------------------------------

function weight_fixture(n)
    model = CartesianDiscreteModel((0.0, 1.0, 0.0, 1.0), (n, n))
    info = grid_info(model)
    points = vec(collect(Gridap.Geometry.get_node_coordinates(model)))
    return model, info, points
end

# A square active region carves a real support boundary (not just the domain
# wall), so some stencils are clipped even though source and target share a mask.
function square_geometry(n)
    model, info, points = weight_fixture(n)
    phi = [max(abs(p[1] - 0.5), abs(p[2] - 0.5)) - 0.3 for p in points]
    geom = EvolvingDiscreteGeometry(vec(phi), model)
    ensure_cut!(geom)
    geom.cache.prev_cut = geom.cache.cut
    return geom, info, points
end

function compressible_map(n; cache_weights::Bool=false)
    geom, info, points = square_geometry(n)
    velocity = [VectorValue(0.3 * (p[1] - 0.5), 0.3 * (p[2] - 0.5)) for p in points]
    return TransportMap(geom, velocity, 0.1; cache_weights=cache_weights), points, info
end

# --- Original uncached algorithm, transcribed from the baseline -----------------
# Kept independent from the cache so a cached/uncached mismatch cannot hide.
function reference_advect(source_data, map)
    target = zeros(Float64, length(source_data))
    if map.is_identity
        target[map.active_indices] .= source_data[map.active_indices]
        return target
    end
    for k in eachindex(map.active_indices)
        val_accum = 0.0
        for x_dep in map.backward_rays[k]
            indices, weights = SL.compute_conservative_weights(
                x_dep, map.grid_meta, map.source_mask)
            for m in 1:16
                s_idx = indices[m]
                w = weights[m]
                if w > 0
                    req = map.demand_map[s_idx]
                    scale = req > 1.0 ? (1.0 / req) : 1.0
                    val_accum += 0.25 * w * scale * source_data[s_idx]
                end
            end
        end
        target[map.active_indices[k]] = val_accum
    end
    for (k, s_idx) in enumerate(map.leakage_indices)
        val = source_data[s_idx]
        if val != 0.0
            req = map.demand_map[s_idx]
            leftover = (1.0 - req) * val
            indices, weights = SL.compute_conservative_weights(
                map.leakage_rays[k], map.grid_meta, map.target_mask)
            for m in 1:16
                target[indices[m]] += weights[m] * leftover
            end
        end
    end
    return target
end

# --- Tests --------------------------------------------------------------------

@testset "Semi-Lagrangian weight cache" begin
    map_on, points, _ = compressible_map(24; cache_weights=true)
    map_off, _, _ = compressible_map(24)          # default = uncached
    n = length(points)

    # The default map carries no cache storage at all.
    @test !map_off.cache_valid
    @test isempty(map_off.weight_offsets)
    @test isempty(map_off.weight_indices)
    @test isempty(map_off.weight_values)
    @test isempty(map_off.cached_source_mask)

    # The opt-in cache is present, packed, and covers exactly the active nodes.
    @test map_on.cache_valid
    @test length(map_on.weight_offsets) == length(map_on.active_indices) + 1
    @test map_on.weight_offsets[1] == 1
    @test map_on.weight_offsets[end] == length(map_on.weight_indices) + 1
    @test length(map_on.weight_indices) == length(map_on.weight_values)
    @test length(map_on.weight_indices) < 64 * length(map_on.active_indices)  # packed, not dense
    @test map_on.cached_source_mask == map_on.source_mask

    # Per-node cached entry counts: the clipped boundary stencils must be
    # strictly sparser than the full-support interior ones.
    node_nnz = [map_on.weight_offsets[k + 1] - map_on.weight_offsets[k]
                for k in 1:length(map_on.active_indices)]
    @test minimum(node_nnz) > 0
    @test minimum(node_nnz) < maximum(node_nnz)   # boundary and interior both present

    # Smooth source supported on the map support.
    source = [map_on.source_mask[i] ? (1.0 + sinpi(p[1]) * cospi(p[2])) : 0.0
              for (i, p) in enumerate(points)]

    # (1) Cached, default and reference agree bit for bit.
    cached = zeros(Float64, n)
    advect!(cached, source, map_on; type=:conservative)
    uncached = zeros(Float64, n)
    advect!(uncached, source, map_off; type=:conservative)
    reference = reference_advect(source, map_on)
    @test cached == reference
    @test uncached == reference
    @test cached == uncached

    # (2) Multiple scalar fields share one map without cross-talk.
    second = [map_on.source_mask[i] ? (0.25 - 2.0 * p[1] + 0.5 * p[2]) : 0.0
              for (i, p) in enumerate(points)]
    cached2 = zeros(Float64, n)
    advect!(cached2, second, map_on; type=:conservative)
    @test cached2 == reference_advect(second, map_on)
    @test cached2 != cached

    # (3) Conservation over the source support (linear map, exact to round-off).
    @test sum(cached) ≈ sum(source) rtol = 1e-12 atol = 1e-12
    @test sum(cached2) ≈ sum(second) rtol = 1e-12 atol = 1e-12

    # (4) Leakage / over-demand are actually exercised by the fixture.
    @test any(map_on.demand_map .> 1.0)
    @test !isempty(map_on.leakage_indices)

    # (5) Changing the mutable support invalidates the packed cache instead of
    #     silently reusing stale weights: drop one active source node and compare
    #     with the reference that sees the same mutation.
    flip = findfirst(map_on.source_mask)
    map_on.source_mask[flip] = false
    @test map_on.source_mask != map_on.cached_source_mask
    stale_test = zeros(Float64, n)
    advect!(stale_test, source, map_on; type=:conservative)
    @test stale_test == reference_advect(source, map_on)  # fallback recomputes
    @test stale_test != cached                            # mutation genuinely mattered
    map_on.source_mask[flip] = true                       # restore for later tests
    @test map_on.source_mask == map_on.cached_source_mask

    # (6) The raw 9-field positional constructor stays available and uncached.
    map9 = TransportMap(map_on.active_indices, map_on.backward_rays, map_on.demand_map,
        map_on.leakage_indices, map_on.leakage_rays, map_on.grid_meta,
        map_on.source_mask, map_on.target_mask, map_on.is_identity)
    @test !map9.cache_valid
    @test isempty(map9.weight_indices)
    nine = zeros(Float64, n)
    advect!(nine, source, map9; type=:conservative)
    @test nine == reference

    # Cache storage must not alias a field. Resize it to the field length to
    # exercise rejection before any cached read, even after external mutation.
    alias_map, alias_points, _ = compressible_map(8; cache_weights=true)
    resize!(alias_map.weight_values, length(alias_points))
    aliased = alias_map.weight_values
    saved_weights = copy(aliased)
    alias_target = fill(-7.0, length(alias_points))
    @test_throws ArgumentError advect!(aliased, ones(length(alias_points)), alias_map)
    @test aliased == saved_weights
    @test_throws ArgumentError advect!(alias_target, aliased, alias_map)
    @test all(==(-7.0), alias_target)

    # (7) Stationary map stays an identity and never needs the cache.  On a cut
    #     support the identity branch copies the support and zeroes the rest.
    geom_id, _, pts_id = square_geometry(8)
    stationary = TransportMap(geom_id, fill(VectorValue(0.0, 0.0), length(pts_id)), 0.1;
        cache_weights=true)
    @test stationary.is_identity
    src_id = collect(1.0:length(pts_id))
    tgt_id = fill(-1.0, length(pts_id))
    advect!(tgt_id, src_id, stationary; type=:conservative)
    @test tgt_id[stationary.active_indices] == src_id[stationary.active_indices]
    @test all(iszero, tgt_id[.!stationary.source_mask])
end

println("TestTransportWeightCache.jl: all checks passed")
