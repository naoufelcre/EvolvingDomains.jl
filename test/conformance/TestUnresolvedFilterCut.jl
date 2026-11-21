using Test
using StaticArrays: SVector
using Main.EvolvingDomains: Grid, Discretization, Field, cut, filterUnresolvedFeatures!,
    reinitialize!, eulerianTransport!, cutCellIds, CUT_IN, CUT_OUT, CUT_CUT, WENO5Cache

# Conformance of `filterUnresolvedFeatures!` against the native cut producer: the cut must be
# rebuilt consistently after filtering, and the filter must behave as published on a dynamic
# single-vortex benchmark (Luo, Hu & Adams 2016, §4.1).
#
# Helper names are prefixed to stay clear of the other conformance files, which share `Main`.

# Area of the subcells carrying `state`: triangles in 2D, so |det|/2 of the three vertices.
function _filteredPhaseArea(snapshot, state)
    total = 0.0
    for cell in eachindex(snapshot.subcell_states)
        snapshot.subcell_states[cell] == state || continue
        ids = Int(snapshot.point_pointers[cell]):(Int(snapshot.point_pointers[cell + 1]) - 1)
        a, b, c = Int.(snapshot.point_ids[ids])
        p, q, r = snapshot.physical_points[[a, b, c]]
        total += abs((q[1] - p[1]) * (r[2] - p[2]) - (q[2] - p[2]) * (r[1] - p[1])) / 2
    end
    return total
end

# Sum of the interface facet lengths (the discrete measure of the zero level set).
function _filteredInterfaceLength(snapshot)
    total = 0.0
    surface = snapshot.boundary
    for facet in eachindex(surface.facet_to_bgcell)
        i = Int(surface.point_ids[surface.point_pointers[facet]])
        j = Int(surface.point_ids[surface.point_pointers[facet] + 1])
        p, q = surface.physical_points[[i, j]]
        total += hypot(q[1] - p[1], q[2] - p[2])
    end
    return total
end

# Background states recomputed from the nodal signs of the current field, independently of the
# producer. Linear index: first axis fastest, matching `cutCellIds`.
function _filteredFieldStates(phi)
    nc1, nc2 = size(phi, 1) - 1, size(phi, 2) - 1
    states = Vector{Int8}(undef, nc1 * nc2)
    k = 0
    for j in 1:nc2, i in 1:nc1
        k += 1
        v = (phi[i, j], phi[i + 1, j], phi[i, j + 1], phi[i + 1, j + 1])
        states[k] = all(>(0), v) ? CUT_OUT : (all(<(0), v) ? CUT_IN : CUT_CUT)
    end
    return states
end

# Independent re-implementation of Eqs. (7)–(10): the corner nodes of every cut cell whose
# 3×3 neighbourhood lacks a +eps crossing or lacks a -eps crossing. eps = 0.75h.
function _filteredUnresolvedNodes(phi, eps)
    nc1, nc2 = size(phi, 1) - 1, size(phi, 2) - 1
    straddles(level) = [begin
        v = (phi[i, j], phi[i + 1, j], phi[i, j + 1], phi[i + 1, j + 1])
        any(x -> x < level, v) && any(x -> x > level, v)
    end for i in 1:nc1, j in 1:nc2]
    S, P, N = straddles(0.0), straddles(eps), straddles(-eps)
    nodes = falses(size(phi))
    for i in 1:nc1, j in 1:nc2
        S[i, j] || continue
        box = (max(1, i - 1):min(nc1, i + 1), max(1, j - 1):min(nc2, j + 1))
        any(P[box...]) && any(N[box...]) && continue
        nodes[i:i+1, j:j+1] .= true
    end
    return nodes
end

# Every background cell of the snapshot is classified exactly as the field's corner signs say,
# and every cut cell carries a boundary facet. A stale cut would fail the first check.
function _filteredClassificationMatches(snapshot, phi)
    snapshot.bg_states == _filteredFieldStates(phi) || return false
    owners = Set(Int.(snapshot.boundary.facet_to_bgcell))
    return all(cell -> snapshot.bg_states[cell] != CUT_CUT || cell in owners,
               eachindex(snapshot.bg_states))
end

@testset "unresolved-feature filter with the native cut" begin

    @testset "filtering leaves no stale cut geometry" begin
        ax = collect(range(0.0, 1.0; length=65))
        grid = Grid((copy(ax), copy(ax)))
        d = Discretization(Dict(:init => grid))
        h = 1.0 / (length(ax) - 1)
        eps = 0.75h
        # Resolved disk (negative phase inside) with one isolated positive drop inside it,
        # max|phi| = 0.0052 < eps = 0.0117, so the drop is a sub-resolution island.
        phi = [max(hypot(x - 0.5, y - 0.5) - 0.3, 0.010 - hypot(x - 0.4, y - 0.4))
               for x in ax, y in ax]
        field = Field(d, Dict(:init => phi))
        before = cut(field, :init)
        @test _filteredClassificationMatches(before, phi)

        islandcells = [c for c in cutCellIds(before)
                       if hypot(ax[(c - 1) % 64 + 1] - 0.4, ax[(c - 1) ÷ 64 + 1] - 0.4) < 0.1]
        dropmax = maximum(phi[i, j] for i in 1:65, j in 1:65 if hypot(ax[i] - 0.4, ax[j] - 0.4) < 0.1)
        @info "cut before filtering" cutcells=length(cutCellIds(before)) islandcells=length(islandcells) maxdrop=dropmax eps=eps
        @test 0 < dropmax < eps                        # sub-resolution island
        @test !isempty(islandcells)                    # the island is cut geometry
        @test length(cutCellIds(before)) > 100         # the disk is cut geometry too

        filterUnresolvedFeatures!(field, :init)
        after = cut(field, :init)
        # No stale geometry: the rebuilt classification agrees with the mutated field.
        @test _filteredClassificationMatches(after, field.data[:init])
        @test _filteredFieldStates(field.data[:init]) != before.bg_states   # the field did change
        # The island disappeared from the cut, and only the island cells did.
        remaining = [c for c in cutCellIds(after)
                     if hypot(ax[(c - 1) % 64 + 1] - 0.4, ax[(c - 1) ÷ 64 + 1] - 0.4) < 0.1]
        @test isempty(remaining)
        @test length(cutCellIds(after)) == length(cutCellIds(before)) - length(islandcells)
        # No facet is left near the island either: the geometry, not just the classification,
        # was rebuilt.
        surface = after.boundary
        neardrop = [f for f in eachindex(surface.facet_to_bgcell)
                    if hypot(surface.physical_points[Int(surface.point_ids[surface.point_pointers[f]])][1] - 0.4,
                             surface.physical_points[Int(surface.point_ids[surface.point_pointers[f]])][2] - 0.4) < 0.1]
        @test isempty(neardrop)
        @test all(isfinite, field.data[:init])
    end

    @testset "single-vortex benchmark: 97², t = 2" begin
        # §4.1: circle R = 0.15 at (0.5, 0.75) in the unit square, velocity Eq. (16). The paper
        # uses 512²; 97² is enough to produce a filament and keeps the run at a few seconds
        # (measured 2–3 s here, advection + 63 reinitialisations; the filter calls add ~3 %).
        n = 97
        ax = collect(range(0.0, 1.0; length=n))
        grid = Grid((copy(ax), copy(ax)))
        d = Discretization(Dict(:init => grid))
        h = 1.0 / (n - 1)
        disk(x, y) = hypot(x - 0.5, y - 0.75) - 0.15
        u(x, y) = 2sinpi(x)^2 * sinpi(y) * cospi(y)
        v(x, y) = -2sinpi(y)^2 * sinpi(x) * cospi(x)
        velocity = vec([SVector(u(x, y), v(x, y)) for x in ax, y in ax])
        cfl = 0.4
        dt = cfl * h / maximum(abs(w[1]) + abs(w[2]) for w in velocity)
        steps = ceil(Int, 2.0 / dt)
        every = 10                                   # filter+reinitialise cadence (as published:
                                                     # filter immediately before reinitialisation)

        function advect(filt::Bool)
            field = Field(d, Dict(:init => [disk(x, y) for x in ax, y in ax]))
            cache = WENO5Cache()
            for step in 1:steps
                eulerianTransport!(field, velocity, dt; level_id=:init, cache)
                if step % every == 0
                    filt && filterUnresolvedFeatures!(field, :init)
                    reinitialize!(field, cut(field, :init))
                end
            end
            return field
        end

        initial = Field(d, Dict(:init => [disk(x, y) for x in ax, y in ax]))
        snapshot0 = cut(initial, :init)
        area0 = _filteredPhaseArea(snapshot0, CUT_IN)
        @test 0.9 * pi * 0.15^2 < area0 < 1.1 * pi * 0.15^2      # subcell area of the disk

        elapsed = @elapsed begin
            plain = advect(false)
            kept = advect(true)
        end

        cplain, ckept = cut(plain, :init), cut(kept, :init)
        ratio_plain = _filteredPhaseArea(cplain, CUT_IN) / area0
        ratio_kept = _filteredPhaseArea(ckept, CUT_IN) / area0
        length_plain, length_kept = _filteredInterfaceLength(cplain), _filteredInterfaceLength(ckept)
        @info "single vortex (Eq. 16, 97², t=2)" steps=steps dt=dt seconds=elapsed area_ratio_no_filter=ratio_plain area_ratio_filtered=ratio_kept interface_no_filter=length_plain interface_filtered=length_kept cutcells_no_filter=length(cutCellIds(cplain)) cutcells_filtered=length(cutCellIds(ckept))

        # Eq. (17) is not conserved by a level-set method, filtered or not: the paper says so
        # explicitly (§4.1). Values are reported above; the checks below only exclude nonsense.
        @test 0 < ratio_plain <= 1.05
        @test 0 < ratio_kept <= 1.05
        @test all(isfinite, plain.data[:init]) && all(isfinite, kept.data[:init])
        # The vertex flow stretches the circle into a filament: the interface grows.
        @test length_plain > 2 * _filteredInterfaceLength(snapshot0)

        # Filtering removes sub-resolution filament cells: fewer cut cells, shorter interface,
        # less area than the same advection without the filter (measured 0.54 vs 0.67).
        @test length(cutCellIds(ckept)) < length(cutCellIds(cplain))
        @test length_kept < length_plain
        @test ratio_kept < ratio_plain

        # What one filter application on the unfiltered field removes, checked against an
        # independent re-implementation of Eqs. (7)-(10): only nodes of unresolved cells move.
        eps = 0.75h
        unresolved = _filteredUnresolvedNodes(plain.data[:init], eps)
        reserved = copy(plain.data[:init])
        cells_before = length(cutCellIds(cut(plain, :init)))
        filterUnresolvedFeatures!(plain, :init)
        cells_after = length(cutCellIds(cut(plain, :init)))
        changed = [CartesianIndex(i, j) for i in 1:n, j in 1:n if plain.data[:init][i, j] != reserved[i, j]]
        @info "one filter at t=2" changed_nodes=length(changed) cutcells=(cells_before, cells_after) unresolved_nodes=count(unresolved)
        @test cells_after < cells_before                  # filament cells are removed
        @test !isempty(changed)
        @test all(unresolved[I] for I in changed)          # nothing resolved is rewritten
        @test all(isfinite, plain.data[:init])

        # Repeated filtering does not grow the interface: two further filter+reinitialise
        # passes (the published loop) shorten it or leave it alone, and never lengthen it.
        lengths = [_filteredInterfaceLength(cut(kept, :init))]
        areas = [_filteredPhaseArea(cut(kept, :init), CUT_IN)]
        for _ in 1:2
            filterUnresolvedFeatures!(kept, :init)
            reinitialize!(kept, cut(kept, :init))
            push!(lengths, _filteredInterfaceLength(cut(kept, :init)))
            push!(areas, _filteredPhaseArea(cut(kept, :init), CUT_IN))
        end
        @info "repeated filtering" interface=lengths area=areas
        @test lengths[2] <= lengths[1]
        @test lengths[3] <= lengths[2]
        @test areas[2] <= areas[1] * (1 + 1e-12)
        @test areas[3] <= areas[2] * (1 + 1e-12)
        # It converges: the third pass is the identity on a reinitialised field.
        @test lengths[3] == lengths[2] && areas[3] == areas[2]
        @test all(isfinite, kept.data[:init])
    end
end
