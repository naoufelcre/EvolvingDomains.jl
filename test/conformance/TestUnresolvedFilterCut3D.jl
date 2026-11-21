using Test
using Main.EvolvingDomains: Grid, Discretization, Field, cut, filterUnresolvedFeatures!,
    cutCellIds, CUT_IN, CUT_OUT, CUT_CUT

# Conformance of `filterUnresolvedFeatures!` in 3D against the native cut producer (Luo, Hu &
# Adams 2016, §3; eps = 0.75h, 3x3x3 neighbourhood, 9x9x9 search, project island fallback).
# Resolved geometry must survive bit-identically, an isolated sub-eps island must disappear
# from the *rebuilt* cut, a thin attached feature must behave as the published branches do,
# and a stale `Cut` must be rejected before the field is mutated.
#
# Helper names are prefixed to stay clear of the other conformance files, which share `Main`.

_filtered3DAxis(n) = collect(range(0.0, 1.0; length=n))

# Uniform grid on the unit cube, plus a Field holding one level `:init`.
function _filtered3DGrid(n)
    ax = _filtered3DAxis(n)
    return Grid((copy(ax), copy(ax), copy(ax)))
end

function _filtered3DField(n, f)
    grid = _filtered3DGrid(n)
    ax = grid.axes[1]
    return Field(Discretization(Dict(:init => grid)),
                 Dict(:init => [f(x, y, z) for x in ax, y in ax, z in ax]))
end

# Volume of the subcells carrying `state`: tetrahedra in 3D, |det|/6 of the four vertices.
function _filtered3DSubcellVolume(snapshot, state)
    total = 0.0
    for cell in eachindex(snapshot.subcell_states)
        snapshot.subcell_states[cell] == state || continue
        ids = Int(snapshot.point_pointers[cell]):(Int(snapshot.point_pointers[cell + 1]) - 1)
        p, q, r, s = snapshot.physical_points[Int.(snapshot.point_ids[ids])]
        a, b, c = ntuple(d -> q[d] - p[d], 3), ntuple(d -> r[d] - p[d], 3),
                 ntuple(d -> s[d] - p[d], 3)
        total += abs(a[1] * (b[2] * c[3] - b[3] * c[2]) - a[2] * (b[1] * c[3] - b[3] * c[1]) +
                     a[3] * (b[1] * c[2] - b[2] * c[1])) / 6
    end
    return total
end

# Sum of the interface facet areas (the discrete measure of the zero level set).
function _filtered3DInterfaceArea(snapshot)
    total = 0.0
    surface = snapshot.boundary
    for facet in eachindex(surface.facet_to_bgcell)
        ids = Int(surface.point_pointers[facet]):(Int(surface.point_pointers[facet + 1]) - 1)
        p, q, r = surface.physical_points[Int.(surface.point_ids[ids])]
        a, b = ntuple(d -> q[d] - p[d], 3), ntuple(d -> r[d] - p[d], 3)
        total += hypot(a[2] * b[3] - a[3] * b[2], a[3] * b[1] - a[1] * b[3],
                       a[1] * b[2] - a[2] * b[1]) / 2
    end
    return total
end

_filtered3DFacetCentres(snapshot) = [begin
    ids = Int(snapshot.boundary.point_pointers[f]):(Int(snapshot.boundary.point_pointers[f + 1]) - 1)
    p, q, r = snapshot.boundary.physical_points[Int.(snapshot.boundary.point_ids[ids])]
    ntuple(d -> (p[d] + q[d] + r[d]) / 3, 3)
end for f in eachindex(snapshot.boundary.facet_to_bgcell)]

# Background states recomputed from the nodal signs of the current field, independently of the
# producer. Linear index: first axis fastest, matching `cutCellIds`.
function _filtered3DStates(phi)
    nc = size(phi) .- 1
    states = Vector{Int8}(undef, prod(nc))
    k = 0
    for c in 1:nc[3], b in 1:nc[2], a in 1:nc[1]
        k += 1
        v = (phi[a, b, c], phi[a + 1, b, c], phi[a, b + 1, c], phi[a + 1, b + 1, c],
             phi[a, b, c + 1], phi[a + 1, b, c + 1], phi[a, b + 1, c + 1], phi[a + 1, b + 1, c + 1])
        states[k] = all(>(0), v) ? CUT_OUT : (all(<(0), v) ? CUT_IN : CUT_CUT)
    end
    return states
end

# Every background cell of the snapshot is classified exactly as the field's corner signs say,
# and every cut cell carries a boundary facet. A stale cut would fail the first check.
function _filtered3DClassificationMatches(snapshot, phi)
    snapshot.bg_states == _filtered3DStates(phi) || return false
    owners = Set(Int.(snapshot.boundary.facet_to_bgcell))
    return all(cell -> snapshot.bg_states[cell] != CUT_CUT || cell in owners,
               eachindex(snapshot.bg_states))
end

# Physical centre of a background cell from its linear id (first axis fastest).
function _filtered3DCellCentre(ax, id, cell_dims)
    i = (id - 1) % cell_dims[1] + 1
    j = (id - 1) ÷ cell_dims[1] % cell_dims[2] + 1
    k = (id - 1) ÷ (cell_dims[1] * cell_dims[2]) + 1
    h = ax[2] - ax[1]
    return (ax[i] + h / 2, ax[j] + h / 2, ax[k] + h / 2)
end

_distance3(a, b) = hypot(a[1] - b[1], a[2] - b[2], a[3] - b[3])

@testset "3D unresolved-feature filter with the native cut" begin

    # Sub-resolution positive island at a node inside a resolved ball of the negative phase,
    # and its exact negative mirror. eps = 0.75/16 = 0.046875 and a = 0.03 < eps, so the
    # island's phase never reaches the auxiliary level set.
    island_h = 1.0 / 16
    island_eps = 0.75island_h
    island_a = 0.03
    island_c = (0.4375, 0.5, 0.5)          # a node of the 17^3 grid, r = 0.0625 < 0.3
    island_ball(x, y, z) = hypot(x - 0.5, y - 0.5, z - 0.5) - 0.3
    island_drop(x, y, z) = hypot(x - island_c[1], y - island_c[2], z - island_c[3])
    island_pos(x, y, z) = max(island_ball(x, y, z), island_a - island_drop(x, y, z))
    island_neg(x, y, z) = min(-island_ball(x, y, z), island_drop(x, y, z) - island_a)

    @testset "a resolved sphere is not written" begin
        # Every cut cell of the sphere has both auxiliary crossings (phi = +-eps) inside its
        # 3x3x3 neighbourhood, so it is S2 and neither Eq. (15) nor the fallback touches it.
        # No node is written: the field, and therefore the cut built from it, must agree bit
        # for bit (exact equality, no tolerance).
        for n in (17, 25)
            field = _filtered3DField(n, (x, y, z) -> hypot(x - 0.5, y - 0.5, z - 0.5) - 0.3)
            before = copy(field.data[:init])
            initial = cut(field, :init)

            @test filterUnresolvedFeatures!(field, :init) === field
            filtered = cut(field, :init)
            area = (_filtered3DInterfaceArea(initial), _filtered3DInterfaceArea(filtered))
            @info "resolved sphere" n cutcells=length(cutCellIds(initial)) facets=length(initial.boundary.facet_to_bgcell) area=area changed=count(!=(0.0), field.data[:init] .- before)

            @test field.data[:init] == before
            @test maximum(abs.(field.data[:init] .- before)) == 0.0     # the reported number
            @test filtered.bg_states == initial.bg_states
            @test length(filtered.boundary.facet_to_bgcell) == length(initial.boundary.facet_to_bgcell)
            @test area[1] == area[2]
            @test filtered.boundary.physical_points == initial.boundary.physical_points

            # The sphere really is cut geometry, and the measure is the sphere's surface.
            @test length(cutCellIds(initial)) > 400
            @test isapprox(area[1], 4pi * 0.3^2; rtol=0.02)
        end
    end

    @testset "an isolated island of either sign is removed" begin
        h, eps = island_h, island_eps
        @test island_a < eps
        cell_dims = (16, 16, 16)
        ax = _filtered3DAxis(17)
        near_centre(id) = _distance3(_filtered3DCellCentre(ax, id, cell_dims), island_c) < 0.1

        positive = _filtered3DField(17, island_pos)
        negative = _filtered3DField(17, island_neg)
        before_pos = copy(positive.data[:init])
        before_neg = copy(negative.data[:init])
        initial = cut(positive, :init)

        island_cells = [id for id in cutCellIds(initial) if near_centre(id)]
        island_nodes = [node for node in CartesianIndices(size(before_pos))
                        if _distance3(ntuple(d -> ax[node[d]], 3), island_c) < 0.1]
        @info "island before filtering" eps island_max=maximum(before_pos[island_nodes]) island_positives=count(c -> before_pos[c] > 0, island_nodes) cutcells=length(cutCellIds(initial)) island_cutcells=length(island_cells)
        @test maximum(before_pos[island_nodes]) == island_a           # the node holding the island
        @test count(c -> before_pos[c] > 0, island_nodes) == 1
        @test count(c -> before_neg[c] < 0, island_nodes) == 1        # the mirror holds one node
        @test !isempty(island_cells)                                  # it is cut geometry
        @test !isempty([c for c in _filtered3DFacetCentres(initial) if _distance3(c, island_c) < 0.1])

        filterUnresolvedFeatures!(positive, :init)
        filterUnresolvedFeatures!(negative, :init)
        after_pos = positive.data[:init]
        filtered = cut(positive, :init)

        # The island is gone, from the field and from the rebuilt cut, and only its cells went.
        @test count(c -> after_pos[c] > 0, island_nodes) == 0
        @test maximum(after_pos[island_nodes]) < 0
        @test isempty([id for id in cutCellIds(filtered) if near_centre(id)])
        @test length(cutCellIds(filtered)) == length(cutCellIds(initial)) - length(island_cells)
        @test isempty([c for c in _filtered3DFacetCentres(filtered) if _distance3(c, island_c) < 0.1])
        @test all(isfinite, after_pos) && all(isfinite, negative.data[:init])

        # No stale geometry: the rebuilt classification agrees with the mutated field.
        @test _filtered3DClassificationMatches(filtered, after_pos)
        @test _filtered3DStates(after_pos) != initial.bg_states       # the field did change

        # The published Eq. (15) branch did this, not the island fallback: the corrected nodes
        # are pushed below -eps, whereas the fallback copies a neighbouring value of the
        # surrounding phase, which here is only about -0.01 to -0.08.
        changed = [c for c in CartesianIndices(size(before_pos)) if after_pos[c] != before_pos[c]]
        @info "island removed" changed_nodes=length(changed) after_values=extrema([after_pos[c] for c in changed]) fallback_range=extrema([before_pos[c] for c in changed])
        @test !isempty(changed)
        @test all(c -> after_pos[c] < -eps, changed)

        # The two phases behave symmetrically: negating the input negates the output.
        @test negative.data[:init] == -after_pos

        # Volume/measure: filtering only removes material here. The delta is the island's own
        # subcell volume, so the sphere's volume is preserved to that accuracy.
        v0 = _filtered3DSubcellVolume(initial, CUT_OUT)
        v1 = _filtered3DSubcellVolume(filtered, CUT_OUT)
        @info "island volume" h3=(h^3) in_phase_before_h3=v0 / h^3 in_phase_after_h3=v1 / h^3 delta_h3=(v1 - v0) / h^3 interface_before=_filtered3DInterfaceArea(initial) interface_after=_filtered3DInterfaceArea(filtered)
        @test v1 <= v0
        @test v1 < v0
        @test abs(v1 - v0) < 2h^3                # measured -0.25 h^3: only the island went
        @test _filtered3DInterfaceArea(filtered) <= _filtered3DInterfaceArea(initial)
    end

    @testset "a one-node-thick rod attached to a resolved ball" begin
        # A capsule of the negative phase (half-width 0.5h, so |phi| <= 0.5h < eps across it)
        # protrudes from a resolved ball of radius 0.2; it ends at x = 0.9, i.e. it sticks out
        # by 3 cells, and it is one node thick (its axis is the node line y = z = 0.5).
        h = island_h
        eps = island_eps
        rod(x, y, z) = min(hypot(x - 0.5, y - 0.5, z - 0.5) - 0.2,
                           hypot(max(abs(x - 0.7) - 0.2, 0.0), hypot(y - 0.5, z - 0.5)) - 0.5h)
        field = _filtered3DField(17, rod)
        before = copy(field.data[:init])
        initial = cut(field, :init)
        filterUnresolvedFeatures!(field, :init)
        after = field.data[:init]
        filtered = cut(field, :init)

        changed = [c for c in CartesianIndices(size(before)) if after[c] != before[c]]
        xmin = minimum(c -> _filtered3DAxis(17)[c[1]], changed)
        @info "thin rod" changed_nodes=length(changed) changed_x_min=xmin cutcells=(length(cutCellIds(initial)), length(cutCellIds(filtered))) negatives=(count(<(0), before), count(<(0), after)) axis=[(i, before[i, 9, 9], after[i, 9, 9]) for i in 12:17]
        @test _filtered3DClassificationMatches(filtered, after)

        # Measured: partially corrected, never deleted. The protruding cells are S1 with
        # side = -1 — their 3x3x3 neighbourhood misses the phi = -eps crossing because the rod
        # is one node thick — so the published Eq. (15) re-distances their nodes from the
        # ball's -eps shell found inside the 9x9x9 search region. The result is positive
        # there (d > eps), so the protruding part keeps only its innermost node.
        @test length(changed) == 27
        @test xmin >= 0.8                                   # only the outer part is rewritten
        @test all(c -> abs(_filtered3DAxis(17)[c[2]] - 0.5) <= 0.0625 &&
                       abs(_filtered3DAxis(17)[c[3]] - 0.5) <= 0.0625, changed)
        @test count(c -> before[c] < 0, changed) == 2
        @test all(c -> after[c] > 0, [c for c in changed if before[c] < 0])
        @test all(c -> after[c] > eps, [c for c in changed if before[c] < 0])

        # What survives: the two rod nodes closest to the ball keep their sub-eps value,
        # because their cells still see both auxiliary crossings and are S2 (Eq. 10) — the
        # same resolved-stub effect the 2D conformance test records at a junction.
        @test after[13, 9, 9] == before[13, 9, 9] == -0.5h
        @test after[12, 9, 9] == before[12, 9, 9]
        # The tip node (x = 1.0) lies outside the corrected cells and keeps its value.
        @test after[17, 9, 9] == before[17, 9, 9] > 0
        @test count(<(0), after) == count(<(0), before) - 2
        @test count(<(0), after) > 0                         # a stub is still there

        # Measure: 8 cut cells fewer, and the in-phase volume does not grow.
        @test length(cutCellIds(filtered)) == length(cutCellIds(initial)) - 8
        v0 = _filtered3DSubcellVolume(initial, CUT_IN)
        v1 = _filtered3DSubcellVolume(filtered, CUT_IN)
        @info "thin rod volume" in_phase_before_h3=v0 / h^3 in_phase_after_h3=v1 / h^3 delta_h3=(v1 - v0) / h^3 interface_before=_filtered3DInterfaceArea(initial) interface_after=_filtered3DInterfaceArea(filtered)
        @test v1 <= v0
        @test _filtered3DInterfaceArea(filtered) <= _filtered3DInterfaceArea(initial)
    end

    @testset "a field that never reaches eps is rejected" begin
        # 3D, finite, Float64, with a genuine zero crossing (a cut exists, so the single-signed
        # early return is not what fires), but |phi| = 0.01 < eps = 0.09375 everywhere: no
        # auxiliary level set exists, so the call must refuse instead of classifying the whole
        # phase as unresolved. No cut is needed for this check.
        field = _filtered3DField(9, (x, y, z) -> x < 0.5 ? -0.01 : 0.01)
        eps = 0.75 * (1.0 / 8)
        @test eltype(field.data[:init]) === Float64
        @test all(isfinite, field.data[:init])
        @test maximum(abs, field.data[:init]) < eps
        @test !isempty(cutCellIds(cut(field, :init)))
        snapshot = copy(field.data[:init])
        @test_throws ArgumentError filterUnresolvedFeatures!(field, :init)
        @test field.data[:init] == snapshot
    end
end
