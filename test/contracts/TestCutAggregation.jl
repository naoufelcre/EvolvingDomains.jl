using Test
using Main.EvolvingDomains

function _aggregation_cut(grid, bg_states; cell_facet_states, areas=Float64[])
    points = NTuple{2,Float64}[]
    point_ids = Int32[]
    point_pointers = Int32[1]
    parents = Int32[]
    subcell_states = Int8[]
    for (cell, area) in enumerate(areas)
        i = (cell - 1) % (length(grid.axes[1]) - 1) + 1
        j = (cell - 1) ÷ (length(grid.axes[1]) - 1) + 1
        x, y = grid.axes[1][i], grid.axes[2][j]
        leg = sqrt(2area)
        append!(points, [(x, y), (x + leg, y), (x, y + leg)])
        append!(point_ids, Int32.(length(points) - 2:length(points)))
        push!(point_pointers, Int32(length(point_ids) + 1))
        push!(parents, Int32(cell))
        push!(subcell_states, Int8(-1))
    end
    reference = NTuple{2,Float64}[]
    for (cell, area) in enumerate(areas)
        i = (cell - 1) % (length(grid.axes[1]) - 1) + 1
        j = (cell - 1) ÷ (length(grid.axes[1]) - 1) + 1
        hx = grid.axes[1][i + 1] - grid.axes[1][i]
        hy = grid.axes[2][j + 1] - grid.axes[2][j]
        leg = sqrt(2area)
        append!(reference, [(0.0, 0.0), (leg / hx, 0.0), (0.0, leg / hy)])
    end
    EvolvingDomains.Cut(grid, :test; bg_states, cell_facet_states,
        subcell_states, cell_to_bgcell=parents, point_ids, point_pointers,
        physical_points=points, reference_points=reference)
end

@testset "Cut-based AgFEM-compatible aggregation" begin
    grid = EvolvingDomains.Grid(([0.0, 1.0, 2.0, 3.0], [0.0, 1.0]))
    # Coherent split cells without facet states: aggregation must reject the missing
    # facet data, not the geometry.
    missing = _aggregation_cut(grid, Int8[0, 0, 0]; cell_facet_states=nothing,
        areas=[0.5, 0.5, 0.5])
    @test_throws ArgumentError EvolvingDomains.aggregateCutCells(missing, grid)
    @test_throws ArgumentError EvolvingDomains.Cut(grid, :test;
        bg_states=fill(Int8(0), 3), cell_facet_states=fill((0, 0, 0), 3))
    @test_throws ArgumentError EvolvingDomains.Cut(grid, :test;
        bg_states=fill(Int16(256), 3))
    @test_throws ArgumentError EvolvingDomains.Cut(grid, :test;
        bg_states=fill(Int8(0), 3), cell_facet_states=fill((Int16(256), 0, 0, 0), 3))
    @test_throws ArgumentError EvolvingDomains.Cut(grid, :test;
        bg_states=fill(Int8(0), 3),
        cell_facet_states=[(1, 0, 1, 1), (1, 1, 1, 1), (1, 1, 1, 1)])
    @test_throws ArgumentError EvolvingDomains.Cut(grid, :test;
        bg_states=fill(Int8(0), 3),
        cell_facet_states=fill((Int8(0), Int8(0), Int8(0), Int8(2)), 3))
    other_grid = EvolvingDomains.Grid(grid.axes)
    with_facets = _aggregation_cut(grid, Int8[0, 0, 0];
        cell_facet_states=fill((Int8(1), Int8(1), Int8(1), Int8(1)), 3), areas=[0.5, 0.5, 0.5])
    @test_throws ArgumentError EvolvingDomains.aggregateCutCells(with_facets, other_grid)

    # Propagation starts at the IN cell and walks only over CUT/IN faces.
    propagation = _aggregation_cut(grid, Int8[-1, 0, 0];
        cell_facet_states=[(1, -1, 1, 1), (-1, 0, 1, 1), (0, 1, 1, 1)])
    @test EvolvingDomains.aggregateCutCells(propagation, grid) == Int32[1, 1, 1]
end

@testset "AgFEM ranking ties and all-cut physical-area fallback" begin
    grid = EvolvingDomains.Grid(([0.0, 1.0, 2.0, 3.0], [0.0, 1.0, 2.0, 3.0]))
    states = fill(Int8(1), 9)
    states[4] = states[6] = -1
    states[5] = 0
    facets = fill((Int8(1), Int8(1), Int8(1), Int8(1)), 9)
    facets[4] = (Int8(1), Int8(-1), Int8(1), Int8(1))
    facets[5] = (Int8(-1), Int8(-1), Int8(1), Int8(1))
    facets[6] = (Int8(-1), Int8(1), Int8(1), Int8(1))
    tied = _aggregation_cut(grid, states; cell_facet_states=facets)
    roots = EvolvingDomains.aggregateCutCells(tied, grid)
    @test roots[5] == 4 # left precedes right when distances tie

    vertical_tie_states = fill(Int8(1), 9)
    vertical_tie_states[4] = vertical_tie_states[8] = -1
    vertical_tie_states[5] = 0
    vertical_tie_facets = fill((Int8(1), Int8(1), Int8(1), Int8(1)), 9)
    vertical_tie_facets[4] = (Int8(1), Int8(-1), Int8(1), Int8(1))
    vertical_tie_facets[5] = (Int8(-1), Int8(1), Int8(1), Int8(-1))
    vertical_tie_facets[8] = (Int8(1), Int8(1), Int8(-1), Int8(1))
    vertical_tie = _aggregation_cut(grid, vertical_tie_states;
        cell_facet_states=vertical_tie_facets)
    @test EvolvingDomains.aggregateCutCells(vertical_tie, grid)[5] == 8

    line = EvolvingDomains.Grid(([0.0, 1.0, 2.0], [0.0, 1.0]))
    all_cut = _aggregation_cut(line, Int8[0, 0];
        cell_facet_states=[(1, 0, 1, 1), (0, 1, 1, 1)], areas=[0.5, 0.125])
    roots = @test_logs (:warn, r"no interior cell") begin
        EvolvingDomains.aggregateCutCells(all_cut, line)
    end
    @test roots == Int32[1, 1]
end

@testset "single-signed cuts: 0 is unassigned, not a root" begin
    # Reference convention (GridapEmbedded's cell_to_cellin::Vector{Int32}): 0 means "no
    # aggregate". A cut cell may never be 0 (the operator errors instead); cells outside
    # the phase are 0; with no interface at all nothing is assigned, and that is warned.
    grid = EvolvingDomains.Grid(([0.0, 1.0, 2.0, 3.0], [0.0, 1.0]))

    outside = EvolvingDomains.Cut(grid, :test; bg_states=fill(CUT_OUT, 3),
        cell_facet_states=fill((CUT_OUT, CUT_OUT, CUT_OUT, CUT_OUT), 3))
    roots = @test_logs (:warn, r"no cell belongs to an aggregate") begin
        EvolvingDomains.aggregateCutCells(outside, grid)
    end
    @test roots == zeros(Int32, 3)
    @test isempty(unique(r for r in roots if r > 0))   # the aggregate set, not unique(roots)

    inside = EvolvingDomains.Cut(grid, :test; bg_states=fill(CUT_IN, 3),
        cell_facet_states=fill((CUT_IN, CUT_IN, CUT_IN, CUT_IN), 3))
    @test @test_logs EvolvingDomains.aggregateCutCells(inside, grid) == Int32[1, 2, 3]
end
