using Test
using Main.EvolvingDomains: Cut, Grid, Discretization, Field, cut, subcells, cut_cell_ids

function _signed_subcell_areas(cut)
    return [begin
        ids = Int(cut.point_pointers[cell]):(Int(cut.point_pointers[cell + 1]) - 1)
        a, b, c = Int.(cut.point_ids[ids])
        p, q, r = cut.physical_points[[a, b, c]]
        ((q[1] - p[1]) * (r[2] - p[2]) - (q[2] - p[2]) * (r[1] - p[1])) / 2
    end for cell in eachindex(cut.subcell_states)]
end

function _boundary_lengths(cut)
    return [begin
        first_id = Int(cut.boundary.point_ids[cut.boundary.point_pointers[facet]])
        last_id = Int(cut.boundary.point_ids[cut.boundary.point_pointers[facet] + 1])
        p, q = cut.boundary.physical_points[[first_id, last_id]]
        hypot(q[1] - p[1], q[2] - p[2])
    end for facet in eachindex(cut.boundary.facet_to_bgcell)]
end

function _cut_area(cut, state)
    area = 0.0
    for cell in eachindex(cut.subcell_states)
        cut.subcell_states[cell] == state || continue
        ids = Int(cut.point_pointers[cell]):(Int(cut.point_pointers[cell + 1]) - 1)
        a, b, c = Int.(cut.point_ids[ids])
        p, q, r = cut.physical_points[[a, b, c]]
        area += abs((q[1] - p[1]) * (r[2] - p[2]) - (q[2] - p[2]) * (r[1] - p[1])) / 2
    end
    return area
end

@testset "Native CutCells backend" begin
    xs, ys = collect(0.0:1.0:4.0), collect(0.0:1.0:4.0)
    grid = Grid((xs, ys))
    phi = [x - 1.5 for x in xs, y in ys]
    field = Field(Discretization(Dict(:init => grid)), Dict(:init => phi))
    snapshot = cut(field, :init)

    @test unique(snapshot.bg_states) == Int8[-1, 0, 1]
    @test isapprox(_cut_area(snapshot, Int8(-1)), 6.0; atol=1e-12)
    @test isapprox(_cut_area(snapshot, Int8(1)), 10.0; atol=1e-12)
    @test length(snapshot.boundary.facet_to_bgcell) == 4

    # Native producer: the precomputed index agrees with a scan, every cut cell owns at
    # least one volume piece and every single-signed cell owns exactly its full cell.
    ids = cut_cell_ids(snapshot)
    @test ids == sort(ids)
    @test all(c -> !isempty(subcells(snapshot, c)), ids)
    @test sum(c -> length(subcells(snapshot, c)), 1:length(snapshot.bg_states)) ==
          length(snapshot.subcell_states)
    for cell in 1:length(snapshot.bg_states)
        @test collect(subcells(snapshot, cell)) ==
              Int32.(findall(==(Int32(cell)), snapshot.cell_to_bgcell))
        # a single-signed cell owns only pieces of its own phase (its full-cell pieces)
        snapshot.bg_states[cell] == 0 ||
            @test all(==(snapshot.bg_states[cell]), snapshot.subcell_states[subcells(snapshot, cell)])
    end
    @test all(isapprox(n[1], 1.0; atol=1e-12) && isapprox(n[2], 0.0; atol=1e-12)
              for n in snapshot.boundary.normals)

    fill!(phi, NaN)
    @test all(isfinite, first.(snapshot.physical_points))
    @test_throws ArgumentError cut(field, :init)

    aligned_xs = [0.0, 0.5, 1.0]
    aligned_grid = Grid((aligned_xs, aligned_xs))
    aligned_phi = [x - 0.5 for x in aligned_xs, y in aligned_xs]
    aligned_field = Field(Discretization(Dict(:aligned => aligned_grid)), Dict(:aligned => aligned_phi))
    aligned = cut(aligned_field, :aligned)
    @test length(aligned.boundary.facet_to_bgcell) == 2
    @test isapprox(sum(_boundary_lengths(aligned)), 1.0; atol=1e-12)
    @test all(aligned.bg_states[Int(parent)] == 0 for parent in aligned.boundary.facet_to_bgcell)
    @test count(==(Int8(1)), aligned.bg_states) == 2
    @test all(isapprox(p[1], 0.5; atol=1e-12) for p in aligned.boundary.physical_points)
    @test all(isapprox(n[1], 1.0; atol=1e-12) for n in aligned.boundary.normals)

    positive_zero_edge = Field(Discretization(Dict(:positive => Grid(([0.0, 1.0], [0.0, 1.0])))),
        Dict(:positive => [0.0 0.0; 1.0 1.0]))
    positive_zero = cut(positive_zero_edge, :positive)
    @test positive_zero.bg_states == Int8[1]
    @test all(==(Int8(1)), positive_zero.subcell_states)
    @test isempty(positive_zero.boundary.facet_to_bgcell)

    positive_zero_point = Field(Discretization(Dict(:positive_point => Grid(([0.0, 1.0], [0.0, 1.0])))),
        Dict(:positive_point => [0.0 1.0; 1.0 1.0]))
    positive_point = cut(positive_zero_point, :positive_point)
    @test positive_point.bg_states == Int8[1]
    @test all(==(Int8(1)), positive_point.subcell_states)
    @test isempty(positive_point.boundary.facet_to_bgcell)

    negative_zero_point = Field(Discretization(Dict(:negative_point => Grid(([0.0, 1.0], [0.0, 1.0])))),
        Dict(:negative_point => [-1.0 -1.0; -1.0 0.0]))
    negative_point = cut(negative_zero_point, :negative_point)
    @test negative_point.bg_states == Int8[0]
    @test isapprox(_cut_area(negative_point, Int8(-1)), 1.0; atol=1e-12)
    @test isempty(negative_point.boundary.facet_to_bgcell)

    for a in (-1.0, 0.0, 1.0), b in (-1.0, 0.0, 1.0),
        c in (-1.0, 0.0, 1.0), d in (-1.0, 0.0, 1.0)
        samples = (a, b, c, d)
        phi = [a d; b c]
        edge_field = Field(Discretization(Dict(:edge => Grid(([0.0, 1.0], [0.0, 1.0])))),
                           Dict(:edge => phi))
        if all(iszero, samples)
            @test_throws ArgumentError cut(edge_field, :edge)
        else
            snapshot = cut(edge_field, :edge)
            @test isapprox(_cut_area(snapshot, Int8(-1)) + _cut_area(snapshot, Int8(1)),
                           1.0; atol=1e-12)
            @test all(>(0.0), _signed_subcell_areas(snapshot))
            @test all(p -> all(isfinite, p), snapshot.physical_points)
            @test all(n -> all(isfinite, n), snapshot.boundary.normals)
        end
    end

    tangent_values = (1.0, -1.0, 0.0, -1.0)
    tangent_phi = [tangent_values[1] tangent_values[4];
                   tangent_values[2] tangent_values[3]]
    tangent_field = Field(Discretization(Dict(:tangent => Grid(([0.0, 1.0], [0.0, 1.0])))),
                          Dict(:tangent => tangent_phi))
    tangent = cut(tangent_field, :tangent)
    @test isapprox(_cut_area(tangent, Int8(-1)) + _cut_area(tangent, Int8(1)), 1.0; atol=1e-12)
    @test !isempty(tangent.boundary.facet_to_bgcell)
    @test all(>(0.0), _boundary_lengths(tangent))

    for samples in ((-1.0, 1.0, 1.0, -1.0), (-1.0, 1.0, -1.0, 1.0))
        grid = Grid(([0.0, 1.0], [0.0, 1.0]))
        make_field = values -> Field(Discretization(Dict(:scaled => grid)),
                                     Dict(:scaled => [values[1] values[4]; values[2] values[3]]))
        reference = cut(make_field(samples), :scaled)
        for scale in (1e-100, 1e100)
            scaled = cut(make_field(scale .* samples), :scaled)
            @test isapprox(_cut_area(scaled, Int8(-1)), _cut_area(reference, Int8(-1)); atol=1e-12)
            @test isapprox(_cut_area(scaled, Int8(1)), _cut_area(reference, Int8(1)); atol=1e-12)
            @test length(scaled.boundary.normals) == length(reference.boundary.normals)
            @test all(isapprox.(first.(scaled.boundary.normals), first.(reference.boundary.normals); atol=1e-12))
            @test all(isapprox.(last.(scaled.boundary.normals), last.(reference.boundary.normals); atol=1e-12))
        end
    end

    ambiguous_values = [-1.0, 1.0, -1.0, 1.0]
    ambiguous_phi = zeros(2, 2)
    ambiguous_phi[1, 1], ambiguous_phi[2, 1] = ambiguous_values[1], ambiguous_values[2]
    ambiguous_phi[2, 2], ambiguous_phi[1, 2] = ambiguous_values[3], ambiguous_values[4]
    ambiguous_field = Field(Discretization(Dict(:ambiguous => Grid(([0.0, 1.0], [0.0, 1.0])))),
        Dict(:ambiguous => ambiguous_phi))
    ambiguous = cut(ambiguous_field, :ambiguous)
    @test length(ambiguous.boundary.facet_to_bgcell) == 2
    @test isapprox(_cut_area(ambiguous, Int8(-1)) + _cut_area(ambiguous, Int8(1)), 1.0; atol=1e-12)
    @test all(n -> all(isfinite, n), ambiguous.boundary.normals)

    nonuniform = Grid(([0.0, 0.3, 1.0], [0.0, 0.5, 1.0]))
    bad_field = Field(Discretization(Dict(:bad => nonuniform)), Dict(:bad => zeros(3, 3)))
    @test_throws ArgumentError cut(bad_field, :bad)
end
