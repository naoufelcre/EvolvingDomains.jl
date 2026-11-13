using Test
using Main.EvolvingDomains: Grid, Cut, Discretization, Field, cut, subcells, cut_cell_ids

@testset "3D Cut snapshot ownership and validation" begin
    grid = Grid(([0.0, 2.0], [0.0, 3.0], [0.0, 4.0]))
    points = [(0.0,0.0,0.0), (1.0,0.0,0.0), (0.0,1.0,0.0), (0.0,0.0,1.0)]
    facets = [(-1,1,0,0,-1,1)]
    cut = Cut(grid, :initial; bg_states=Int8[0], subcell_states=Int8[-1],
        cell_to_bgcell=Int32[1], point_ids=Int32[1,2,3,4], point_pointers=Int32[1,5],
        physical_points=points, reference_points=points, cell_facet_states=facets,
        boundary_facet_to_bgcell=Int32[1], boundary_point_ids=Int32[1,2,3],
        boundary_point_pointers=Int32[1,4],
        boundary_physical_points=points[1:3], boundary_reference_points=points[1:3],
        boundary_normals=[(1.0,0.0,0.0)])
    points[1] = (9.0,9.0,9.0)
    facets[1] = (1,1,1,1,1,1)
    @test cut isa Cut{Float64,3}
    @test cut.physical_points[1] == (0.0,0.0,0.0)
    @test cut.cell_facet_states[1] == (-1,1,0,0,-1,1)
    @test cut.boundary.normals == [(1.0,0.0,0.0)]
    @test cut.boundary.point_pointers == Int32[1,4]
    @test cut_cell_ids(cut) == Int32[1]
    @test collect(subcells(cut, 1)) == Int32[1]
    @test cut.cell_subcell_ptrs == Int32[1, 2]
    _checkSubcellIndex(cut)

    @test_throws ArgumentError Cut(grid, :bad; bg_states=Int8[0],
        subcell_states=Int8[-1], cell_to_bgcell=Int32[1], point_ids=Int32[1,2,3],
        point_pointers=Int32[1,4], physical_points=points[1:3],
        reference_points=points[1:3])
    @test_throws ArgumentError Cut(grid, :bad; bg_states=Int8[0],
        cell_facet_states=[(0,0,0,0)])
    @test_throws ArgumentError Cut(Grid(([0.0,1.0],)), :bad; bg_states=Int8[0])

    multi = Grid(([0.0,1.0,2.0], [0.0,1.0], [0.0,1.0]))
    common = (bg_states=Int8[0,0], cell_facet_states=[(-1,-1,0,0,0,0),(-1,1,0,0,0,0)])
    @test Cut(multi, :ok; common...).cell_facet_states[2][1] == -1
    @test_throws ArgumentError Cut(multi, :bad; bg_states=Int8[0,0],
        cell_facet_states=[(-1,-1,0,0,0,0),(1,1,0,0,0,0)])
end

@testset "3D cut refuses a non-uniform axis before any native call" begin
    axis = collect(range(0.0, 1.0; length=9))
    z = copy(axis); z[5] += 0.01   # strictly increasing, non-uniform
    field = Field(Discretization(Dict(:init => Grid((axis, copy(axis), z)))),
                  Dict(:init => fill(-1.0, 9, 9, 9)))   # single-signed: no straddling cell
    @test_throws ArgumentError cut(field, :init)
end
