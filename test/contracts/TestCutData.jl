using Test
using Main.EvolvingDomains

"""Check the precomputed background-cell -> subcell index against a fresh scan."""
function _checkSubcellIndex(cut)
    @test cut.cell_subcell_ptrs[1] == 1
    @test length(cut.cell_subcell_ptrs) == length(cut.bg_states) + 1
    @test cut.cell_subcell_ptrs[end] == length(cut.subcell_states) + 1
    @test issorted(cut.cell_subcell_ptrs)
    for cell in eachindex(cut.bg_states)
        @test collect(EvolvingDomains.subcells(cut, cell)) ==
              Int32.(findall(==(Int32(cell)), cut.cell_to_bgcell))
    end
    @test EvolvingDomains.cutCellIds(cut) ==
          Int32.(findall(==(EvolvingDomains.CUT_CUT), cut.bg_states))
    return nothing
end

@testset "Cut subcell index" begin
    grid = EvolvingDomains.Grid(([0.0, 1.0, 2.0], [0.0, 1.0, 2.0]))
    phys = [(0.0, 0.0), (1.0, 0.0), (0.0, 1.0)]
    nsub = 4
    cut = EvolvingDomains.Cut(grid, :init;
        bg_states=Int8[0, 0, 1, -1], subcell_states=Int8[-1, 1, -1, -1],
        cell_to_bgcell=Int32[1, 1, 2, 4],
        point_ids=Int32.(repeat(1:3, nsub)), point_pointers=Int32.(1:3:(3nsub + 1)),
        physical_points=phys, reference_points=phys)
    @test EvolvingDomains.cutCellIds(cut) == Int32[1, 2]
    @test cut.cell_subcell_ptrs == Int32[1, 3, 4, 4, 5]
    @test collect(EvolvingDomains.subcells(cut, 1)) == Int32[1, 2]
    @test collect(EvolvingDomains.subcells(cut, 3)) == Int32[]      # CUT cell owning none
    @test collect(EvolvingDomains.subcells(cut, 4)) == Int32[4]
    _checkSubcellIndex(cut)

    # The cut-cell view touches only the CSR index, not `cell_to_bgcell`.
    owned = [subcell for cell in EvolvingDomains.cutCellIds(cut)
             for subcell in EvolvingDomains.subcells(cut, cell)]
    @test owned == Int32[1, 2, 3]
    @test owned == findall(<=(2), Int.(cut.cell_to_bgcell))

    # Grouping is what makes the CSR row a subcell id, so the index builder is the single
    # owner of the invariant: a decreasing parent list raises instead of mis-indexing.
    @test_throws ArgumentError EvolvingDomains.Cut(grid, :init;
        bg_states=Int8[0, 0, 1, -1], subcell_states=Int8[-1, 1, -1, -1],
        cell_to_bgcell=Int32[2, 1, 2, 4],
        point_ids=Int32.(repeat(1:3, nsub)), point_pointers=Int32.(1:3:(3nsub + 1)),
        physical_points=phys, reference_points=phys)

    # A boundary-only Cut: cells are classified CUT but own no volume piece.
    boundary_only = EvolvingDomains.Cut(grid, :init; bg_states=Int8[0, 0, 1, -1])
    @test EvolvingDomains.cutCellIds(boundary_only) == Int32[1, 2]
    @test all(==(1), boundary_only.cell_subcell_ptrs)
    _checkSubcellIndex(boundary_only)
end

@testset "Cut snapshot ownership and structure" begin
    grid = EvolvingDomains.Grid(([0.0, 1.0], [0.0, 1.0]))
    bg = Int8[0]
    physical = [(0.0, 0.0), (1.0, 0.0), (0.0, 1.0)]
    cut = EvolvingDomains.Cut(grid, :init;
        bg_states=bg, subcell_states=Int8[-1], cell_to_bgcell=Int32[1],
        point_ids=Int32[1, 2, 3], point_pointers=Int32[1, 4],
        physical_points=physical, reference_points=copy(physical),
        boundary_facet_to_bgcell=Int32[1], boundary_point_ids=Int32[1, 2],
        boundary_point_pointers=Int32[1, 3],
        boundary_physical_points=[(0.0, 0.0), (1.0, 0.0)],
        boundary_reference_points=[(0.0, 0.0), (1.0, 0.0)],
        boundary_normals=[(0.0, -1.0)])

    bg[1] = 1
    physical[1] = (9.0, 9.0)
    @test cut.mesh_id === grid.mesh_id
    @test cut.level_id == :init
    @test cut.bg_states == Int8[0]
    @test cut.physical_points[1] == (0.0, 0.0)
    @test cut.boundary.facet_to_bgcell == Int32[1]
    @test EvolvingDomains.cutCellIds(cut) == Int32[1]
    @test collect(EvolvingDomains.subcells(cut, 1)) == Int32[1]
    _checkSubcellIndex(cut)

    # A Cut without boundary facets. An outside cell owning its full-cell volume piece is
    # what the native producer emits for a cell entirely outside the phase.
    not_cut = EvolvingDomains.Cut(grid, :init;
        bg_states=Int8[1], subcell_states=Int8[-1], cell_to_bgcell=Int32[1],
        point_ids=Int32[1, 2, 3], point_pointers=Int32[1, 4],
        physical_points=cut.physical_points, reference_points=cut.reference_points)
    original = (copy(not_cut.bg_states), copy(not_cut.subcell_states),
        copy(not_cut.cell_to_bgcell), copy(not_cut.point_ids), copy(not_cut.point_pointers),
        copy(not_cut.physical_points), copy(not_cut.reference_points))
    @test (not_cut.bg_states, not_cut.subcell_states, not_cut.cell_to_bgcell,
        not_cut.point_ids, not_cut.point_pointers, not_cut.physical_points,
        not_cut.reference_points) == original

    @test_throws ArgumentError EvolvingDomains.Cut(grid, :init; bg_states=Int8[])
    @test_throws ArgumentError EvolvingDomains.Cut(grid, :init; bg_states=Int8[2])
    @test_throws ArgumentError EvolvingDomains.Cut(grid, :init; bg_states=Int8[0],
        subcell_states=Int8[-1], cell_to_bgcell=Int32[1], point_ids=Int32[1, 2, 3],
        point_pointers=Int32[1, 2], physical_points=[(0.0, 0.0), (1.0, 0.0), (0.0, 1.0)],
        reference_points=[(0.0, 0.0), (1.0, 0.0), (0.0, 1.0)])
    @test_throws ArgumentError EvolvingDomains.Cut(grid, :init; bg_states=Int8[0],
        physical_points=[(0.0, 0.0, 1.0)], reference_points=[(0.0, 0.0)])
    @test_throws ArgumentError EvolvingDomains.Cut(grid, :init; bg_states=Int8[0],
        boundary_facet_to_bgcell=Int32[1], boundary_point_ids=Int32[1, 2],
        boundary_point_pointers=Int32[1, 3],
        boundary_physical_points=[(0.0, 0.0), (1.0, 0.0)],
        boundary_reference_points=[(0.0, 0.0), (1.0, 0.0)],
        boundary_normals=[(0.0, 1.0, 0.0)])
    @test_throws ArgumentError EvolvingDomains.Cut(
        EvolvingDomains.Grid(([0.0, 1.0],)), :init; bg_states=Int8[0])
end
