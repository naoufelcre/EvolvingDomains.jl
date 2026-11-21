using Test
using Main.EvolvingDomains

@testset "Cut aggregation in 3D" begin
    grid=EvolvingDomains.Grid(([0.0,1,2],[0.0,1],[0.0,1,2]))
    # x-fast cells: [x1,z1], [x2,z1], [x1,z2], [x2,z2]
    states=Int8[-1,0,1,1]
    facets=fill((Int8(1),Int8(1),Int8(1),Int8(1),Int8(1),Int8(1)),4)
    facets[1]=(1,0,1,1,1,1); facets[2]=(0,1,1,1,1,0); facets[4]=(1,1,1,1,0,1)
    cut=EvolvingDomains.Cut(grid,:test;bg_states=states,cell_facet_states=facets)
    @test EvolvingDomains.aggregateCutCells(cut,grid)==Int32[1,1,0,0]
    @test_throws ArgumentError EvolvingDomains.aggregateCutCells(cut,EvolvingDomains.Grid(grid.axes))

    # A CUT cell connects to IN only by its +z face; z faces are sides 5/6.
    zgrid=EvolvingDomains.Grid(([0.0,1],[0.0,1],[0.0,1,2]))
    zfacets=[(1,1,1,1,1,0),(1,1,1,1,0,1)]
    zcut=EvolvingDomains.Cut(zgrid,:test;bg_states=Int8[0,-1],cell_facet_states=zfacets)
    @test EvolvingDomains.aggregateCutCells(zcut,zgrid)==Int32[2,2]

    # Equal-distance eligible roots preserve the existing first-facet tie break (-x first).
    tiegrid=EvolvingDomains.Grid(([0.0,1,2,3],[0.0,1],[0.0,1]))
    tf=fill((Int8(1),Int8(1),Int8(1),Int8(1),Int8(1),Int8(1)),3)
    tf[1]=(1,0,1,1,1,1); tf[2]=(0,0,1,1,1,1); tf[3]=(0,1,1,1,1,1)
    tie=EvolvingDomains.Cut(tiegrid,:test;bg_states=Int8[-1,0,-1],cell_facet_states=tf)
    @test EvolvingDomains.aggregateCutCells(tie,tiegrid)==Int32[1,1,3]

    # All-cut connected component picks greatest physical tetrahedral volume.
    f=[(1,1,1,1,1,0),(1,1,1,1,0,1)]
    fallbackgrid=EvolvingDomains.Grid(([0.0,1],[0.0,1],[0.0,1,3]))
    points=NTuple{3,Float64}[(0,0,0),(1,0,0),(0,1,0),(0,0,1), (0,0,1),(1,0,1),(0,1,1),(0,0,3)]
    fallback=EvolvingDomains.Cut(fallbackgrid,:test;bg_states=Int8[0,0],cell_facet_states=f,
        subcell_states=Int8[-1,-1],cell_to_bgcell=Int32[1,2],
        point_ids=Int32[1,2,3,4,5,6,7,8],point_pointers=Int32[1,5,9],
        physical_points=points,reference_points=points)
    roots=@test_logs (:warn,r"no interior cell") EvolvingDomains.aggregateCutCells(fallback,fallbackgrid)
    @test roots==Int32[2,2]
end
