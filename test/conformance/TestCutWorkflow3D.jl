using Test
using Main.EvolvingDomains: Grid, Discretization, Field, cut, eulerianTransport!, lagrangianTransport!,
    TransportMap, aggregateCutCells, reinitialize!
using StaticArrays: SVector

@testset "3D native Cut / Field combined workflow" begin
    axis = collect(range(-1.0,1.0;length=13))
    g = Grid((axis,copy(axis),copy(axis)))
    d = Discretization(Dict(:init=>g))
    levelset = Field(d,Dict(:init=>[sqrt((x-.03)^2+(y+.02)^2+(z-.01)^2)-.47 for x in axis,y in axis,z in axis]))
    previous = cut(levelset,:init)
    saved_points = copy(previous.physical_points)
    velocity = [SVector(.02sin(π*(x+1)/2),-.015sin(π*(y+1)/2),.01sin(π*(z+1)/2))
        for z in axis for y in axis for x in axis]
    dt = .005
    eulerianTransport!(levelset,velocity,dt;level_id=:init)
    current = cut(levelset,:init)
    @test previous.mesh_id === current.mesh_id === g.mesh_id
    @test previous.physical_points == saved_points
    @test !isempty(current.boundary.facet_to_bgcell)

    source = Field(d,Dict(:init=>[exp(-x^2-y^2-z^2) for x in axis,y in axis,z in axis]))
    target = Field(d,Dict(:init=>zeros(13,13,13)))
    map = TransportMap(previous,current,g,velocity,dt;cache_weights=true)
    @test lagrangianTransport!(target,source,map) === target
    @test sum(target.data[:init]) ≈ sum(vec(source.data[:init])[map.source_mask]) rtol=2e-13
    @test all(iszero,vec(target.data[:init])[.!map.target_mask])

    intensive = Field(d,Dict(:init=>zeros(13,13,13)))
    @test lagrangianTransport!(intensive,source,velocity,dt) === intensive
    @test all(isfinite,intensive.data[:init])
    @test intensive.data[:init] != source.data[:init]

    roots = aggregateCutCells(current,g)
    @test all(roots[current.bg_states .!= 1] .> 0)
    @test all(iszero,roots[current.bg_states .== 1])
    signs = sign.(levelset.data[:init])
    surface_points = copy(current.boundary.physical_points)
    @test reinitialize!(levelset,current) === levelset
    @test all(isfinite,levelset.data[:init])
    @test sign.(levelset.data[:init]) == signs
    @test current.boundary.physical_points == surface_points
end
