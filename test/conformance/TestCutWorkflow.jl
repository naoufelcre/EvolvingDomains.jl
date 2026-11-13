using Test
using Main.EvolvingDomains: Grid, Discretization, Field, Cut, cut, eulerianTransport!,
    TransportMap, lagrangianTransport!, aggregateCutCells, reinitialize!
using StaticArrays: SVector

@testset "Native Cut / Field single-step workflow" begin
    axis = collect(range(-1.0, 1.0; length=25))
    grid = Grid((axis, copy(axis)))
    d = Discretization(Dict(:init => grid))
    levelset = Field(d, Dict(:init => [hypot(x - 0.03, y + 0.02) - 0.47 for x in axis, y in axis]))
    previous = cut(levelset, :init)
    previous_points = copy(previous.physical_points)
    velocity = fill(SVector(0.05, 0.0), length(axis)^2)
    dt = 0.01
    eulerianTransport!(levelset, velocity, dt)
    current = cut(levelset, :init)
    @test current.mesh_id === previous.mesh_id === grid.mesh_id
    @test previous.physical_points == previous_points

    source = Field(d, Dict(:init => [exp(-x^2 - y^2) for x in axis, y in axis]))
    target = Field(d, Dict(:init => zeros(length(axis), length(axis))))
    map = TransportMap(previous, current, grid, velocity, dt; cache_weights=true)
    @test lagrangianTransport!(target, source, map) === target
    @test sum(target.data[:init]) ≈ sum(vec(source.data[:init])[map.source_mask]) rtol=1e-12
    @test all(iszero, vec(target.data[:init])[.!map.target_mask])

    roots = aggregateCutCells(current, grid)
    @test all(roots[current.bg_states .!= 1] .> 0)
    @test all(iszero, roots[current.bg_states .== 1])
    boundary_points = copy(current.boundary.physical_points)
    reinitialize!(levelset, current)
    @test all(isfinite, levelset.data[:init])
    @test current.boundary.physical_points == boundary_points
end

@testset "Native producer and checking constructor agree" begin
    # D5: the binding returns Cut{T,N}(Val(:validated), ...) without re-validating or
    # re-copying. Rebuilding the same snapshot through the checking constructor must yield
    # identical arrays, so the fast path cannot drift from the checked one.
    axis = collect(range(0.0, 1.0; length=13))
    grid = Grid((axis, copy(axis)))
    levelset = Field(Discretization(Dict(:init => grid)),
        Dict(:init => [hypot(x - 0.5, y - 0.5) - 0.3 for x in axis, y in axis]))
    native = cut(levelset, :init)
    checked = Cut(grid, :init;
        bg_states=native.bg_states, subcell_states=native.subcell_states,
        cell_to_bgcell=native.cell_to_bgcell, point_ids=native.point_ids,
        point_pointers=native.point_pointers, physical_points=native.physical_points,
        reference_points=native.reference_points,
        cell_facet_states=native.cell_facet_states,
        boundary_facet_to_bgcell=native.boundary.facet_to_bgcell,
        boundary_point_ids=native.boundary.point_ids,
        boundary_point_pointers=native.boundary.point_pointers,
        boundary_physical_points=native.boundary.physical_points,
        boundary_reference_points=native.boundary.reference_points,
        boundary_normals=native.boundary.normals)
    for field in (:bg_states, :subcell_states, :cell_to_bgcell, :cell_subcell_ptrs,
                  :point_ids, :point_pointers, :physical_points, :reference_points,
                  :cell_facet_states)
        @test getfield(native, field) == getfield(checked, field)
    end
    @test native.boundary.facet_to_bgcell == checked.boundary.facet_to_bgcell
    @test native.boundary.point_ids == checked.boundary.point_ids
    @test native.boundary.physical_points == checked.boundary.physical_points
    @test native.boundary.normals == checked.boundary.normals
end
