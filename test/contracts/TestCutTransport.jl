using Test
using Main.EvolvingDomains
using Main.EvolvingDomains.Kinematic: TransportMap, lagrangianTransport!
using StaticArrays: SVector

@testset "Transport maps from explicit Cuts" begin
    grid2d = Grid((range(0.0, 1.0; length=5), range(0.0, 1.0; length=5)))
    d = Discretization(Dict(:init => grid2d))
    levelset = Field(d, Dict(:init => fill(-1.0, 5, 5)))
    cut = Cut(grid2d, :init; bg_states=fill(Int8(-1), 16))
    zero_velocity = fill(SVector(0.0, 0.0), 25)

    fill!(levelset.data[:init], NaN)
    map = TransportMap(cut, cut, grid2d, zero_velocity, 0.1)
    source = collect(1.0:25.0)
    target = similar(source)
    lagrangianTransport!(target, source, map)
    @test map.is_identity
    @test target == source

    source_field = Field(d, Dict(:init => reshape(copy(source), 5, 5)))
    target_field = Field(d, Dict(:init => fill(0.0, 5, 5)))
    @test lagrangianTransport!(target_field, source_field, map) === target_field
    @test target_field.data[:init] == source_field.data[:init]

    shifted_velocity = fill(SVector(0.2, 0.0), 25)
    moving_map = TransportMap(cut, cut, grid2d, shifted_velocity, 0.1)
    lagrangianTransport!(target, source, moving_map)
    @test sum(target) ≈ sum(source) rtol=5e-15
    lagrangianTransport!(target_field, source_field, moving_map)
    @test sum(target_field.data[:init]) ≈ sum(source_field.data[:init]) rtol=5e-15
    @test_throws ArgumentError lagrangianTransport!(source_field, source_field, map)

    other_grid = Grid((range(0.0, 1.0; length=5), range(0.0, 1.0; length=5)))
    other_cut = Cut(other_grid, :init; bg_states=fill(Int8(-1), 16))
    other_level = Field(Discretization(Dict(:init => other_grid)), Dict(:init => fill(-1.0, 5, 5)))
    other_target = Field(Discretization(Dict(:init => other_grid)), Dict(:init => fill(0.0, 5, 5)))
    @test_throws ArgumentError lagrangianTransport!(other_target, source_field, map)

    multi_level_d = Discretization(Dict(:init => grid2d, :other => grid2d))
    level_mismatch_target = Field(multi_level_d, Dict(:other => fill(0.0, 5, 5)))
    level_mismatch_source = Field(multi_level_d, Dict(:other => copy(source_field.data[:init])))
    @test_throws ArgumentError lagrangianTransport!(level_mismatch_target, level_mismatch_source, map)
    @test_throws ArgumentError TransportMap(cut, other_cut, grid2d, zero_velocity, 0.1)
    @test_throws ArgumentError TransportMap(cut, cut, other_grid, zero_velocity, 0.1)

    other_level_id = Cut(grid2d, :other; bg_states=fill(Int8(-1), 16))
    @test_throws ArgumentError TransportMap(cut, other_level_id, grid2d, zero_velocity, 0.1)

    nonuniform = Grid(([0.0, 0.2, 0.5, 0.75, 1.0], range(0.0, 1.0; length=5)))
    nonuniform_cut = Cut(nonuniform, :init; bg_states=fill(Int8(-1), 16))
    nonuniform_field = Field(Discretization(Dict(:init => nonuniform)), Dict(:init => fill(-1.0, 5, 5)))
    @test_throws ArgumentError TransportMap(nonuniform_cut, nonuniform_cut, nonuniform,
        fill(SVector(0.0, 0.0), 25), 0.1)

    @test_throws ArgumentError TransportMap(cut, cut, grid2d, zero_velocity, -0.1)
    @test_throws DimensionMismatch TransportMap(cut, cut, grid2d, zero_velocity[1:end-1], 0.1)
    invalid_velocity = copy(zero_velocity)
    invalid_velocity[1] = SVector(NaN, 0.0)
    @test_throws ArgumentError TransportMap(cut, cut, grid2d, invalid_velocity, 0.1)
    wrong_shape = Field(d, Dict(:init => zeros(25, 1)))
    @test_throws DimensionMismatch lagrangianTransport!(wrong_shape, source_field, map)
    @test_throws MethodError lagrangianTransport!(target_field, source_field, map; type=:intensive)
    source_field.data[:init][1] = NaN
    @test_throws ArgumentError lagrangianTransport!(target_field, source_field, map)
    source_field.data[:init][1] = source[1]
    grid2d.axes[1] .+= 1.0
    @test_throws ArgumentError lagrangianTransport!(target_field, source_field, map)
end

@testset "conservative operator: paper invariants and the boundary limit" begin
    # Lentine, Grétarsson & Fedkiw (2011) §2: every source must contribute exactly what it
    # has (row sums = 1), and the interpolation is a partition of unity per target wherever
    # the stencil is not truncated. The second property is what makes a constant field stay
    # constant; it therefore holds in the interior and fails in the boundary columns, where
    # the four-wide stencil runs off the grid (the paper's ghost-cell case).
    n = 17
    h = 1.0 / (n - 1)
    ax = collect(range(0.0, 1.0; length = n))
    grid = EvolvingDomains.Grid((ax, copy(ax)))
    level = Field(Discretization(Dict(:init => grid)), Dict(:init => fill(-1.0, n, n)))
    v = vec([SVector(0.5, 0.0) for x in ax, y in ax])
    map = TransportMap(cut(level, :init), cut(level, :init), grid, v, 0.4h / 0.5)

    src = ones(n * n)
    dst = zeros(n * n)
    lagrangianTransport!(dst, src, map)

    @test sum(dst) ≈ n^2 rtol = 1e-12                    # total conserved
    M = reshape(dst, n, n)
    # Partition of unity to machine precision from seven nodes inward (measured 2.2e-16).
    @test all(isapprox(M[i, j], 1.0; atol = 1e-12) for i in 7:n-6, j in 7:n-6)
    # Along an axis the layer thins fast: 3.3e-3 three nodes in, 4.4e-8 five nodes in.
    @test abs(M[3, 9] - 1) > 1e-6
    @test abs(M[5, 9] - 1) < 1e-6
    # Known limit, pinned so that fixing it is a deliberate act: the boundary columns are
    # not a partition of unity (measured 0.713 at inflow, 1.403 at outflow; and 0.884, then
    # exactly 1, at the two nodes inside — the band-wide far corner is what supplies them).
    @test M[1, 9] < 0.9
    @test M[end, 9] > 1.1
end
