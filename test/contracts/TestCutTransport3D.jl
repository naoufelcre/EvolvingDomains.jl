using Test
using Main.EvolvingDomains
using Main.EvolvingDomains.Kinematic: TransportMap, lagrangianTransport!
using StaticArrays: SVector

function _cut3(grid, state_at)
    nx, ny, nz = ntuple(d -> length(grid.axes[d]) - 1, 3)
    states = Int8[state_at(i, j, k) for k in 1:nz for j in 1:ny for i in 1:nx]
    return Cut(grid, :init; bg_states=states)
end

@testset "Conservative Cut transport in 3D" begin
    axes = ntuple(_ -> collect(range(0.0, 1.0; length=7)), 3)
    grid3 = Grid(axes)
    d = Discretization(Dict(:init => grid3))
    levelset = Field(d, Dict(:init => fill(NaN, 7, 7, 7)))
    full = _cut3(grid3, (_, _, _) -> Int8(-1))
    empty = _cut3(grid3, (_, _, _) -> Int8(1))
    zero_velocity = fill(SVector(0.0, 0.0, 0.0), 343)
    source = [1.0 + 0.01i for i in 1:343]
    target = zeros(343)

    identity = TransportMap(full, full, grid3, zero_velocity, 0.2; cache_weights=true)
    @test identity isa EvolvingDomains.Kinematic.SemiLagrangian.TransportMap{3,8}
    @test identity.is_identity
    @test lagrangianTransport!(target, source, identity) === target
    @test target == source

    velocity_z = fill(SVector(0.0, 0.0, 0.12), 343)
    moving_z_cached = TransportMap(full, full, grid3, velocity_z, 0.2; cache_weights=true)
    moving_z_plain = TransportMap(full, full, grid3, velocity_z, 0.2)
    cached = zeros(343); plain = zeros(343)
    lagrangianTransport!(cached, source, moving_z_cached)
    lagrangianTransport!(plain, source, moving_z_plain)
    @test cached == plain
    @test sum(cached) ≈ sum(source) rtol=2e-14
    @test cached != source
    @test any(ray -> ray[3] != 0.0, moving_z_cached.leakage_rays)

    velocity_xyz = fill(SVector(0.04, -0.03, 0.05), 343)
    xyzmap = TransportMap(full, full, grid3, velocity_xyz, 0.2; cache_weights=true)
    xyz = zeros(343)
    lagrangianTransport!(xyz, source, xyzmap)
    @test sum(xyz) ≈ sum(source) rtol=2e-14
    @test xyz != source
    @test any(ray -> ray[1] != 0 && ray[2] != 0 && ray[3] != 0, xyzmap.leakage_rays)

    stencil_support = falses(343)
    for k in 4:5, j in 4:5, i in 4:5
        stencil_support[i + (j - 1) * 7 + (k - 1) * 49] = true
    end
    stencil_indices, stencil_weights =
        EvolvingDomains.Kinematic.SemiLagrangian._computeConservativeWeights(
            SVector(0.5, 0.5, 0.5), moving_z_cached.grid_meta, stencil_support)
    @test all(>=(0.0), stencil_weights)
    @test sum(stencil_weights) ≈ 1.0 rtol=2eps(Float64)
    @test all(stencil_weights[q] == 0.0 for q in eachindex(stencil_weights)
        if !stencil_support[stencil_indices[q]])

    # A distant/disconnected moving pull is rejected rather than losing mass or
    # inventing a global fallback destination for its leakage.
    near_source = _cut3(grid3, (i, _, _) -> i == 1 ? Int8(-1) : Int8(1))
    far_target = _cut3(grid3, (i, _, _) -> i == 6 ? Int8(-1) : Int8(1))
    far_velocity = fill(SVector(0.2, 0.0, 0.0), 343)
    @test_throws ArgumentError TransportMap(near_source, far_target, grid3,
        far_velocity, 0.1)

    src_cut = _cut3(grid3, (i, _, _) -> i <= 3 ? Int8(-1) : Int8(1))
    dst_cut = _cut3(grid3, (i, _, _) -> i <= 4 ? Int8(-1) : Int8(1))
    plane_velocity = fill(SVector(0.5, -0.03, 0.05), 343)
    plane_map = TransportMap(src_cut, dst_cut, grid3, plane_velocity, 0.2; cache_weights=true)
    srcmask = plane_map.source_mask
    plane_data = [srcmask[i] ? source[i] : 0.0 for i in eachindex(source)]
    plane_target = zeros(343)
    lagrangianTransport!(plane_target, plane_data, plane_map)
    @test sum(plane_target) ≈ sum(plane_data) rtol=2e-14
    @test all(plane_target[.!plane_map.target_mask] .== 0.0)

    # A disconnected support is legal when stationary: only active corners carry values.
    disconnected = _cut3(grid3, (i, j, k) ->
        (i <= 2 && j <= 2 && k <= 2) || (i >= 5 && j >= 5 && k >= 5) ? Int8(-1) : Int8(1))
    disconnected_map = TransportMap(disconnected, disconnected, grid3, zero_velocity, 0.0)
    disc_target = zeros(343)
    lagrangianTransport!(disc_target, source, disconnected_map)
    @test disc_target[disconnected_map.source_mask] == source[disconnected_map.source_mask]
    @test all(disc_target[.!disconnected_map.source_mask] .== 0.0)

    empty_map = TransportMap(empty, empty, grid3, zero_velocity, 0.1)
    empty_target = ones(343)
    lagrangianTransport!(empty_target, source, empty_map)
    @test all(iszero, empty_target)
    @test_throws ArgumentError TransportMap(full, empty, grid3, zero_velocity, 0.1)

    source_field = Field(d, Dict(:init => reshape(copy(source), 7, 7, 7)))
    target_field = Field(d, Dict(:init => zeros(7, 7, 7)))
    @test lagrangianTransport!(target_field, source_field, moving_z_cached) === target_field
    @test sum(target_field.data[:init]) ≈ sum(source) rtol=2e-14
    @test_throws ArgumentError lagrangianTransport!(source_field, source_field, moving_z_cached)
    @test_throws MethodError lagrangianTransport!(target_field, source_field, moving_z_cached; type=:intensive)

    bad_velocity = copy(zero_velocity)
    bad_velocity[1] = SVector(0.0, NaN, 0.0)
    @test_throws ArgumentError TransportMap(full, full, grid3, bad_velocity, 0.1)
    @test_throws DimensionMismatch TransportMap(full, full, grid3, zero_velocity[1:end-1], 0.1)
    @test_throws ArgumentError TransportMap(full, full, grid3, zero_velocity, Inf)
    @test_throws ArgumentError TransportMap(full, full, grid3, zero_velocity, -0.1)
    velocity2d = fill(SVector(0.0, 0.0), 343)
    @test_throws MethodError TransportMap(full, full, grid3, velocity2d, 0.1)
    @test_throws DimensionMismatch lagrangianTransport!(zeros(342), source, identity)
    bad_source = copy(source)
    bad_source[1] = NaN
    @test_throws ArgumentError lagrangianTransport!(zeros(343), bad_source, moving_z_cached)
    wrong_shape = Field(d, Dict(:init => zeros(7, 7, 49)))
    @test_throws DimensionMismatch lagrangianTransport!(wrong_shape, source_field, moving_z_cached)
    wrong_level_cut = Cut(grid3, :other; bg_states=fill(Int8(-1), 216))
    @test_throws ArgumentError TransportMap(full, wrong_level_cut, grid3, zero_velocity, 0.1)

    other_grid = Grid(ntuple(_ -> collect(range(0.0, 1.0; length=7)), 3))
    other_d = Discretization(Dict(:init => other_grid))
    wrong_field = Field(other_d, Dict(:init => zeros(7, 7, 7)))
    @test_throws ArgumentError lagrangianTransport!(wrong_field, source_field, identity)

    shifted_target = deepcopy(target_field)
    shifted_target.discretization.levels[:init].axes[3] .+= 1.0
    fill!(shifted_target.data[:init], -1.0)
    @test_throws ArgumentError lagrangianTransport!(shifted_target, source_field, moving_z_cached)
    @test all(==(-1.0), shifted_target.data[:init])

    grid3.axes[1] .+= 1.0
    @test_throws ArgumentError lagrangianTransport!(target_field, source_field, moving_z_cached)
end
