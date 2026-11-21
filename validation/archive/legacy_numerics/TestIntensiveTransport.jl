using Test
using EvolvingDomains
using EvolvingDomains.Geometric
using Gridap
using Gridap.TensorValues

function transport_fixture(n=16)
    model = CartesianDiscreteModel((0.0, 1.0, 0.0, 1.0), (n, n))
    info = grid_info(model)
    points = vec(collect(Gridap.Geometry.get_node_coordinates(model)))
    return info, points
end

@testset "Intensive CIP identity and API" begin
    info, points = transport_fixture()
    source = CartesianMeshField([sinpi(p[1]) * cospi(p[2]) for p in points], info)
    target = CartesianMeshField(zeros(length(points)), info)
    velocity = fill(VectorValue(0.0, 0.0), length(points))

    @test advect!(target, source, velocity, .1) === target
    @test target.data == source.data
    @test_throws ArgumentError advect!(target, source, velocity, 0.1; type=:conservative)
    @test_throws ArgumentError advect!(target, source, velocity, 0.1; type=:unknown)
    @test_throws ArgumentError advect!(source, source, velocity, 0.1)
end

@testset "Compressible flow preserves intensive constants" begin
    info, points = transport_fixture()
    a, b, dt = 0.5, 0.25, 0.1
    velocity = [VectorValue(a * (p[1] - 0.5), b * (p[2] - 0.5)) for p in points]
    constant_source = CartesianMeshField(fill(2.5, length(points)), info)
    target = CartesianMeshField(zeros(length(points)), info)

    advect!(target, constant_source, velocity, dt)
    @test target.data ≈ constant_source.data atol=1e-14

    affine_source = CartesianMeshField([1 + 2p[1] - 0.5p[2] for p in points], info)
    advect!(target, affine_source, velocity, dt)
    expected = [1 + 2(p[1] - dt * v[1]) - 0.5(p[2] - dt * v[2])
                for (p, v) in zip(points, velocity)]
    @test target.data ≈ expected atol=2e-14
end

@testset "Anisotropic shifted grid and field velocity" begin
    model = CartesianDiscreteModel((-2.0, 1.0, 3.0, 5.0), (24, 10))
    info = grid_info(model)
    points = vec(collect(Gridap.Geometry.get_node_coordinates(model)))
    dt = 0.1
    velocity_data = [VectorValue(
        0.2 * (p[1] + 0.5), 0.15 * (p[2] - 4.0)) for p in points]
    velocity = CartesianMeshField(velocity_data, info)
    source = CartesianMeshField([1 + 2p[1] - 0.5p[2] for p in points], info)
    target = CartesianMeshField(zeros(length(points)), info)

    advect!(target, source, velocity, dt)
    expected = [1 + 2(p[1] - dt * v[1]) - 0.5(p[2] - dt * v[2])
                for (p, v) in zip(points, velocity_data)]
    @test target.data ≈ expected atol=3e-14
end

@testset "Split CIP carries cross-gradient terms" begin
    info, points = transport_fixture(24)
    a, b, dt = 0.2, 0.15, 0.4
    velocity = [VectorValue(
        a * p[2] * p[1] * (1 - p[1]),
        b * p[1] * p[2] * (1 - p[2]),
    ) for p in points]
    source = CartesianMeshField([1 + 2p[1] + 3p[2] for p in points], info)
    target = CartesianMeshField(zeros(length(points)), info)
    cache = CIPCache(source)

    advect!(target, source, velocity, dt; cache=cache)
    expected = map(points) do p
        y_departure = p[2] - dt * b * p[1] * p[2] * (1 - p[2])
        x_departure = p[1] - dt * a * y_departure * p[1] * (1 - p[1])
        return 1 + 2x_departure + 3y_departure
    end
    @test target.data ≈ expected atol=3e-13

    nx, ny = info.dims
    cx_error = cy_error = cxy_error = 0.0
    for j in 1:ny, i in 1:nx
        k = i + (j - 1) * nx
        x, y = points[k]
        r = x * (1 - x)
        dr = 1 - 2x
        q = y * (1 - y)
        y_departure = y - dt * b * x * q
        Yx = -dt * b * q
        Yy = 1 - dt * b * x * (1 - 2y)
        Yxy = -dt * b * (1 - 2y)
        Xx = 1 - dt * a * (dr * y_departure + r * Yx)
        Xy = -dt * a * r * Yy
        Xxy = -dt * a * (dr * Yy + r * Yxy)

        cx_error = max(cx_error, abs(cache.cx[k] - (2Xx + 3Yx)))
        cy_error = max(cy_error, abs(cache.cy[k] - (2Xy + 3Yy)))
        cxy_error = max(cxy_error, abs(cache.cxy[k] - (2Xxy + 3Yxy)))
    end
    @test cx_error < 3e-12
    @test cy_error < 3e-12
    @test cxy_error < 3e-11
end

@testset "CIP cache detects external field changes" begin
    info, points = transport_fixture()
    velocity = [VectorValue(0.2 * p[1] * (1 - p[1]), 0.0) for p in points]
    source = CartesianMeshField([sinpi(p[1]) + 0.25cospi(p[2]) for p in points], info)
    transported = CartesianMeshField(zeros(length(points)), info)
    cache = CIPCache(source)

    advect!(transported, source, velocity, 0.2; cache=cache)
    transported.data .+= [0.1p[2] for p in points]

    cached_target = CartesianMeshField(zeros(length(points)), info)
    stateless_target = CartesianMeshField(zeros(length(points)), info)
    advect!(cached_target, transported, velocity, 0.2; cache=cache)
    advect!(stateless_target, transported, velocity, 0.2)
    @test cached_target.data ≈ stateless_target.data atol=1e-14
end

@testset "CIP cache carries the profile across steps" begin
    info, points = transport_fixture()
    a, dt = 0.4, 0.2
    scale = 1 - a * dt
    velocity = [VectorValue(a * p[1], 0.0) for p in points]
    source = CartesianMeshField([p[1]^2 + 0.5p[2] for p in points], info)
    first = CartesianMeshField(zeros(length(points)), info)
    second = CartesianMeshField(zeros(length(points)), info)
    cache = CIPCache(source)

    advect!(first, source, velocity, dt; cache=cache)
    advect!(second, first, velocity, dt; cache=cache)

    expected = [scale^4 * p[1]^2 + 0.5p[2] for p in points]
    expected_cx = [2scale^4 * p[1] for p in points]
    @test second.data ≈ expected atol=3e-14
    @test cache.cx ≈ expected_cx atol=2e-13
    @test cache.cy ≈ fill(0.5, length(points)) atol=2e-13
    @test cache.cxy ≈ zeros(length(points)) atol=2e-12
end

@testset "Multi-cell departures and convergence" begin
    info, points = transport_fixture(32)
    dt = 0.5
    velocity = [VectorValue(
        0.5 * p[1] * (1 - p[1]), 0.5 * p[2] * (1 - p[2])) for p in points]
    @test maximum(dt * v[1] for v in velocity) > info.spacing[1]
    @test maximum(dt * v[2] for v in velocity) > info.spacing[2]

    source = CartesianMeshField([1 + 2p[1] + 3p[2] for p in points], info)
    target = CartesianMeshField(zeros(length(points)), info)
    advect!(target, source, velocity, dt)
    expected = [1 + 2(p[1] - dt * v[1]) + 3(p[2] - dt * v[2])
                for (p, v) in zip(points, velocity)]
    @test target.data ≈ expected atol=2e-14

    cubic_source = CartesianMeshField([p[1]^3 + 0.2p[2] for p in points], info)
    for direction in (-1.0, 1.0)
        cubic_velocity = [VectorValue(
            direction * 0.5 * p[1] * (1 - p[1]), 0.0) for p in points]
        cubic_cache = CIPCache(cubic_source)
        cubic_cache.cx .= [3p[1]^2 for p in points]
        cubic_cache.cy .= 0.2
        cubic_cache.cxy .= 0.0
        advect!(target, cubic_source, cubic_velocity, dt; cache=cubic_cache)
        cubic_expected = [(p[1] - dt * v[1])^3 + 0.2p[2]
                          for (p, v) in zip(points, cubic_velocity)]
        @test target.data ≈ cubic_expected atol=3e-14
    end

    function refinement_error(n)
        grid, nodes = transport_fixture(n)
        initial(p) = sinpi(2p[1]) * cospi(p[2])
        input = CartesianMeshField(initial.(nodes), grid)
        output = CartesianMeshField(zeros(length(nodes)), grid)
        flow = [VectorValue(0.3 * p[1] * (1 - p[1]), 0.0) for p in nodes]
        step = 0.3
        advect!(output, input, flow, step)
        exact = [initial((p[1] - step * v[1], p[2])) for (p, v) in zip(nodes, flow)]
        return maximum(abs.(output.data .- exact))
    end

    e16, e32, e64 = refinement_error.((16, 32, 64))
    @test e32 < e16 / 5
    @test e64 < e32 / 5
end

@testset "CIP rejects unsupported inputs before writing" begin
    info, points = transport_fixture(8)
    source = CartesianMeshField([p[1] + p[2] for p in points], info)
    sentinel = fill(-7.0, length(points))
    target = CartesianMeshField(copy(sentinel), info)

    inflow = fill(VectorValue(1.0, 0.0), length(points))
    @test_throws DomainError advect!(target, source, inflow, 0.1)
    @test target.data == sentinel

    invalid = copy(inflow)
    invalid[1] = VectorValue(NaN, 0.0)
    @test_throws ArgumentError advect!(target, source, invalid, 0.1)
    @test target.data == sentinel
    @test_throws ArgumentError advect!(target, source, inflow, -0.1)
    @test_throws ArgumentError advect!(target, source, inflow, Inf)

    late_inflow = fill(VectorValue(0.0, 0.0), length(points))
    late_inflow[end] = VectorValue(0.0, -1.0)
    @test_throws DomainError advect!(target, source, late_inflow, 0.1)
    @test target.data == sentinel

    folded_info, folded_points = transport_fixture(128)
    folded_source = CartesianMeshField([sinpi(2p[1]) for p in folded_points], folded_info)
    folded_target = CartesianMeshField(fill(-3.0, length(folded_points)), folded_info)
    folded_velocity = [VectorValue(
        p[1] * (1 - p[1]) * sin(16π * p[1]), 0.0) for p in folded_points]
    folded_departures = [p[1] - 0.5v[1] for (p, v) in zip(folded_points, folded_velocity)]
    @test all(x -> 0 <= x <= 1, folded_departures)
    @test_throws DomainError advect!(folded_target, folded_source, folded_velocity, 0.5)
    @test folded_target.data == fill(-3.0, length(folded_points))

    y_folded_velocity = [VectorValue(
        0.0, p[2] * (1 - p[2]) * sin(16π * p[2])) for p in folded_points]
    y_folded_departures = [p[2] - 0.5v[2]
                           for (p, v) in zip(folded_points, y_folded_velocity)]
    @test all(y -> 0 <= y <= 1, y_folded_departures)
    @test_throws DomainError advect!(folded_target, folded_source, y_folded_velocity, 0.5)
    @test folded_target.data == fill(-3.0, length(folded_points))

    cache = CIPCache(source)
    cache.cx[1] = Inf
    recovered = CartesianMeshField(zeros(length(points)), info)
    reference = CartesianMeshField(zeros(length(points)), info)
    safe_velocity = [VectorValue(0.1 * p[1] * (1 - p[1]), 0.0) for p in points]
    advect!(recovered, source, safe_velocity, 0.1; cache=cache)
    advect!(reference, source, safe_velocity, 0.1)
    @test recovered.data ≈ reference.data atol=1e-14

    malformed = CIPCache(source)
    resize!(malformed.cx, 1)
    @test_throws DimensionMismatch advect!(target, source, safe_velocity, 0.1; cache=malformed)

    shifted_info = CartesianGridInfo((1.0, 0.0), info.spacing, info.dims, info.cells)
    shifted_cache = CIPCache(CartesianMeshField(copy(source.data), shifted_info))
    @test_throws ArgumentError advect!(target, source, safe_velocity, 0.1; cache=shifted_cache)

    old_snapshot = copy(cache.snapshot)
    changed_source = CartesianMeshField(source.data .+ 1.0, info)
    @test_throws DomainError advect!(target, changed_source, inflow, 0.1; cache=cache)
    @test cache.snapshot == old_snapshot

    huge_source = CartesianMeshField(fill(floatmax(Float64), length(points)), info)
    huge_cache = CIPCache(huge_source)
    @test all(iszero, huge_cache.cx)
    @test all(iszero, huge_cache.cy)
    @test all(iszero, huge_cache.cxy)
    advect!(target, huge_source, safe_velocity, 0.1; cache=huge_cache)
    @test target.data == huge_source.data

    zero_source = CartesianMeshField(zeros(length(points)), info)
    overflow_cache = CIPCache(zero_source)
    overflow_cache.cy .= floatmax(Float64)
    overflow_target = CartesianMeshField(copy(sentinel), info)
    overflow_velocity = [VectorValue(
        0.0, 0.5 * p[2] * (1 - p[2])) for p in points]
    @test_throws ErrorException advect!(
        overflow_target, zero_source, overflow_velocity, 0.1; cache=overflow_cache)
    @test overflow_target.data == sentinel
    @test isnan(overflow_cache.snapshot[1])
    advect!(overflow_target, zero_source, overflow_velocity, 0.1; cache=overflow_cache)
    @test all(iszero, overflow_target.data)

    large_info = CartesianGridInfo((1e16, 0.0), (100.0, 0.125), (9, 9), (8, 8))
    large_source = CartesianMeshField(zeros(81), large_info)
    large_target = CartesianMeshField(fill(-2.0, 81), large_info)
    large_velocity = fill(VectorValue(100.0, 0.0), 81)
    @test_throws DomainError advect!(large_target, large_source, large_velocity, 1.0)
    @test large_target.data == fill(-2.0, 81)

    duplicate_info = CartesianGridInfo((1e16, 0.0), (1.0, 0.125), (9, 9), (8, 8))
    duplicate_source = CartesianMeshField(zeros(81), duplicate_info)
    duplicate_target = CartesianMeshField(fill(-2.0, 81), duplicate_info)
    @test_throws ArgumentError advect!(
        duplicate_target, duplicate_source, fill(VectorValue(0.0, 0.0), 81), 0.1)
    @test duplicate_target.data == fill(-2.0, 81)
end


@testset "Frozen velocity map API" begin
    info, points = transport_fixture(8)
    grid = CartesianDiscreteModel((0.0, 1.0, 0.0, 1.0), (8, 8))
    geom = EvolvingDiscreteGeometry(fill(-1.0, length(points)), grid)
    ensure_cut!(geom)
    set_levelset!(geom, copy(current_levelset(geom)))

    sampled = fill(VectorValue(0.0, 0.0), length(points))
    sampled_field = CartesianMeshField(sampled, info)
    @test TransportMap(geom, sampled, 0.1).is_identity
    @test TransportMap(geom, sampled_field, 0.1).is_identity
    @test TransportMap(geom, StaticFunctionVelocity(_ -> (0.0, 0.0)), 0.1).is_identity
    @test_throws ArgumentError TransportMap(
        geom, TimeDependentVelocity((_, t) -> VectorValue(t, 0.0)), 0.1)
    @test advance!(geom, sampled_field, 0.0) === geom
    shifted_info = CartesianGridInfo((1.0, 0.0), info.spacing, info.dims, info.cells)
    @test_throws ArgumentError advance!(geom, CartesianMeshField(sampled, shifted_info), 0.0)
end
