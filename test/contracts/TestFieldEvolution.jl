using Test
using Main.EvolvingDomains
using Main.EvolvingDomains.Geometric: _GridMetrics, _gridMetrics
using Main.EvolvingDomains.Kinematic: WENO5Cache, weno5Step!
using StaticArrays: SVector

using Main.EvolvingDomains: Grid, Discretization, Field

function field_fixture(; nx=19, ny=17)
    x = collect(range(-1.0, 1.0; length=nx))
    y = collect(range(-0.75, 0.75; length=ny))
    g0, gfine = Grid((x, y)), Grid((x, y))
    d = Discretization(Dict(:init => g0, :fine => gfine))
    u = [sin(2xv) + cos(3yv) for xv in x, yv in y]
    return Field(d, Dict(:init => copy(u), :fine => copy(u))), x, y
end

@testset "Field WENO5 + SSP-RK3 overload" begin
    f, x, y = field_fixture()
    vel = fill(SVector(0.3, -0.2), length(f.data[:init]))
    actual = copy(f.data[:init])
    expected = Vector(vec(copy(actual)))
    # Derive the metrics exactly as the adapter does, so the comparison pins the kernel,
    # not a second spacing rule.
    info = _gridMetrics(f.discretization.levels[:fine])
    weno5Step!(expected, info, vel, 0.001, WENO5Cache())

    cache = WENO5Cache()
    @test eulerianTransport!(f, vel, 0.001; level_id=:fine, cache=cache) === f
    @test f.data[:init] == actual
    @test vec(f.data[:fine]) == expected
    @test cache.rhs !== nothing

end

@testset "Field WENO5 validation and default level" begin
    f, _, _ = field_fixture()
    vel = fill(SVector(0.0, 0.0), length(f.data[:init]))
    before = copy(f.data[:init])
    @test eulerianTransport!(f, vel, 0.0) === f
    @test f.data[:init] == before

    @test_throws DimensionMismatch eulerianTransport!(f, vel[1:end-1], 0.01; level_id=:init)
    badgrid = Grid(([0.0, 0.2, 0.5], [0.0, 0.5, 1.0]))
    bad = Field(Discretization(Dict(:init => badgrid)), Dict(:init => zeros(3, 3)))
    @test_throws ArgumentError eulerianTransport!(bad,
        fill(SVector(0.0, 0.0), 9), 0.01)
    @test_throws MethodError eulerianTransport!(nothing, vel, 0.01)
    # The uniformity gate is scaled by the axis extent, not the local coordinate: an axis
    # spanning symmetric coordinates must not be rejected near its zero crossing.
    for (lo, hi, n) in ((-2.0, 2.0, 321), (-5.0, 5.0, 161), (-100.0, 100.0, 401))
        axis = collect(range(lo, hi; length=n))
        symmetric = Field(Discretization(Dict(:init => Grid((axis, copy(axis))))),
                          Dict(:init => zeros(n, n)))
        @test eulerianTransport!(symmetric, fill(SVector(0.0, 0.0), n * n), 0.0) === symmetric
    end
    @test_throws ArgumentError eulerianTransport!(f, vel, -0.01)
    @test_throws ArgumentError eulerianTransport!(f, vel, Inf)

    invalid_velocity = copy(vel)
    invalid_velocity[1] = SVector(NaN, 0.0)
    @test_throws ArgumentError eulerianTransport!(f, invalid_velocity, 0.0)
    @test f.data[:init] == before

    f.data[:init][1] = NaN
    invalid_values = copy(f.data[:init])
    @test_throws ArgumentError eulerianTransport!(f, vel, 0.0)
    @test isequal(f.data[:init], invalid_values)
end
