using Test
using Main.EvolvingDomains
using Main.EvolvingDomains.Geometric: _GridMetrics
using StaticArrays: SVector


function fieldcip_fixture(nx=16, ny=12; xlim=(-1.0, 2.0), ylim=(0.0, 3.0))
    axes = (collect(range(xlim...; length=nx)), collect(range(ylim...; length=ny)))
    g = EvolvingDomains.Grid(axes)
    d = EvolvingDomains.Discretization(Dict(:init => g))
    points = [(x, y) for y in axes[2] for x in axes[1]]
    source = EvolvingDomains.Field(d, Dict(:init => [1 + 2x - 0.5y for x in axes[1], y in axes[2]]))
    target = EvolvingDomains.Field(d, Dict(:init => zeros(nx, ny)))
    return g, points, source, target
end

@testset "Field CIP full-grid adapter" begin
    g, points, source, target = fieldcip_fixture()
    source_data = source.data[:init]
    @test vec(source_data) isa Vector{Float64}
    @test Base.mightalias(vec(source_data), source_data)
    velocity = [x == first(g.axes[1]) || x == last(g.axes[1]) ||
                y == first(g.axes[2]) || y == last(g.axes[2]) ?
                SVector(0.0, 0.0) : SVector(0.05, -0.03) for (x, y) in points]
    dt = 0.1

    @test lagrangianTransport!(target, source, velocity, dt) === target
    @test EvolvingDomains.CIPCache(vec(copy(source_data)), g) isa EvolvingDomains.CIPCache{2}
    adapted = copy(target.data[:init])
    info = _GridMetrics((first(g.axes[1]), first(g.axes[2])),
        ((last(g.axes[1]) - first(g.axes[1])) / (length(g.axes[1]) - 1),
         (last(g.axes[2]) - first(g.axes[2])) / (length(g.axes[2]) - 1)),
        (length(g.axes[1]), length(g.axes[2])))
    kernel = zeros(length(source_data))
    EvolvingDomains.Kinematic.CIP._cipLagrangianTransport!(
        kernel, vec(copy(source_data)), velocity, info, dt)
    @test vec(adapted) == kernel

    # Constant-velocity linear translation is exact away from boundary stencils.
    interior = (4:length(g.axes[1])-3, 4:length(g.axes[2])-3)
    expected = [1 + 2(x - dt * 0.05) - 0.5(y + dt * 0.03)
                for x in g.axes[1][interior[1]], y in g.axes[2][interior[2]]]
    @test adapted[interior...] ≈ expected atol=3e-14

    fill!(target.data[:init], -9.0)
    @test lagrangianTransport!(target, source, fill(SVector(0.0, 0.0), length(points)), 0.0) === target
    @test target.data[:init] == source.data[:init]
    @test_throws MethodError lagrangianTransport!(target, source, velocity, dt; type=:conservative)
end

@testset "Field CIP grid and validation boundary" begin
    g, points, source, target = fieldcip_fixture(8, 8)
    velocity = fill(SVector(0.0, 0.0), length(points))
    sentinel = fill(-3.0, size(target.data[:init]))
    target.data[:init] .= sentinel

    other_grid = EvolvingDomains.Grid(g.axes)
    other_d = EvolvingDomains.Discretization(Dict(:init => other_grid))
    other_target = EvolvingDomains.Field(other_d, Dict(:init => zeros(size(sentinel))))
    @test_throws ArgumentError lagrangianTransport!(other_target, source, velocity, 0.1)
    @test other_target.data[:init] == zeros(size(sentinel))

    @test_throws DimensionMismatch lagrangianTransport!(target, source, velocity[1:end-1], 0.1)
    @test target.data[:init] == sentinel

    uneven = EvolvingDomains.Grid((collect(range(0.0, 1.0; length=8)) .^ 2, g.axes[2]))
    uneven_d = EvolvingDomains.Discretization(Dict(:init => uneven))
    uneven_source = EvolvingDomains.Field(uneven_d, Dict(:init => zeros(8, 8)))
    @test_throws ArgumentError lagrangianTransport!(target, uneven_source, velocity, 0.1)
    @test target.data[:init] == sentinel
    wrong_type = EvolvingDomains.Field(source.discretization, Dict(:init => zeros(Float32, 8, 8)))
    @test_throws ArgumentError lagrangianTransport!(target, wrong_type, velocity, 0.1)
end
