using Test
using Main.EvolvingDomains
using Main.EvolvingDomains.Geometric: _GridMetrics, _getInterpolator
using Main.EvolvingDomains.Kinematic: WENO5Cache
using StaticArrays: SVector
using Main.EvolvingDomains: Grid, Discretization, Field

function field_fixture3D(n)
    x = collect(range(-1.0, 1.0; length=n))
    y = collect(range(-0.8, 0.8; length=n))
    z = collect(range(-0.6, 0.6; length=n))
    g = Grid((x, y, z))
    u = [sin(1.2a) + cos(1.1b) + sin(1.3c) for a in x, b in y, c in z]
    return Field(Discretization(Dict(:init => g)), Dict(:init => u)), x, y, z
end

@testset "3D Cartesian field and interpolation" begin
    x, y, z = 0.0:1.0:2.0, -1.0:1.0:1.0, 2.0:1.0:4.0
    info = _GridMetrics((first(x), first(y), first(z)), (1.0, 1.0, 1.0), (3, 3, 3))
    vals = [2 + 3a - 4b + 0.5c for a in x, b in y, c in z]
    itp = _getInterpolator(vec(vals), info)
    @test itp(0.4, -0.3, 2.7) ≈ 2 + 3(0.4) - 4(-0.3) + 0.5(2.7)
    @test itp(-Inf, -0.3, 2.7) ≈ itp(0.0, -0.3, 2.7)
    @test isnan(itp(NaN, 0.0, 0.0))
    @test_throws DimensionMismatch _getInterpolator(zeros(26), info)
    @test_throws ArgumentError _getInterpolator(zeros(27), _GridMetrics((0.,0.,0.),
        (1.,1.,0.), (3,3,3)))
end

@testset "3D WENO5 field evolution" begin
    f, x, y, z = field_fixture3D(15)
    f.data[:init] .= 2.75
    constant_velocity = fill(SVector(0.1, -0.2, 0.3), length(f.data[:init]))
    eulerianTransport!(f, constant_velocity, 0.01; cache=WENO5Cache())
    @test f.data[:init] ≈ fill(2.75, size(f.data[:init]))

    f, x, y, z = field_fixture3D(15)
    original = copy(f.data[:init])
    cache = WENO5Cache()
    zero_velocity = fill(SVector(0.0, 0.0, 0.0), length(original))
    @test eulerianTransport!(f, zero_velocity, 0.0; cache) === f
    @test f.data[:init] == original
    @test cache.rhs === nothing
    @test eulerianTransport!(f, zero_velocity, 0.01; cache) === f
    @test maximum(abs.(f.data[:init] .- original)) <= 8eps(maximum(abs, original))
    @test cache.rhs !== nothing

    # A constant z translation checks the third derivative and SSP-RK3 path.
    f, x, y, z = field_fixture3D(21)
    initial = copy(f.data[:init])
    velocity = fill(SVector(0.0, 0.0, 0.2), length(initial))
    eulerianTransport!(f, velocity, 0.02; cache=WENO5Cache())
    expected = [sin(1.2a) + cos(1.1b) + sin(1.3(c - 0.004)) for a in x, b in y, c in z]
    @test maximum(abs.(f.data[:init][3:end-2, 3:end-2, 3:end-2] .-
        expected[3:end-2, 3:end-2, 3:end-2])) < 2e-4

    # All three velocity components contribute the upwind derivatives.
    evolved = Array{Float64,3}[]
    exact_fields = Array{Float64,3}[]
    for n in (15, 29)
        f, x, y, z = field_fixture3D(n)
        u0 = copy(f.data[:init])
        v = SVector(0.08, -0.06, 0.05)
        eulerianTransport!(f, fill(v, length(u0)), 0.01; cache=WENO5Cache())
        exact = [sin(1.2(a - 0.0008)) + cos(1.1(b + 0.0006)) +
                 sin(1.3(c - 0.0005)) for a in x, b in y, c in z]
        interior = (3:n-2, 3:n-2, 3:n-2)
        err = maximum(abs.(f.data[:init][interior...] .- exact[interior...]))
        push!(evolved, copy(f.data[:init]))
        push!(exact_fields, exact)
        @test isfinite(err)
        @test err < (n == 15 ? 3e-4 : 1.5e-4)
    end
    # Compare both resolutions at identical physical nodes on the same interior.
    coarse_nodes = 3:13
    fine_nodes = 2 .* (coarse_nodes .- 1) .+ 1
    coarse_error = maximum(abs.(evolved[1][coarse_nodes, coarse_nodes, coarse_nodes] .-
        exact_fields[1][coarse_nodes, coarse_nodes, coarse_nodes]))
    fine_error = maximum(abs.(evolved[2][fine_nodes, fine_nodes, fine_nodes] .-
        exact_fields[2][fine_nodes, fine_nodes, fine_nodes]))
    @test fine_error < coarse_error

    f, _, _, _ = field_fixture3D(7)
    axes = f.discretization.levels[:init].axes
    vel = fill(SVector(0.0, 0.0, 0.0), length(f.data[:init]))
    @test_throws ArgumentError eulerianTransport!(f, vel, -0.1)
    @test_throws ArgumentError eulerianTransport!(f, vel, Inf)
    @test_throws DimensionMismatch eulerianTransport!(f, vel[1:end-1], 0.0)
    bad = copy(vel); bad[1] = SVector(NaN, 0.0, 0.0)
    @test_throws ArgumentError eulerianTransport!(f, bad, 0.0)
    before = copy(f.data[:init]); f.data[:init][1] = NaN
    @test_throws ArgumentError eulerianTransport!(f, vel, 0.0)
    @test isnan(f.data[:init][1])
    @test before[2:end] == f.data[:init][2:end]
end
