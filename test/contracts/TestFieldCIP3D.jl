using Test
using Main.EvolvingDomains
using Main.EvolvingDomains.Geometric: _GridMetrics
using Main.EvolvingDomains.Kinematic.CIP: CIPCache, cipLagrangianTransport!
using StaticArrays: SVector

function cip3d_fixture(n=12; bounds=((0.0,1.0),(0.0,1.0),(0.0,1.0)))
    axes = ntuple(d -> collect(range(bounds[d]...; length=n)), 3)
    g = EvolvingDomains.Grid(axes)
    d = EvolvingDomains.Discretization(Dict(:init => g))
    points = vec([(x,y,z) for x in axes[1], y in axes[2], z in axes[3]])
    info = _GridMetrics(ntuple(d -> bounds[d][1], 3),
        ntuple(d -> (bounds[d][2]-bounds[d][1])/(n-1), 3),
        ntuple(_ -> n, 3))
    return g, d, points, info
end
cip3d_field(d, axes, f) = EvolvingDomains.Field(d, Dict(:init =>
    [f(x,y,z) for x in axes[1], y in axes[2], z in axes[3]]))

@testset "3D Field CIP affine motion and cache" begin
    g, d, points, info = cip3d_fixture(12)
    source = cip3d_field(d, g.axes, (x,y,z) -> 1 + 2x - 0.5y + 3z)
    target = cip3d_field(d, g.axes, (_,_,_) -> -1.0)
    velocity = [all(0 < c < 1 for c in p) ? SVector(0.035,-0.02,0.025) :
                SVector(0.0,0.0,0.0) for p in points]
    dt = 0.2
    cache = CIPCache(vec(copy(source.data[:init])), info)
    @test lagrangianTransport!(target, source, velocity, dt; cache=cache) === target
    got = target.data[:init]
    interior = 4:9
    expected = [begin
        index = i + (j-1)*12 + (k-1)*144
        1 + 2(x - dt*velocity[index][1]) - 0.5(y - dt*velocity[index][2]) +
            3(z - dt*velocity[index][3])
    end for (i,x) in enumerate(g.axes[1]), (j,y) in enumerate(g.axes[2]),
        (k,z) in enumerate(g.axes[3])]
    @test got[interior,interior,interior] ≈ expected[interior,interior,interior] atol=5e-13
    @test length(cache.cxyz) == length(points)

    stateless = cip3d_field(d, g.axes, (_,_,_) -> 0.0)
    lagrangianTransport!(stateless, source, velocity, dt)
    @test stateless.data[:init] ≈ got atol=1e-14
    fill!(target.data[:init], -7.0)
    @test lagrangianTransport!(target, source, velocity, 0.0) === target
    @test target.data[:init] == source.data[:init]
    zero_velocity = fill(SVector(0.0,0.0,0.0), length(points))
    target.data[:init] .= -3.0
    lagrangianTransport!(target, source, zero_velocity, dt; cache=cache)
    @test target.data[:init] == source.data[:init]
end

@testset "3D CIP profile and split oracles" begin
    g, d, points, info = cip3d_fixture(15)
    values = [x*y*z for (x,y,z) in points]
    source = copy(values)
    cache = CIPCache(source, info)
    velocity = [SVector(x*y*z,0.0,0.0) for (x,y,z) in points]
    EvolvingDomains.Kinematic.CIP._sweep3!(cache, velocity, 0.05, 1, info)
    xyz_error = 0.0
    for k in 4:12, j in 4:12, i in 4:12
        index = i + (j-1)*15 + (k-1)*225
        xyz_error = max(xyz_error, abs(cache.tmp_cxyz[index] -
            (1 - 4*0.05*points[index][2]*points[index][3])))
    end
    @test xyz_error < 2e-11

    target = zeros(length(values))
    zflow = [SVector(0.0,0.0,0.04*z*(1-z)) for (_,_,z) in points]
    cipLagrangianTransport!(target, source, zflow, info, 0.2; cache=CIPCache(source, info))
    expected_z = [x*y*(z - 0.2*0.04*z*(1-z)) for (x,y,z) in points]
    @test target[4:12] ≈ expected_z[4:12] atol=2e-13

    xyz_flow = [SVector(x*(1-x)*y*z, y*(1-y)*x*z, z*(1-z)*x*y)
                for (x,y,z) in points]
    xyz_target = zeros(length(values)); xyz_reference = zeros(length(values))
    cipLagrangianTransport!(xyz_target, source, xyz_flow, info, 0.02; cache=CIPCache(source, info))
    cipLagrangianTransport!(xyz_reference, source, xyz_flow, info, 0.02)
    @test all(isfinite, xyz_target)
    @test xyz_target ≈ xyz_reference atol=1e-14
    @test xyz_target != values

    polynomial_flow = [SVector(x*y*z, x*y*z, x*y*z) for (x,y,z) in points]
    polynomial_target = zeros(length(values))
    cipLagrangianTransport!(polynomial_target, source, polynomial_flow, info, 0.02;
        cache=CIPCache(source, info))
    split_exact = map(points) do (x,y,z)
        z1 = z - 0.02*x*y*z
        y1 = y - 0.02*x*y*z1
        x1 = x - 0.02*x*y1*z1
        x1*y1*z1
    end
    interior = [i + (j-1)*15 + (k-1)*225 for k in 3:13, j in 3:13, i in 3:13]
    @test maximum(abs.(polynomial_target[interior] .- split_exact[interior])) < 3e-11
end

@testset "3D CIP refinement and cache parity" begin
    function refinement_error(n)
        _, _, points, info = cip3d_fixture(n)
        initial(x,y,z) = sinpi(2x)*cospi(y)*sinpi(z)
        source = [initial(x,y,z) for (x,y,z) in points]
        flow = [SVector(0.12*x*(1-x),0.0,0.0) for (x,_,_) in points]
        output = zeros(length(points)); reference = zeros(length(points))
        cache = CIPCache(source, info)
        cipLagrangianTransport!(output, source, flow, info, 0.2; cache=cache)
        cipLagrangianTransport!(reference, source, flow, info, 0.2)
        @test output ≈ reference atol=1e-14
        source .+= 0.25
        cipLagrangianTransport!(output, source, flow, info, 0.2; cache=cache)
        cipLagrangianTransport!(reference, source, flow, info, 0.2)
        @test output ≈ reference atol=1e-14
        exact = [initial(x - 0.2v[1], y, z) + 0.25 for ((x,y,z),v) in zip(points,flow)]
        ids = [i + (j-1)*n + (k-1)*n*n for k in 3:n-2, j in 3:n-2, i in 3:n-2]
        maximum(abs.(output[ids] .- exact[ids]))
    end
    e12, e24, e48 = refinement_error.((12,24,48))
    @test e24 < e12/2
    @test e48 < e24/2
end

@testset "3D CIP validation boundaries" begin
    g, d, points, info = cip3d_fixture(8)
    source = [sinpi(x)*sinpi(y)*sinpi(z) for (x,y,z) in points]
    target = zeros(length(source)); velocity = fill(SVector(0.0,0.0,0.0), length(source))
    cipLagrangianTransport!(target, source, velocity, info, 0.01)
    @test_throws ArgumentError cipLagrangianTransport!(source, source, velocity, info, 0.01)
    @test_throws DimensionMismatch cipLagrangianTransport!(target, source, velocity[1:end-1], info, 0.01)
    bad = copy(source); bad[1] = NaN
    @test_throws ArgumentError cipLagrangianTransport!(target, bad, velocity, info, 0.01)
    nonfinite = copy(velocity); nonfinite[1] = SVector(0.0,NaN,0.0)
    @test_throws ArgumentError cipLagrangianTransport!(target, source, nonfinite, info, 0.01)
    @test_throws ArgumentError cipLagrangianTransport!(target, source, velocity, info, -0.01)
    other = _GridMetrics((0.,0.,0.), (1.,1.,1.), (8,8,8))
    @test_throws ArgumentError cipLagrangianTransport!(target, source, velocity, other, 0.01;
        cache=CIPCache(source, info))

end
