using Test
using EvolvingDomains
using EvolvingDomains.Geometric
using Gridap.TensorValues
import Gridap
using Interpolations
using Random

# Reference implementation: the pre-ED9 `get_interpolator`, kept here so the
# hand-written bilinear is checked against the exact Interpolations.jl object it
# replaced (BSpline(Linear()) |> scale |> extrapolate(Flat())).
function reference_interpolator(f::CartesianMeshField)
    nx, ny = f.grid.dims
    data_2d = reshape(f.data, nx, ny)
    x0, y0 = f.grid.origin
    dx, dy = f.grid.spacing
    xaxis = range(x0, step=dx, length=nx)
    yaxis = range(y0, step=dy, length=ny)
    sitp = scale(interpolate(data_2d, BSpline(Linear())), xaxis, yaxis)
    return extrapolate(sitp, Flat())
end

grid_meta(origin, spacing, dims) =
    CartesianGridInfo(origin, spacing, dims, (dims[1] - 1, dims[2] - 1))

const TOL = 1.0e-13

@testset "ED9 clamped bilinear" begin
    Random.seed!(20261003)

    @testset "scalar random interior/boundary/exterior on shifted anisotropic grid" begin
        nx, ny = 7, 5
        info = grid_meta((0.3, -1.7), (0.25, 0.4), (nx, ny))
        f = CartesianMeshField(randn(nx * ny), info)
        new = get_interpolator(f)
        old = reference_interpolator(f)

        @test new.data !== f.data   # snapshot, not a borrowed reference
        @test new.data == f.data

        worst = 0.0
        for _ in 1:50_000
            x = info.origin[1] + (rand() * (nx + 1) - 0.5) * info.spacing[1]
            y = info.origin[2] + (rand() * (ny + 1) - 0.5) * info.spacing[2]
            worst = max(worst, abs(new(x, y) - old(x, y)))
        end
        @test worst < TOL

        # Exact grid nodes: must reproduce the stored value within tolerance.
        for j in 0:(ny - 1), i in 0:(nx - 1)
            x = info.origin[1] + i * info.spacing[1]
            y = info.origin[2] + j * info.spacing[2]
            @test new(x, y) ≈ f.data[i + 1 + j * nx] atol = TOL
            @test new(x, y) ≈ old(x, y) atol = TOL
        end

        # Corners and edges, including points just inside/outside.
        corners = [
            (info.origin[1], info.origin[2]),
            (info.origin[1] + (nx - 1) * info.spacing[1], info.origin[2]),
            (info.origin[1], info.origin[2] + (ny - 1) * info.spacing[2]),
            (info.origin[1] + (nx - 1) * info.spacing[1], info.origin[2] + (ny - 1) * info.spacing[2]),
        ]
        for (x, y) in corners
            @test new(x, y) ≈ old(x, y) atol = TOL
            @test new(x, y) ≈ new(clamp(x, info.origin[1], info.origin[1] + (nx - 1) * info.spacing[1]),
                                  clamp(y, info.origin[2], info.origin[2] + (ny - 1) * info.spacing[2])) atol = TOL
        end

        # Near-ULP perturbations of nodes and of the domain bounds.
        for i in 1:(nx - 1), j in 1:(ny - 1)
            bx = info.origin[1] + i * info.spacing[1]
            by = info.origin[2] + j * info.spacing[2]
            for x in (prevfloat(bx), nextfloat(bx)), y in (prevfloat(by), nextfloat(by))
                @test new(x, y) ≈ old(x, y) atol = TOL
            end
        end
        for x in (prevfloat(info.origin[1]), nextfloat(info.origin[1]),
                  prevfloat(info.origin[1] + (nx - 1) * info.spacing[1]),
                  nextfloat(info.origin[1] + (nx - 1) * info.spacing[1]))
            @test new(x, info.origin[2] + info.spacing[2]) ≈
                  old(x, info.origin[2] + info.spacing[2]) atol = TOL
        end

        # Flat extrapolation: far exterior returns the nearest boundary value and
        # must match the reference, not merely be "bounded".
        for (x, y) in [(info.origin[1] - 10.0, info.origin[2] - 10.0),
                       (info.origin[1] + 10.0 * (nx - 1) * info.spacing[1], info.origin[2] + info.spacing[2]),
                       (info.origin[1] + info.spacing[1], info.origin[2] + 10.0 * (ny - 1) * info.spacing[2])]
            @test new(x, y) ≈ old(x, y) atol = TOL
        end
    end

    @testset "coordinate precision and automatic differentiation" begin
        info = grid_meta((0.0, 0.0), (0.5, 0.5), (3, 3))
        f = CartesianMeshField([1.0 + 2i * 0.5 - 3j * 0.5
                               for j in 0:2 for i in 0:2], info)
        new, old = get_interpolator(f), reference_interpolator(f)
        FD = Gridap.Fields.ForwardDiff
        for point in ([0.2, 0.3], [-0.2, 0.3], [0.2, 1.3])
            @test FD.gradient(z -> new(z...), point) ≈
                  FD.gradient(z -> old(z...), point) atol=TOL
        end
        @test new(big"0.2", big"0.3") isa BigFloat
        @test new(big"0.2", big"0.3") ≈ old(big"0.2", big"0.3") atol=TOL
    end

    @testset "sweep of grid sizes, origins, spacings" begin
        worst = 0.0
        for _ in 1:40
            nx = rand(2:12); ny = rand(2:12)
            info = grid_meta((randn(), randn()),
                             (10.0^(rand() * 3 - 2), 10.0^(rand() * 3 - 2)), (nx, ny))
            f = CartesianMeshField(randn(nx * ny), info)
            new = get_interpolator(f); old = reference_interpolator(f)
            for _ in 1:400
                x = info.origin[1] + (rand() * (nx + 2) - 1.0) * info.spacing[1]
                y = info.origin[2] + (rand() * (ny + 2) - 1.0) * info.spacing[2]
                worst = max(worst, abs(new(x, y) - old(x, y)))
            end
        end
        @test worst < TOL
    end

    @testset "vector-valued data (Gridap VectorValue)" begin
        nx, ny = 6, 9
        info = grid_meta((-2.0, 5.0), (0.5, 0.125), (nx, ny))
        vdata = [VectorValue(randn(), randn()) for _ in 1:(nx * ny)]
        f = CartesianMeshField(vdata, info)
        new = get_interpolator(f)
        old = reference_interpolator(f)

        @test new(0.0, 5.0) isa VectorValue{2,Float64}
        worst = 0.0
        for _ in 1:20_000
            x = info.origin[1] + (rand() * (nx + 1) - 0.5) * info.spacing[1]
            y = info.origin[2] + (rand() * (ny + 1) - 0.5) * info.spacing[2]
            a = new(x, y); b = old(x, y)
            worst = max(worst, abs(a[1] - b[1]), abs(a[2] - b[2]))
        end
        @test worst < TOL

        # Exact nodes for vector data.
        for j in 0:(ny - 1), i in 0:(nx - 1)
            x = info.origin[1] + i * info.spacing[1]
            y = info.origin[2] + j * info.spacing[2]
            @test new(x, y) ≈ vdata[i + 1 + j * nx] atol = TOL
        end
    end

    @testset "snapshot semantics (input mutation does not leak)" begin
        info = grid_meta((0.0, 0.0), (1.0, 1.0), (4, 3))
        f = CartesianMeshField(zeros(12), info)
        new = get_interpolator(f)
        old = reference_interpolator(f)

        @test new.data !== f.data
        @test new.data == f.data
        @test new(1.0, 1.0) == 0.0
        @test old(1.0, 1.0) == 0.0

        # Mutating the input after construction must not change either interpolant.
        f.data[1 + 1 * 4 + 1] = 7.5
        f.data .= 3.0
        @test new(1.0, 1.0) == 0.0
        @test old(1.0, 1.0) == 0.0

        # The private snapshot itself is still evaluable in place.
        new.data[1 + 1 * 4 + 1] = 9.0
        @test new(1.0, 1.0) == 9.0
    end

    @testset "nonfinite coordinates and grid metadata" begin
        info = grid_meta((0.0, 0.0), (0.5, 0.25), (5, 4))
        data = Float64.(1:20)
        f = CartesianMeshField(data, info)
        new = get_interpolator(f); old = reference_interpolator(f)

        @test isnan(new(NaN, 0.5))
        @test isnan(new(0.5, NaN))
        @test isnan(new(NaN, NaN))

        # ±Inf clamps to the flat boundary; compare with the reference.
        for (x, y) in [(Inf, 0.5), (-Inf, 0.5), (0.5, Inf), (0.5, -Inf), (Inf, -Inf)]
            @test isfinite(new(x, y))
            @test new(x, y) ≈ old(x, y) atol = TOL
        end

        # Vector NaN propagates componentwise.
        vf = CartesianMeshField([VectorValue(1.0, 2.0) for _ in 1:20], info)
        vnew = get_interpolator(vf)
        v = vnew(NaN, 0.5)
        @test isnan(v[1]) && isnan(v[2])

        # Invalid spacing stays rejected, matching Interpolations.jl `scale`.
        for bad in (0.0, -1.0, NaN, Inf)
            @test_throws ArgumentError get_interpolator(
                CartesianMeshField(zeros(20), grid_meta((0.0, 0.0), (bad, 0.25), (5, 4))))
            @test_throws ArgumentError get_interpolator(
                CartesianMeshField(zeros(20), grid_meta((0.0, 0.0), (0.5, bad), (5, 4))))
        end

        # Data length must equal nx*ny (the old `reshape` threw the same way).
        @test_throws DimensionMismatch get_interpolator(
            CartesianMeshField(zeros(19), grid_meta((0.0, 0.0), (0.5, 0.25), (5, 4))))

        @test_throws OverflowError get_interpolator(CartesianMeshField(
            Float64[], grid_meta((0.0, 0.0), (1.0, 1.0), (typemax(Int), 2))))

        # Origins and upper bounds are validated before evaluation. Legacy
        # Interpolations.jl accepted these and returned NaN instead of throwing;
        # the deviation is deliberate (see ED9_wave1.md).
        for bad_origin in ((NaN, 0.0), (Inf, 0.0), (0.0, NaN), (0.0, -Inf))
            @test_throws ArgumentError get_interpolator(
                CartesianMeshField(zeros(20), grid_meta(bad_origin, (0.5, 0.25), (5, 4))))
        end
        @test_throws ArgumentError get_interpolator(CartesianMeshField(
            zeros(12), grid_meta((1.0e308, 0.0), (1.0e308, 0.25), (3, 4))))

        # Positive dimensions, then the legacy singleton rejection.
        for dims in ((0, 3), (3, 0), (0, 0))
            @test_throws ArgumentError get_interpolator(CartesianMeshField(
                zeros(dims[1] * dims[2]), grid_meta((0.0, 0.0), (0.5, 0.25), dims)))
        end
    end

    @testset "singleton grids rejected (legacy behaviour)" begin
        for dims in ((1, 3), (3, 1), (1, 1))
            f = CartesianMeshField(ones(dims[1] * dims[2]),
                                   grid_meta((2.0, -1.0), (0.5, 0.25), dims))
            @test_throws ArgumentError get_interpolator(f)
            @test_throws ArgumentError reference_interpolator(f)
        end
    end

    @testset "tensor-valued data (Gridap TensorValue)" begin
        nx, ny = 5, 7
        info = grid_meta((0.0, 0.0), (0.5, 0.25), (nx, ny))
        tdata = [TensorValue(randn(), randn(), randn(), randn()) for _ in 1:(nx * ny)]
        f = CartesianMeshField(tdata, info)
        new = get_interpolator(f); old = reference_interpolator(f)
        @test new(0.0, 0.0) isa TensorValue

        worst = 0.0
        for _ in 1:5_000
            x = info.origin[1] + (rand() * (nx + 1) - 0.5) * info.spacing[1]
            y = info.origin[2] + (rand() * (ny + 1) - 0.5) * info.spacing[2]
            a = new(x, y); b = old(x, y)
            for k in 1:4
                worst = max(worst, abs(a[k] - b[k]))
            end
        end
        @test worst < TOL

        tn = new(NaN, 0.5)
        @test isnan(tn[1]) && isnan(tn[2]) && isnan(tn[3]) && isnan(tn[4])
    end
end
