module TestCurvature

using Test
using EvolvingDomains
using Gridap: CartesianDiscreteModel, Measure, CellField, ∫
using Gridap.Geometry: get_node_coordinates
using GridapEmbedded: EmbeddedBoundary

# Curvature on the embedded interface.
#
# The headline check is the total-turning invariant: for any *smooth* simple closed
# curve ∮κ ds = 2π whatever its shape. It exercises the curvature, the surface measure
# and the sign convention in one assertion, and being shape-independent it must give
# the same number for the circle, the ellipse and the sign-changing flower.

function _geom(f, n)
    model = CartesianDiscreteModel((0.0, 1.0, 0.0, 1.0), (n - 1, n - 1))
    coords = vec(collect(get_node_coordinates(model)))
    EvolvingDiscreteGeometry([f((p[1], p[2])) for p in coords], model)
end

circle(p, c, r) = hypot(p[1] - c[1], p[2] - c[2]) - r
ellipse(p) = sqrt(((p[1] - 0.5) / 0.25)^2 + ((p[2] - 0.5) / 0.10)^2) - 1
function flower(p)
    x, y = p[1] - 0.5, p[2] - 0.5
    hypot(x, y) - 0.20 * (1 + 0.25cos(5atan(y, x)))
end

"∮κ ds over the embedded boundary"
function total_turning(geom)
    Γ = EmbeddedBoundary(ensure_cut!(geom))
    sum(∫(CellField(interface_curvature(geom), Γ))Measure(Γ, 2))
end

@testset "Curvature" begin

    @testset "circle: accuracy and sign" begin
        for (n, tol) in ((101, 0.06), (201, 0.03), (401, 0.02))
            geom = _geom(p -> circle(p, (0.5, 0.5), 0.15), n)
            κ = interface_curvature(geom)
            # κ = ∇·(∇φ/|∇φ|) is positive for a droplet with φ < 0 inside
            @test all(>(0), κ)
            @test maximum(abs, κ .- 1 / 0.15) / (1 / 0.15) < tol
        end
    end

    @testset "total turning is 2π" begin
        for (f, n, tol) in ((p -> circle(p, (0.5, 0.5), 0.15), 201, 0.01),
                            (ellipse, 201, 0.03),
                            (flower,  401, 0.05))
            @test isapprox(total_turning(_geom(f, n)), 2π; rtol=tol)
        end
    end

    @testset "kink: two drops one cell apart" begin
        # Differencing φ here gives κ spikes ~Δ/h that grow under refinement, because
        # the medial axis between the drops falls inside the stencil. Fitting on Γ
        # never touches it.
        for n in (201, 401)
            h = 1 / (n - 1); r = 0.15
            geom = _geom(p -> min(circle(p, (0.5, 0.5 - r - h / 2), r),
                                  circle(p, (0.5, 0.5 + r + h / 2), r)), n)
            @test maximum(abs, interface_curvature(geom)) < 2 / r
            @test isapprox(total_turning(geom), 4π; rtol=0.05)
        end
    end

    @testset "topology change: drops merging" begin
        # Nothing here reasons about topology -- `cut` handles that -- so the only
        # requirement is bounded curvature and a total turning matching the geometry
        # on both sides of the transition.
        #
        # The merged shape is NOT smooth: the union of two circles has two reflex
        # corners at the intersection points, so ∮κ ds is not 2π there. Each arc spans
        # 2(π−θ) with θ = acos(d/2r), giving ∮κ ds = 4(π−θ), which tends to 4π for a
        # shallow overlap since the corners then carry nearly all the deficit.
        n = 201; h = 1 / (n - 1); r = 0.15
        for gap_over_h in (4.0, 1.0, 0.25, -0.5, -2.0, -8.0)
            gap = gap_over_h * h
            geom = _geom(p -> min(circle(p, (0.5, 0.5 - r - gap / 2), r),
                                  circle(p, (0.5, 0.5 + r + gap / 2), r)), n)
            κ = interface_curvature(geom)
            @test all(isfinite, κ)
            @test maximum(abs, κ) < 20 / r

            if gap_over_h >= 1.0
                @test isapprox(total_turning(geom), 4π; rtol=0.05)     # two components
            elseif gap_over_h <= -0.5
                d = 2r + gap
                @test isapprox(total_turning(geom), 4 * (π - acos((d / 2) / r)); rtol=0.05)
            end
        end
    end

    @testset "curvature_at off the interface" begin
        geom = _geom(p -> circle(p, (0.5, 0.5), 0.15), 201)
        s = interface_samples(geom); h = 1 / 200
        @test isapprox(curvature_at(s, 0.5, 0.35 + h; radius=4h), 1 / 0.15; rtol=0.1)
        @test curvature_at(s, 0.5, 0.5; radius=4h) == 0.0      # far from Γ
    end

    @testset "cache invalidation" begin
        geom = _geom(p -> circle(p, (0.5, 0.5), 0.15), 101)
        κ1 = interface_curvature(geom)
        set_levelset!(geom, [circle((p[1], p[2]), (0.5, 0.5), 0.25)
                             for p in vec(collect(get_node_coordinates(geom.grid)))])
        κ2 = interface_curvature(geom)
        @test length(κ1) != length(κ2)
        @test isapprox(maximum(κ2), 1 / 0.25; rtol=0.1)
    end

end

end # module
