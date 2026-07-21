module TestHeleShawDispersion

using Test
using EvolvingDomains
using Gridap
using Gridap.Geometry: get_node_coordinates
using GridapEmbedded
using Printf

# One-phase Hele-Shaw with surface tension: linear dispersion relation.
#
# Reference: I. Lavi, N. Meunier, O. Pantz, "Implicit like time discretization for the
# one-phase Hele-Shaw problem with surface tension", ESAIM M2AN 59 (2025) 419-448,
# Section 8.1 and Appendix A.
#
# The model is
#
#   u + ∇p = 0,   ∇·u = 0   in Ω(t),      p = σκ   on Γ(t),      Vₙ = u·n
#
# discretized in the mixed form of the paper's eq. (4.1),
#
#   ∫_Ω (u·v - p ∇·v) dx + σ∫_Γ κ (n·v) ds = 0     ∀v
#   ∫_Ω (∇·u) q dx = 0                              ∀q
#
# The point of the mixed form is that p = σκ enters as a NATURAL condition through the
# surface term: there is no essential boundary condition, hence no Nitsche and no
# penalty parameter. Just as important, u is a primary unknown, so Vₙ = u·n is read
# straight off the solution instead of being post-processed from ∇p — and Vₙ is the
# only quantity this benchmark measures.
#
# (Solving the equivalent scalar problem Δp = 0 with p = σκ imposed by Nitsche instead
# was measurably worse: the m=2 rate did not converge at all, oscillating around 2%,
# versus 0.01% here. Getting accurate fluxes is exactly what the mixed form is for.)
#
# Perturbing a disc of radius R₀ as R(θ) = R₀ + ε cos(mθ), linear theory (Appendix A)
# predicts each Fourier mode decays as R_m(t) = R_m(0)exp(s_m t) with
#
#   s_m = -σ m (m² - 1) / R₀³
#
# For σ = 0.5, R₀ = 1 this is s₂ = -3, s₃ = -12, s₄ = -30, s₅ = -60; the paper measures
# -2.995, -11.996, -29.83, -59.50 with a body-fitted FreeFem++ discretization.
#
# We measure s_m from a SINGLE solve rather than a time loop. Explicit time stepping of
# this problem carries a Δt ≲ 2h³/(σπ³) surface-tension restriction, which at the
# resolution needed to even see ε = 0.05 would require O(10⁵-10⁶) steps to reach the
# decay timescale. The initial rate needs none of that: Vₙ projected onto cos(mθ) is
# dR_m/dt, and s_m = (dR_m/dt)/ε.

@testset "Hele-Shaw surface tension: dispersion relation" begin

    σ  = 0.5
    R₀ = 1.0
    ε  = 0.05
    n  = 200
    L  = 2.0
    grid = CartesianDiscreteModel((-L, L, -L, L), (n, n))
    h = 2L / n

    """Solve one perturbation mode; return (measured sₘ, area, ∮Vₙ ds, |u_cm|)."""
    function run_mode(m)
        φ = map(get_node_coordinates(grid)) do x
            hypot(x[1], x[2]) - (R₀ + ε * cos(m * atan(x[2], x[1])))
        end
        geom = EvolvingDiscreteGeometry(vec(collect(φ)), grid)
        cutgeo = ensure_cut!(geom)

        Ω  = Triangulation(cutgeo, PHYSICAL)
        Ωa = Triangulation(cutgeo, ACTIVE)
        Γ  = EmbeddedBoundary(cutgeo)
        nΓ = get_normal_vector(Γ)

        # Taylor-Hood P2/P1: inf-sup stable, and AgFEM handles the small cut cells
        aggs = aggregate(AggregateAllCutCells(), cutgeo)
        V = AgFEMSpace(TestFESpace(Ωa, ReferenceFE(lagrangian, VectorValue{2,Float64}, 2),
                                   conformity=:H1), aggs)
        Q = AgFEMSpace(TestFESpace(Ωa, ReferenceFE(lagrangian, Float64, 1),
                                   conformity=:H1), aggs)
        X = MultiFieldFESpace([TrialFESpace(V), TrialFESpace(Q)])
        Y = MultiFieldFESpace([V, Q])

        dΩ = Measure(Ω, 4)
        dΓ = Measure(Γ, 4)

        κ = CellField(interface_curvature(geom; radius=4h), Γ)
        a((u, p), (v, q)) = ∫(u ⋅ v - p * (∇ ⋅ v) + q * (∇ ⋅ u))dΩ
        l((v, q)) = ∫((-σ * κ) * (v ⋅ nΓ))dΓ

        uh, ph = solve(AffineFEOperator(a, l, X, Y))

        Vₙ = uh ⋅ nΓ                                         # primary unknown
        cosm = CellField(x -> cos(m * atan(x[2], x[1])), Γ)
        xcf  = CellField(x -> x[1], Γ)
        ycf  = CellField(x -> x[2], Γ)

        area   = sum(∫(1)dΩ)
        dRm    = sum(∫(Vₙ * cosm)dΓ) / π                     # ≈ dR_m/dt
        dArea  = sum(∫(Vₙ)dΓ)                                # mode 0: area rate
        ucm    = hypot(sum(∫(xcf * Vₙ)dΓ), sum(∫(ycf * Vₙ)dΓ)) / area   # mode 1
        return dRm / ε, area, dArea, ucm
    end

    exact(m) = -σ * m * (m^2 - 1) / R₀^3

    println("\nHele-Shaw dispersion relation  (σ=$σ, R₀=$R₀, ε=$ε, n=$n, h=$(round(h,digits=4)))")
    @printf("%4s %12s %12s %10s %14s %12s\n",
            "m", "s_m exact", "s_m meas", "err %", "area rate", "|u_cm|")

    for m in 2:5
        s, area, dArea, ucm = run_mode(m)
        e = exact(m)
        @printf("%4d %12.3f %12.3f %10.2f %14.2e %12.2e\n",
                m, e, s, 100 * (s - e) / abs(e), dArea, ucm)

        # Velocity scale of the mode. The area rate is a velocity times a length, so it
        # has to be compared against v_scale times the perimeter, not against v_scale.
        v_scale = abs(e) * ε

        # The tolerance is dominated by a finite-amplitude effect, not by the
        # discretization: we compare against LINEAR theory while identifying Vₙ with
        # dR/dt, and that identification is only accurate to O((εm)²). Measured errors
        # track (εm)²/2 in both size and m-dependence, and shrink with ε. At m=2, where
        # that term is smallest, the scheme converges to 0.01%.
        @test sign(s) == sign(e)                          # decaying, not growing
        @test isapprox(s, e; rtol=0.02 + (ε * m)^2)       # dispersion relation
        @test abs(dArea) < 0.05 * (2π * R₀) * v_scale     # area conserved (mode 0)
        @test ucm < 0.05 * v_scale                        # no self-propulsion (mode 1)
    end
end

end # module
