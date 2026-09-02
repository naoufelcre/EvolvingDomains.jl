module TestHeleShawST

#Reference: I. Lavi, N. Meunier, O. Pantz, "Implicit like time discretization for the
# one-phase Hele-Shaw problem with surface tension", ESAIM M2AN 59 (2025) 419-448,

# Relaxation of an arbitrary shape toward a circle under Hele-Shaw surface tension.
#
#   u + ∇p = 0,  ∇·u = 0  in Ω,   p = σκ on Γ,   Vₙ = u·n
#
# Run:  julia --project=. test/TestHeleShawST.jl
#
# The perimeter is a Lyapunov functional for this flow and must decrease.

using EvolvingDomains
using EvolvingDomains.Geometric, EvolvingDomains.Kinematic, EvolvingDomains.Transfer
using EvolvingDomains.Transfer: extend
using TransferOperator: prolong
using Gridap, GridapEmbedded, Printf, Statistics
include("Helpers/AggregationStrategy.jl")
using Gridap.Geometry: get_node_coordinates
using Gridap.TensorValues
using Gridap.CellData: get_array

const σ, R₀, L, n = 0.5, 1.0, 2.0, 200
const h  = 2L / n
const Δt = 2.5e-5
const NSTEP = 1000
const STRIDE = 10

shape(θ) = R₀ * (1 + 0.2cos(7θ) + 0.1sin(3θ))

grid = CartesianDiscreteModel((-L, L, -L, L), (n, n))
info = grid_info(grid)
XY   = vec(collect(get_node_coordinates(grid)))

const νc = 0.1     # tangential-filter strength (dimensionless; 0 disables)
                   # damping now comes from the package: tangential_smooth!(geom; strength, band)

# Removing a wiggle of amplitude a from r=R+a·cos(mθ) drops enclosed area by πa²/2 unless R
# grows to compensate — the incompressible physics does that growth automatically, the
# geometric filter does not. Restore it with one global normal shift: for an SDF, φ -= δ
# moves Γ out by δ, changing area by δ·P, so δ = (A − A_target)/P. Sussman–Fatemi (1999).
function current_area(geom)
    cg = ensure_cut!(geom)
    sum(∫(1)Measure(Triangulation(cg, PHYSICAL), 2))
end

function restore_area!(geom, A_target)
    cg = ensure_cut!(geom)
    A = sum(∫(1)Measure(Triangulation(cg, PHYSICAL), 2))
    P = sum(∫(1)Measure(EmbeddedBoundary(cg), 2))
    P > 1e-12 || return
    # φ -= δ moves Γ out by δ, raising area by δ·P. To shed excess area (A > A_target) we
    # move Γ in, i.e. φ += (A − A_target)/P. (The opposite sign runs the domain away.)
    geom.levelset .+= (A - A_target) / P
    invalidate!(geom.cache)
end

function solve_step(geom)
    ai = get_active_indices(geom, :current); cg = geom.cache.cut
    Ω, Ωa, Γ = Triangulation(cg, PHYSICAL), Triangulation(cg, ACTIVE), EmbeddedBoundary(cg)
    nΓ = get_normal_vector(Γ)
    aggs = aggregate(AggregateCutCellsByThreshold(1.0), cg)
    V = AgFEMSpace(TestFESpace(Ωa, ReferenceFE(lagrangian, VectorValue{2,Float64}, 2), conformity=:H1), aggs)
    Q = AgFEMSpace(TestFESpace(Ωa, ReferenceFE(lagrangian, Float64, 1), conformity=:H1), aggs)
    X = MultiFieldFESpace([TrialFESpace(V), TrialFESpace(Q)])
    Y = MultiFieldFESpace([V, Q])
    dΩ, dΓ = Measure(Ω, 4), Measure(Γ, 4)

    s   = interface_samples(geom)
    κv  = get_curvature(s; radius=8h)
    κ   = CellField(κv, Γ)
    a((u, p), (v, q)) = ∫(u ⋅ v - p * (∇ ⋅ v) + q * (∇ ⋅ u))dΩ
    l((v, q)) = ∫((-σ * κ) * (v ⋅ nΓ))dΓ
    uh, ph = solve(AffineFEOperator(a, l, X, Y))

    op   = GridMeshTransfer(info, V, ai, geom.levelset)
    vgr  = prolong(op, uh)
    vext = extend(geom, vgr)
    pgr  = prolong(GridMeshTransfer(info, Q, ai, geom.levelset), ph)
    (Ω=Ω, Γ=Γ, dΩ=dΩ, dΓ=dΓ, s=s, κ=κv, uh=uh, ph=ph, vext=vext, pgr=pgr,
     P=sum(∫(1)dΓ), A=sum(∫(1)dΩ))
end

function main()
    geom = EvolvingDiscreteGeometry([hypot(x[1], x[2]) - shape(atan(x[2], x[1])) for x in XY], grid)
    reinitialize!(geom)
    A0 = current_area(geom)        # conserved target: the initial enclosed area

    hist = (t=Float64[], P=Float64[], A=Float64[], vmax=Float64[], kmax=Float64[])
    t = 0.0

    plot(geom; label="Hele-Shaw step 0 / $NSTEP   t = 0.0")
    started = time_ns()
    for i in 1:NSTEP
        st = solve_step(geom)
        push!(hist.t, t); push!(hist.P, st.P); push!(hist.A, st.A)
        push!(hist.vmax, maximum(hypot(z[1], z[2]) for z in st.vext.data))
        push!(hist.kmax, maximum(abs, st.κ))
        mean_ms = (time_ns() - started) / (1e6 * i)
        if (i - 1) % STRIDE == 0
            plot(geom; label="Hele-Shaw step $i / $NSTEP   t = $(round(t, digits=4))   mean = $(round(mean_ms, digits=2)) ms/iteration")
        end
        advance!(geom, st.vext.data, Δt)
        filter_small_phase_islands!(geom.levelset, info; phase=:negative,
                                    min_component_size=8, connectivity=:edge)
        νc > 0 && tangential_smooth!(geom; strength=νc, band=3)
        reinitialize!(geom)
        restore_area!(geom, A0)     # undo the area a wiggle-removal necessarily sheds
        t += Δt
    end
    imin = argmin(hist.P)
    @printf("\nperimeter: start %.5f  min %.5f at t=%.4f (step %d)  end %.5f   [circle = %.5f]\n",
            hist.P[1], hist.P[imin], hist.t[imin], imin, hist.P[end], 2π * R₀)
    @printf("area drift: %.3e (%.3f %%)\n",
            maximum(abs, hist.A .- hist.A[1]), 100maximum(abs, hist.A .- hist.A[1]) / hist.A[1])
    return hist
end

end # module

TestHeleShawST.main()
