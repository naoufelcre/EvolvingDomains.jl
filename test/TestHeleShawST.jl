module TestHeleShawST

#Reference: I. Lavi, N. Meunier, O. Pantz, "Implicit like time discretization for the
# one-phase Hele-Shaw problem with surface tension", ESAIM M2AN 59 (2025) 419-448,

# Relaxation of an arbitrary shape toward a circle under Hele-Shaw surface tension,
# rendered as an animation with all fields.
#
#   u + ∇p = 0,  ∇·u = 0  in Ω,   p = σκ on Γ,   Vₙ = u·n
#
# Run:  julia --project=. test/TestHeleShawST.jl
#
# Every frame shows what the solver actually used: the level set and interface, the
# pressure on Ω, the extended velocity field, and κ on the subfacets that fed the
# surface-tension term. A diagnostics panel tracks the perimeter, which for this flow
# is a Lyapunov functional and must decrease — so the frame where it turns marks where
# discretisation error overtakes the physics.

using EvolvingDomains
using EvolvingDomains.Geometric, EvolvingDomains.Kinematic, EvolvingDomains.Transfer
using EvolvingDomains.Transfer: extend
using TransferOperator: prolong
using Gridap, GridapEmbedded, CairoMakie, Printf, Statistics
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
xs   = range(-L, L, length=n + 1)

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

function frame!(fig, st, geom, t, i, prng, hist)
    empty!(fig)
    φ  = geom.levelset
    φ2 = reshape(φ, n + 1, n + 1)
    lim = (-1.45, 1.45)
    px = [p[1] for p in st.s.positions]; py = [p[2] for p in st.s.positions]

    Label(fig[0, 1:4], @sprintf("Hele-Shaw relaxation   t = %.4f   step %d   P = %.5f",
                                t, i, st.P), fontsize=20)

    ax1 = Axis(fig[1, 1], aspect=DataAspect(), title="level set φ  +  Γ")
    heatmap!(ax1, xs, xs, φ2, colormap=:balance, colorrange=(-0.8, 0.8))
    contour!(ax1, xs, xs, φ2, levels=[0.0], color=:black, linewidth=2.5)
    limits!(ax1, lim..., lim...)

    ax2 = Axis(fig[1, 2], aspect=DataAspect(), title="pressure  p = σκ on Γ")
    pm = fill(NaN, length(φ))
    for k in eachindex(φ); φ[k] < 0 && (pm[k] = st.pgr.data[k]); end
    hm = heatmap!(ax2, xs, xs, reshape(pm, n + 1, n + 1), colormap=:blues, colorrange=prng)
    contour!(ax2, xs, xs, φ2, levels=[0.0], color=:black, linewidth=1.5)
    Colorbar(fig[2, 2], hm, vertical=false, height=10)
    limits!(ax2, lim..., lim...)

    ax3 = Axis(fig[1, 3], aspect=DataAspect(), title="velocity  u = −∇p")
    X_, Y_, U_, V_ = Float64[], Float64[], Float64[], Float64[]
    for j in 1:3:(n+1), ii in 1:3:(n+1)
        k = ii + (j - 1) * (n + 1)
        abs(φ[k]) <= 5h || continue
        push!(X_, xs[ii]); push!(Y_, xs[j])
        push!(U_, st.vext.data[k][1]); push!(V_, st.vext.data[k][2])
    end
    vmag = isempty(U_) ? 1.0 : maximum(hypot.(U_, V_))
    arrows!(ax3, X_, Y_, U_, V_, lengthscale=0.18 / max(vmag, 1e-12),
            arrowsize=5, linewidth=0.9, color=hypot.(U_, V_), colormap=:plasma)
    contour!(ax3, xs, xs, φ2, levels=[0.0], color=:black, linewidth=1.5)
    ax3.subtitle = @sprintf("max|u| = %.2f", vmag)
    limits!(ax3, lim..., lim...)

    ax4 = Axis(fig[1, 4], aspect=DataAspect(), title="curvature κ on Γ")
    sc = scatter!(ax4, px, py, color=st.κ, colormap=:turbo, markersize=6,
                  colorrange=(0.0, 3.0))
    Colorbar(fig[2, 4], sc, vertical=false, height=10)
    limits!(ax4, lim..., lim...)

    ax5 = Axis(fig[2, 1], title="perimeter", xlabel="t", height=110)
    lines!(ax5, hist.t, hist.P, color=:black)
    hlines!(ax5, [2π * R₀], color=:red, linestyle=:dash)
    ax6 = Axis(fig[2, 3], title="area", xlabel="t", height=110)
    lines!(ax6, hist.t, hist.A, color=:black)
    return nothing
end

function main()
    out = joinpath(@__DIR__, "output_relax")
    # clear stale frames only — wiping the directory also destroys the assembled gif
    mkpath(out)
    foreach(f -> rm(joinpath(out, f)), filter(f -> startswith(f, "f") && endswith(f, ".png"), readdir(out)))

    geom = EvolvingDiscreteGeometry([hypot(x[1], x[2]) - shape(atan(x[2], x[1])) for x in XY], grid)
    reinitialize!(geom)
    A0 = current_area(geom)        # conserved target: the initial enclosed area

    hist = (t=Float64[], P=Float64[], A=Float64[], vmax=Float64[], kmax=Float64[])
    fig = Figure(size=(1800, 780))
    t = 0.0; nf = 0
    st = solve_step(geom)
    prng = (minimum(σ .* st.κ), maximum(σ .* st.κ))

    @printf("shape → circle.  Δt=%.0e  %d steps  frames every %d\n", Δt, NSTEP, STRIDE)
    for i in 1:NSTEP
        st = solve_step(geom)
        push!(hist.t, t); push!(hist.P, st.P); push!(hist.A, st.A)
        push!(hist.vmax, maximum(hypot(z[1], z[2]) for z in st.vext.data))
        push!(hist.kmax, maximum(abs, st.κ))
        if (i - 1) % STRIDE == 0
            frame!(fig, st, geom, t, i, prng, hist)
            save(joinpath(out, @sprintf("f%04d.png", nf)), fig); nf += 1
        end
        i % 50 == 0 && @printf("  step %4d  t=%.4f  P=%.5f  A=%.5f  max|u|=%7.2f  max|κ|=%.2f\n",
                               i, t, st.P, st.A, hist.vmax[end], hist.kmax[end])
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
    println("frames: $nf in $out")
    return out, nf
end

end # module

TestHeleShawST.main()
