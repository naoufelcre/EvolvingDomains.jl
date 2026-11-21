#!/usr/bin/env julia
# PlanarStefan -- one-phase Stefan melting in 2-D: EvolvingDomains level-set geometry coupled
# to an independent cut finite element thermal solver (Gridap + GridapEmbedded).
#
#   liquid       Omega(t) = {(x,y) : 0 < x < s(t)},  level set  phi = x - s(t),  liquid = phi < 0
#   temperature  rho cp T_t = k Lap T in Omega(t),  T = Th at x = 0,  T = Tm at x = s(t),
#                dT/dy = 0 at y = 0,1
#   interface    rho L ds/dt = -k dT/dx(s-, t)
#
# Run (from the package root):
#   julia --project=validation/PlanarStefan validation/PlanarStefan/run.PlanarStefan.jl
#   julia --project=validation/PlanarStefan validation/PlanarStefan/run.PlanarStefan.jl \
#         --n 129 --dt 0.002 --gd 100 --out validation/PlanarStefan/output/fine --no-viz
#
# The native CutCells library used by EvolvingDomains `cut` is taken from
# $EVOLVINGDOMAINS_CUTCELLS_LIBRARY, or from the locally built copy in ./build
# (see scripts/build_cxx_bindings.sh in the package root).
#
# Component separation (no opaque framework type):
#   geometry       EvolvingDomains Grid/Field (phi and T samples) + `cut`; phi is advanced by
#                  the package's WENO5 + SSP-RK3 `eulerianTransport!`
#   discretization GridapEmbedded cut FEM on a background Cartesian mesh; the level-set array
#                  handed to GridapEmbedded *is* the EvolvingDomains Field
#   time stepping  implicit Euler for T, first-order explicit splitting for the interface
#   diagnostics    every reported number is computed from the numerical state

using Gridap, GridapEmbedded, StaticArrays, LinearAlgebra, SparseArrays, Printf, SpecialFunctions
using GridapEmbedded.LevelSetCutters: DiscreteGeometry
import Gridap.Arrays: collect1d
import EvolvingDomains as ED

# =============================================================================
# Parameters
# =============================================================================

function defaultParams(; n=65, dt=0.002, t0=0.05, tf=0.20, gd=100.0, gg=0.1, out="output")
    return (ρ=1.0, cp=1.0, k=0.5, L=1.0, Th=1.0, Tm=0.0,
        n=n, order=1, dt=dt, t0=t0, tf=tf, γd=gd, γg=gg, out=out)
end

alpha(p) = p.k / (p.ρ * p.cp)
stefanNumber(p) = p.cp * (p.Th - p.Tm) / p.L

# =============================================================================
# Analytical similarity solution
# =============================================================================

"""λ > 0 solving Ste = sqrt(pi) λ exp(λ^2) erf(λ), by bisection on [0, 2]."""
function similarityRoot(Ste; lo=0.0, hi=2.0, tol=1e-15)
    f(lam) = sqrt(π) * lam * exp(lam^2) * erf(lam) - Ste
    flo = f(lo)
    for _ in 1:500
        mid = (lo + hi) / 2
        (f(mid) < 0) == (flo < 0) ? (lo = mid) : (hi = mid)
        hi - lo < tol * max(1.0, hi) && break
    end
    return (lo + hi) / 2
end

sExact(λ, α, t) = 2 * λ * sqrt(α * t)

"""T_exact(x,t) = Th - (Th-Tm) erf(x/(2 sqrt(α t))) / erf(λ)."""
TExact(λ, α, p, x, t) = p.Th - (p.Th - p.Tm) * erf(x / (2 * sqrt(α * t))) / erf(λ)

"""∂x T_exact(x,t): the liquid-side derivative (negative while melting)."""
dTdxExact(λ, α, p, x, t) =
    -(p.Th - p.Tm) / erf(λ) * exp(-(x / (2 * sqrt(α * t)))^2) / (sqrt(π) * sqrt(α * t))

# =============================================================================
# EvolvingDomains state: background grid, level set, temperature samples
# =============================================================================

function configureCutCellsLibrary!()
    path = get(ENV, "EVOLVINGDOMAINS_CUTCELLS_LIBRARY", "")
    isempty(path) || return path
    ext = Sys.iswindows() ? "dll" : Sys.isapple() ? "dylib" : "so"
    local_lib = joinpath(@__DIR__, "build", "libevolvingdomains_cut_cells." * ext)
    isfile(local_lib) || error("no CutCells library: set EVOLVINGDOMAINS_CUTCELLS_LIBRARY, " *
                               "or build one with scripts/build_cxx_bindings.sh and copy it into " *
                               joinpath(@__DIR__, "build"))
    ENV["EVOLVINGDOMAINS_CUTCELLS_LIBRARY"] = local_lib
    return local_lib
end

"""Grid, level-set Field and temperature Field on one uniform Cartesian mesh."""
function initialState(p, λ)
    α = alpha(p)
    x = collect(range(0.0, 1.0; length=p.n))
    s0 = sExact(λ, α, p.t0)
    φ = [xi - s0 for xi in x, _ in x]                            # phi < 0: liquid
    T = [xi < s0 ? TExact(λ, α, p, xi, p.t0) : p.Tm for xi in x, _ in x]
    grid = ED.Grid((x, x))
    disc = ED.Discretization(Dict(:init => grid))
    return (x=x, h=x[2] - x[1], grid=grid, disc=disc,
        level=ED.Field(disc, Dict(:init => φ)),
        T=ED.Field(disc, Dict(:init => T)))
end

# =============================================================================
# Geometry: EvolvingDomains Field -> GridapEmbedded level set -> cut discretization
# =============================================================================

"""Gridap node id -> (i,j) of the EvolvingDomains grid, checked to be exact."""
function nodeIndexMap(state, model)
    h = state.h
    X = collect1d(Gridap.ReferenceFEs.get_node_coordinates(model))
    map = [(round(Int, pt[1] / h) + 1, round(Int, pt[2] / h) + 1) for pt in X]
    dev = maximum(max(abs(X[q][1] - state.x[map[q][1]]), abs(X[q][2] - state.x[map[q][2]]))
                  for q in eachindex(X))
    dev == 0.0 || error("Gridap nodes and the EvolvingDomains grid disagree by $dev")
    return X, map
end

"""Facet mask selecting x = 0 (hot wall); y = 0,1 stay natural (insulated)."""
function wallFacetMask(model, X)
    topo = Gridap.Geometry.get_grid_topology(model)
    D = Gridap.Geometry.num_cell_dims(model)
    facets = Gridap.Geometry.get_faces(topo, D - 1, 0)
    return [all(abs(X[q][1]) < 1e-12 for q in facets[f]) for f in eachindex(facets)]
end

"""
The adapter: hand the EvolvingDomains level-set samples to GridapEmbedded as a discrete
geometry on the background model, then cut that model. `cutgeo` (GridapEmbedded) and
`ED.cut(level)` (native cut cells) are then two independent cuttings of the same array.
"""
function buildGeometry(state, model, X, nodeindex)
    φnode = [state.level.data[:init][nodeindex[q]...] for q in eachindex(X)]
    geo = DiscreteGeometry(φnode, X)
    cutgeo = cut(model, geo)
    return (geo=geo, cutgeo=cutgeo,
        bgstate=compute_bgcell_to_inoutcut(cutgeo, geo),
        Ω=Triangulation(cutgeo, PHYSICAL),
        Γ=EmbeddedBoundary(cutgeo),
        S=GhostSkeleton(cutgeo))
end

"""Background Q1 space, wall facet triangulation and dof->node map (all geometry independent)."""
function buildDiscretization(state, p, model, nodeindex)
    X = collect1d(Gridap.ReferenceFEs.get_node_coordinates(model))
    reffe = ReferenceFE(lagrangian, Float64, p.order)
    V = TestFESpace(Triangulation(model), reffe; conformity=:H1)
    U = TrialFESpace(V)
    deg = 2 * p.order + 2
    Λ = BoundaryTriangulation(model, wallFacetMask(model, X))
    dofnode = dofNodeMap(V, model, state, nodeindex)
    return (V=V, U=U, Λ=Λ, dΛ=Measure(Λ, deg), nΛ=get_normal_vector(Λ), deg=deg,
        dofnode=dofnode)
end

"""Gridap dof -> EvolvingDomains grid node, verified against the dof coordinates."""
function dofNodeMap(V, model, state, nodeindex)
    cd = get_cell_dof_ids(V)
    cn = Gridap.Geometry.get_cell_node_ids(model)
    dofnode = zeros(Int, num_free_dofs(V))
    for c in 1:Gridap.Geometry.num_cells(model), (d, nn) in zip(cd[c], cn[c])
        dofnode[d] = nn
    end
    all(!iszero, dofnode) || error("some dof was not matched to a background node")
    Xg = collect1d(Gridap.FESpaces.get_free_dof_coordinates(V))
    dev = maximum(max(abs(Xg[d][1] - state.x[nodeindex[dofnode[d]][1]]),
        abs(Xg[d][2] - state.x[nodeindex[dofnode[d]][2]])) for d in eachindex(Xg))
    dev == 0.0 || error("dof/node mismatch of $dev")
    return dofnode
end

# =============================================================================
# Thermal step: assemble M and K, solve (M + dt K) u = M u_old + dt r
# =============================================================================

"""
Assemble the discrete thermal operators of rho cp T_t = k Lap T on the current liquid domain:

    M  = ∫ rho cp v u dΩ                          heat-capacity mass matrix
    K  = ∫ k ∇v·∇u dΩ                            conductivity stiffness matrix
         + Nitsche terms of the interface (T = Tm) and of the wall (T = Th)
         + ghost penalty on the facets adjacent to cut cells
    r  = wall boundary contribution

Both Dirichlet conditions are imposed weakly (symmetric Nitsche), so no equation is eliminated
algebraically and M and K stay inspectable.
"""
function assembleThermalOperators(p, geom, dis, dt)
    k, m = p.k, p.ρ * p.cp
    h = 1.0 / (p.n - 1)
    dΩ = Measure(geom.Ω, dis.deg)
    dΓ = Measure(geom.Γ, dis.deg)
    dS = Measure(geom.S, dis.deg)
    nΓ = get_normal_vector(geom.Γ)
    nS = get_normal_vector(geom.S)
    dΛ, nΛ = dis.dΛ, dis.nΛ

    mass(u, v) = ∫(m * v * u) * dΩ
    stiff(u, v) = ∫(k * (∇(v) ⋅ ∇(u))) * dΩ +
                  ∫(-k * (nΓ ⋅ ∇(u)) * v - k * (nΓ ⋅ ∇(v)) * u + (p.γd * k / h) * u * v) * dΓ +
                  ∫(-k * (nΛ ⋅ ∇(u)) * v - k * (nΛ ⋅ ∇(v)) * u + (p.γd * k / h) * u * v) * dΛ +
                  ∫(p.γg * k * h * jump(nS ⋅ ∇(u)) * jump(nS ⋅ ∇(v))) * dS
    wall(v) = ∫((p.γd * k / h) * p.Th * v - k * (nΛ ⋅ ∇(v)) * p.Th) * dΛ

    M = assemble_matrix(mass, dis.U, dis.V)
    K = assemble_matrix(stiff, dis.U, dis.V)
    r = assemble_vector(wall, dis.V)
    return (M=M, K=K, r=r, A=M + dt * K, dΩ=dΩ, dΓ=dΓ,
        cellvolume=get_cell_measure(geom.Ω),
        masssum=sum(M), symdevM=symmetryDeviation(M), symdevK=symmetryDeviation(K))
end

"""max |A - Aᵀ| / max |A|: assembly-order roundoff, not a structural asymmetry."""
symmetryDeviation(A) = maximum(abs.(A - transpose(A))) / maximum(abs.(A))

"""
One implicit Euler thermal step on Ω(t), written back into the EvolvingDomains Field.

Transfer of the previous state: the mesh does not move and the dofs are the grid nodes, so the
previous nodal state is reused unchanged (no interpolation).  A node whose cell holds no liquid
has an empty equation; it is held at the solid temperature Tm (= the interface value) instead of
being padded with zeros, and a node that the liquid has just gained therefore starts from the
state it had in the solid, i.e. Tm.
"""
function solveTemperature(state, p, geom, dis, dofnode, nodeindex; dt)
    ops = assembleThermalOperators(p, geom, dis, dt)
    Tnodal = state.T.data[:init]
    uold = [Tnodal[nodeindex[dofnode[d]]...] for d in eachindex(dofnode)]
    b = ops.M * uold + dt * ops.r
    A = ops.A

    diagA = abs.(diag(A))
    active = findall(>(1e-12 * maximum(diagA)), diagA)
    inactmask = trues(length(diagA))
    inactmask[active] .= false
    breg = copy(b)
    breg[inactmask] .= p.Tm
    Areg = A + spdiagm(0 => Float64.(inactmask))
    u = Areg \ breg

    unew = fill(p.Tm, p.n, p.n)
    for d in eachindex(u)
        unew[nodeindex[dofnode[d]]...] = u[d]
    end
    state.T.data[:init] .= unew
    return (u=u, uold=uold, ops=ops, Areg=Areg, active=active, A=A, b=breg,
        ndof=length(u), nactive=length(active),
        residual=norm(Areg * u - breg) / max(norm(breg), eps()),
        mindiag=minimum(diagA[active]), maxdiag=maximum(diagA),
        mincut=minimum(ops.cellvolume) / (1.0 / (p.n - 1))^2)
end

# =============================================================================
# Interface heat flux, Stefan velocity, extension, geometry advance
# =============================================================================

"""Liquid-side normal heat flux at the interface: k <∇T·n>_Γ, n the outward liquid normal."""
function interfaceHeatFlux(therm, geom, p, dis)
    uh = FEFunction(dis.U, therm.u)
    nΓ = get_normal_vector(geom.Γ)
    area = sum(∫(1.0) * therm.ops.dΓ)
    q = sum(∫(p.k * (∇(uh) ⋅ nΓ)) * therm.ops.dΓ) / area
    nmean = sum(∫(nΓ) * therm.ops.dΓ) / area
    return q, area, nmean
end

"""One-phase Stefan condition rho L ds/dt = -k dT/dx(s-): interface speed from the flux."""
dynamics(q, p) = -q / (p.ρ * p.L)

"""
Extend the interface speed into a Cartesian nodal velocity field (constant extension).
The interface is the straight line x = s(t) with normal (1,0), so filling (V, 0) is the exact
extension of the normal speed; no PDE extension problem is solved.
"""
velocityFromRates(V, state) = vec([SVector(V, 0.0) for _ in state.x, _ in state.x])

"""Advance the level set with the package's WENO5 + SSP-RK3 Eulerian transport."""
function advanceGeometry!(state, velocity, dt)
    ED.eulerianTransport!(state.level, velocity, dt; level_id=:init, cache=ED.WENO5Cache())
    return nothing
end

# =============================================================================
# Diagnostics computed from the numerical state
# =============================================================================

"""Interface position from the level-set samples: zero crossing interpolated at y = 0.5."""
function interfacePosition(state)
    φ = state.level.data[:init]
    n = size(φ, 1)
    j = (n + 1) ÷ 2
    i = findfirst(>(0), view(φ, :, j))
    (i === nothing || i == 1) && return NaN
    a, b = φ[i - 1, j], φ[i, j]
    return state.x[i - 1] + a / (a - b) * state.h
end

"""Temperature error and planarity metrics on the nodes currently inside the liquid."""
function temperatureErrors(state, λ, p, t)
    α = alpha(p)
    φ = state.level.data[:init]
    T = state.T.data[:init]
    n = size(φ, 1)
    peak, sq, cnt = 0.0, 0.0, 0
    for j in 1:n, i in 1:n
        φ[i, j] < 0 || continue
        e = T[i, j] - TExact(λ, α, p, state.x[i], t)
        peak = max(peak, abs(e))
        sq += e^2
        cnt += 1
    end
    yspread = 0.0
    for i in 1:n
        vals = [T[i, j] for j in 1:n if φ[i, j] < 0]
        length(vals) > 1 && (yspread = max(yspread, maximum(vals) - minimum(vals)))
    end
    return (max=peak, rms=cnt == 0 ? 0.0 : sqrt(sq / cnt), n=cnt, yspread=yspread)
end

"""
Cross-check of the two independent cuttings: compare the per-cell state of `ED.cut` (native cut
cells) with `compute_bgcell_to_inoutcut` (GridapEmbedded), matched by cell centre so the two
cell orderings do not matter.
"""
function compareCuttings(state, model, geom, snapshot)
    n = size(state.level.data[:init], 1)
    h = state.h
    cd = Gridap.Geometry.get_cell_node_ids(model)
    X = collect1d(Gridap.ReferenceFEs.get_node_coordinates(model))
    clin = LinearIndices((n - 1, n - 1))
    mismatch, nin, ncut, nout = 0, 0, 0, 0
    for c in 1:Gridap.Geometry.num_cells(model)
        xc = 0.0
        yc = 0.0
        for node in cd[c]
            xc += X[node][1]
            yc += X[node][2]
        end
        nc = length(cd[c])
        i = clamp(floor(Int, xc / nc / h) + 1, 1, n - 1)
        j = clamp(floor(Int, yc / nc / h) + 1, 1, n - 1)
        edstate = snapshot.bg_states[clin[i, j]]
        gdstate = geom.bgstate[c]
        mismatch += (Int8(edstate) != Int8(gdstate))
        gdstate == Int8(-1) && (nin += 1)
        gdstate == Int8(0) && (ncut += 1)
        gdstate == Int8(1) && (nout += 1)
    end
    return (mismatch=mismatch, nin=nin, ncut=ncut, nout=nout)
end

"""Spectral condition number of the reduced active system (small problems only)."""
function conditioningEstimate(therm)
    na = length(therm.active)
    na > 1200 && return NaN                      # dense SVD would be the run's whole budget
    Ablk = Matrix(therm.Areg[therm.active, therm.active])
    return cond(Ablk)
end

"""Dump the mass and stiffness triplets of the final step for independent inspection."""
function dumpOperators(therm, outdir)
    open(joinpath(outdir, "operators.csv"), "w") do io
        println(io, "matrix,dof_i,dof_j,value")
        for (name, A) in (("M", therm.ops.M), ("K", therm.ops.K)), j in 1:size(A, 2)
            for k in nzrange(A, j)
                row = A.rowval[k]
                row <= j || continue
                @printf(io, "%s,%d,%d,%.16e\n", name, row, j, A.nzval[k])
            end
        end
    end
    return nothing
end

# =============================================================================
# Simulation
# =============================================================================

function simulate(p)
    λ = similarityRoot(stefanNumber(p))
    @assert abs(sqrt(π) * λ * exp(λ^2) * erf(λ) - stefanNumber(p)) < 1e-12 "similarity root not converged"
    α = alpha(p)
    mkpath(p.out)
    configureCutCellsLibrary!()
    state = initialState(p, λ)
    model = CartesianDiscreteModel((0.0, 1.0, 0.0, 1.0), (p.n - 1, p.n - 1))
    X, nodeindex = nodeIndexMap(state, model)
    dis = buildDiscretization(state, p, model, nodeindex)
    nsteps = round(Int, (p.tf - p.t0) / p.dt)

    @printf("PlanarStefan: n=%d (h=%.6f) order=%d dt=%.5f steps=%d γd=%.1f γg=%.2f\n",
        p.n, state.h, p.order, p.dt, nsteps, p.γd, p.γg)
    @printf("Ste=%.4f α=%.4f λ=%.12f s(t0)=%.6f s(tf)=%.6f\n", stefanNumber(p), α, λ,
        sExact(λ, α, p.t0), sExact(λ, α, p.tf))

    hist = NamedTuple[]
    profiles = NamedTuple[]
    snaps = NamedTuple[]
    snap_every = max(1, div(nsteps, 7))
    prev_active = Set{Int}()
    worst = (cut=Inf, cond=NaN, condt=NaN)
    tstart = time()
    recordProfile!(t) = push!(profiles, (t=t, x=copy(state.x),
        T=copy(state.T.data[:init][:, (p.n + 1) ÷ 2]),
        Tex=[TExact(λ, α, p, xi, t) for xi in state.x],
        φ=copy(state.level.data[:init][:, (p.n + 1) ÷ 2])))
    recordSnapshot!(t) = push!(snaps, (t=t, x=copy(state.x),
        T=copy(state.T.data[:init]), φ=copy(state.level.data[:init])))

    recordSnapshot!(p.t0)
    recordProfile!(p.t0)

    t = p.t0
    prev = (area=NaN, line=NaN, q=NaN, V=NaN)
    for step in 0:(nsteps - 1)
        # 1. current geometry and cut
        geom = buildGeometry(state, model, X, nodeindex)
        snapshot = ED.cut(state.level, :init)
        cmp = compareCuttings(state, model, geom, snapshot)
        volume = get_cell_measure(geom.Ω)
        s_area = sum(volume)                 # liquid area = interface position (height 1)
        s_line = interfacePosition(state)    # midline zero crossing of the level set

        # 2. temperature on Ω(t); 3. interface flux; 4. Stefan speed
        therm = solveTemperature(state, p, geom, dis, dis.dofnode, nodeindex; dt=p.dt)
        q, ifacearea, nmean = interfaceHeatFlux(therm, geom, p, dis)
        Vstefan = dynamics(q, p)
        tT = t + p.dt

        # diagnostics of the thermal solve
        err = temperatureErrors(state, λ, p, tT)
        newactive = step == 0 ? Int[] : [d for d in therm.active if !(d in prev_active)]
        initdev = isempty(newactive) ? 0.0 : maximum(abs(therm.uold[d] - p.Tm) for d in newactive)
        prev_active = Set(therm.active)
        ds_area = isnan(prev.area) ? NaN : s_area - prev.area
        ds_line = isnan(prev.line) ? NaN : s_line - prev.line
        rS = isnan(prev.q) ? NaN : p.ρ * p.L * (ds_area / p.dt) + prev.q   # rho L V_n + k dT/dx
        s_ex_num = sExact(λ, α, t)
        q_ex_num = p.k * dTdxExact(λ, α, p, s_area, tT)
        condest = (step == 0 || step == nsteps - 1) ? conditioningEstimate(therm) : NaN
        worst = (cut=min(worst.cut, therm.mincut),
            cond=isnan(condest) ? worst.cond : condest,
            condt=isnan(condest) ? worst.condt : t)
        step == nsteps - 1 && dumpOperators(therm, p.out)

        massdev = therm.ops.masssum / (p.ρ * p.cp * s_area) - 1.0
        @assert abs(massdev) < 1e-8 "mass matrix violates the partition of unity: $massdev"
        @assert max(therm.ops.symdevM, therm.ops.symdevK) < 1e-12 "assembled operator is not symmetric"
        push!(hist, (step=step, t=t, t_T=tT, s_num=s_area, s_line=s_line,
            ls_dev=s_line - s_area, s_exact=s_ex_num, s_err=s_area - s_ex_num,
            q_num=q, q_exact=q_ex_num, q_rel=(q - q_ex_num) / abs(q_ex_num),
            V_stefan=Vstefan, ds_area=ds_area, ds_line=ds_line,
            ds_flux=isnan(prev.V) ? NaN : prev.V * p.dt, r_stefan=rS,
            T_err_max=err.max, T_err_rms=err.rms, T_yspread=err.yspread,
            ndof=therm.ndof, nactive=therm.nactive, nnew=length(newactive), initdev=initdev,
            residual=therm.residual, dimratio=therm.mindiag / therm.maxdiag,
            mincut=therm.mincut, ifacearea=ifacearea, normal_x=nmean[1],
            cutmismatch=cmp.mismatch, cutcells=cmp.ncut,
            massdev=massdev,
            symdevM=therm.ops.symdevM, symdevK=therm.ops.symdevK))
        (step % snap_every == 0 || step == nsteps - 1) && recordSnapshot!(tT)
        (step % snap_every == 0 || step == nsteps - 1) && recordProfile!(tT)

        # 5. extend the speed;  6. advance the level set;  7. the state moves with the dofs
        velocity = velocityFromRates(Vstefan, state)
        advanceGeometry!(state, velocity, p.dt)
        prev = (area=s_area, line=s_line, q=q, V=Vstefan)
        t = tT
    end
    elapsed = time() - tstart

    # geometry at t_final, for the final interface comparison
    finalgeom = buildGeometry(state, model, X, nodeindex)
    s_final = sum(get_cell_measure(finalgeom.Ω))
    s_final_line = interfacePosition(state)
    s_final_exact = sExact(λ, α, p.tf)
    finalerr = temperatureErrors(state, λ, p, p.tf)

    writeOutputs(p, hist, profiles, snaps, (λ=λ, α=α, nsteps=nsteps, elapsed=elapsed,
        s_final=s_final, s_final_line=s_final_line, s_final_exact=s_final_exact,
        finalerr=finalerr, worst=worst))

    return (p=p, λ=λ, α=α, hist=hist, profiles=profiles, snaps=snaps, elapsed=elapsed,
        s_final=s_final, s_final_line=s_final_line, s_final_exact=s_final_exact,
        finalerr=finalerr, worst=worst)
end

# =============================================================================
# Output
# =============================================================================

function writeOutputs(p, hist, profiles, snaps, summary)
    open(joinpath(p.out, "history.csv"), "w") do io
        println(io, join(string.(keys(hist[1])), ","))
        for h in hist
            println(io, join((@sprintf("%.12e", getfield(h, k)) for k in keys(h)), ","))
        end
    end
    open(joinpath(p.out, "profiles.csv"), "w") do io
        println(io, "t,x,T_num,T_exact,phi")
        for (r, pr) in enumerate(profiles), i in eachindex(pr.x)
            @printf(io, "%.12e,%.12e,%.12e,%.12e,%.12e\n", pr.t, pr.x[i], pr.T[i], pr.Tex[i], pr.φ[i])
        end
    end
    open(joinpath(p.out, "snapshots.csv"), "w") do io
        println(io, "t,i,j,x,y,T,phi")
        for sn in snaps, j in eachindex(sn.x), i in eachindex(sn.x)
            @printf(io, "%.12e,%d,%d,%.12e,%.12e,%.12e,%.12e\n", sn.t, i, j, sn.x[i], sn.x[j],
                sn.T[i, j], sn.φ[i, j])
        end
    end

    h = hist
    txt = IOBuffer()
    println(txt, "PlanarStefan one-phase Stefan melting, EvolvingDomains geometry + GridapEmbedded cut FEM")
    @printf(txt, "grid n=%d (h=%.6f)  Q1 cut FEM  order=%d  dt=%.5f  steps=%d\n",
        p.n, 1.0 / (p.n - 1), p.order, p.dt, summary.nsteps)
    @printf(txt, "rho=%.3f cp=%.3f k=%.3f L=%.3f Th=%.3f Tm=%.3f  Ste=%.4f  alpha=%.4f\n",
        p.ρ, p.cp, p.k, p.L, p.Th, p.Tm, stefanNumber(p), summary.α)
    @printf(txt, "similarity root lambda=%.12f (residual %.2e)\n", summary.λ,
        sqrt(π) * summary.λ * exp(summary.λ^2) * erf(summary.λ) - stefanNumber(p))
    @printf(txt, "Nitsche gamma_d=%.1f, ghost gamma_g=%.2f, implicit Euler, explicit interface splitting\n",
        p.γd, p.γg)
    @printf(txt, "runtime %.2f s (%.1f ms/step)\n", summary.elapsed, 1e3 * summary.elapsed / summary.nsteps)
    @printf(txt, "\ninterface (liquid area of the cut domain):  s_num(tf)=%.8f  s_exact(tf)=%.8f  error=%.3e (%.3f%%)\n",
        summary.s_final, summary.s_final_exact, summary.s_final - summary.s_final_exact,
        100 * (summary.s_final - summary.s_final_exact) / summary.s_final_exact)
    @printf(txt, "  midline zero crossing of the level set at tf = %.8f (deviation from area value %.3e)\n",
        summary.s_final_line, summary.s_final_line - summary.s_final)
    @printf(txt, "  max |s error| over run = %.3e;  max level-set distortion (line vs area) = %.3e\n",
        maximum(abs.(getfield.(h, :s_err))), maximum(abs.(getfield.(h, :ls_dev))))
    @printf(txt, "temperature: max nodal error = %.3e (over run), rms at final state = %.3e\n",
        maximum(getfield.(h, :T_err_max)), summary.finalerr.rms)
    @printf(txt, "planarity: max y-variation = %.3e (exact problem is 1-D)\n",
        maximum(getfield.(h, :T_yspread)))
    @printf(txt, "heat flux: |q_num| rel. error mean %.3f%%, max %.3f%% (first step); range [%.4f, %.4f]\n",
        100 * sum(abs.(getfield.(h, :q_rel))) / length(h), 100 * maximum(abs.(getfield.(h, :q_rel))),
        minimum(getfield.(h, :q_num)), maximum(getfield.(h, :q_num)))
    rsv = [r for r in getfield.(h, :r_stefan) if !isnan(r)]
    @printf(txt, "Stefan residual rho L V_n + k dT/dx (V_n = d(liquid area)/dt): mean %.3e, max %.3e\n",
        sum(abs.(rsv)) / length(rsv), maximum(abs.(rsv)))
    dsm = [a - b for (a, b) in zip(getfield.(h, :ds_area), getfield.(h, :ds_flux)) if !isnan(a) && !isnan(b)]
    @printf(txt, "measured vs flux-predicted motion: max |ds_area - ds_flux| = %.3e over the run\n",
        maximum(abs.(dsm)))
    @printf(txt, "linear systems: max residual %.3e, active dofs %d..%d of %d, min diag ratio %.3e\n",
        maximum(getfield.(h, :residual)), minimum(getfield.(h, :nactive)),
        maximum(getfield.(h, :nactive)), last(h).ndof, minimum(getfield.(h, :dimratio)))
    @printf(txt, "cut cells: %d per step, smallest cut-cell volume fraction %.3e, new active dofs max %d (init |T-Tm| max %.1e)\n",
        last(h).cutcells, summary.worst.cut, maximum(getfield.(h, :nnew)),
        maximum(getfield.(h, :initdev)))
    @printf(txt, "spectral cond of the active system: %.3e (t=%.3f), interface area %.6f\n",
        summary.worst.cond, summary.worst.condt, last(h).ifacearea)
    @printf(txt, "discrete operators: |sum(M)/(rho cp liquid area) - 1| <= %.3e (partition of unity),\n",
        maximum(abs.(getfield.(h, :massdev))))
    @printf(txt, "  symmetry |A-A'|_max/|A|_max: M %.2e, K %.2e (roundoff); triplets of the last step in operators.csv\n",
        maximum(getfield.(h, :symdevM)), maximum(getfield.(h, :symdevK)))
    @printf(txt, "cutting cross-check: max per-cell state mismatch vs ED.cut = %d of %d cells\n",
        maximum(getfield.(h, :cutmismatch)), (p.n - 1)^2)
    print(txt, "\nfiles: history.csv, profiles.csv, snapshots.csv, summary.txt\n")
    write(joinpath(p.out, "summary.txt"), String(take!(txt)))
    print(read(joinpath(p.out, "summary.txt"), String))
    return nothing
end

# =============================================================================
# Entry point
# =============================================================================

function parseArgs(args)
    opts = Dict{String,Any}()
    i = 1
    while i <= length(args)
        a = args[i]
        startswith(a, "--") || error("unexpected argument $a")
        key = a[3:end]
        if key == "no-viz"
            opts[key] = true
            i += 1
        else
            i + 1 <= length(args) || error("missing value for $a")
            opts[key] = args[i+1]
            i += 2
        end
    end
    return opts
end

function main(args)
    opts = parseArgs(args)
    p = defaultParams(
        n=parse(Int, get(opts, "n", "65")),
        dt=parse(Float64, get(opts, "dt", "0.002")),
        t0=parse(Float64, get(opts, "t0", "0.05")),
        tf=parse(Float64, get(opts, "tf", "0.20")),
        gd=parse(Float64, get(opts, "gd", "100")),
        gg=parse(Float64, get(opts, "gg", "0.1")),
        out=get(opts, "out", joinpath(@__DIR__, "output")))
    result = simulate(p)
    if !get(opts, "no-viz", false)
        viz = joinpath(@__DIR__, "visualize.PlanarStefan.jl")
        if isfile(viz)
            include(viz)
            Base.invokelatest(visualizePlanarStefan, p.out)
        end
    end
    return result
end

if abspath(PROGRAM_FILE) == @__FILE__
    main(ARGS)
end
