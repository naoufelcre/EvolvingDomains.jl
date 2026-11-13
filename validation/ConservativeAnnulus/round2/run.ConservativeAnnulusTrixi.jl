# ConservativeAnnulusTrixi — four-method transport comparison on a translating cut-cell annulus.
#
# Campaign II, Experiment 1. One evolving geometry (a level-set `Field` advanced by WENO5),
# four independent scalar transport methods on the same background grid:
#
#   1. EvolvingDomains WENO5         eulerianTransport!               (Eulerian, in place)
#   2. EvolvingDomains CIP           lagrangianTransport!(t, s, v, dt) (semi-Lagrangian, intensive)
#   3. EvolvingDomains TransportMap  lagrangianTransport!(t, s, map)   (conservative cut-to-cut)
#   4. Trixi.jl DGSEM degree 3       SemidiscretizationHyperbolic + SimpleSSPRK33
#
# All four advance with the same dt and the same SSP-RK3 family. Every physical quantity is
# measured with one common rule: the material region is the *reconstructed* cut Ω_h(t) and a
# degree-5 triangle rule is applied to the same triangles for every method.
#
# Physical mass and the unweighted nodal sum are kept strictly apart:
#   M_h(t) = Σ_i w_i(t) ρ_i(t),  w_i(t) = ∫_{Ω_h(t)} N_i dx   — the physical integral
#   Σ_i ρ_i                                                   — the map's algebraic invariant
#
# Run:  julia validation/ConservativeAnnulus/round2/run.ConservativeAnnulusTrixi.jl [pilot]

using Pkg
Pkg.activate(joinpath(@__DIR__, "..", "..", "env"); io=devnull)

# Native cut backend: this (unreleased) checkout publishes no artifact, so a library built
# from the repository's own scripts/build_cxx_bindings.sh is shipped in lib/.
const LOCAL_CUT_LIBRARY = joinpath(@__DIR__, "lib",
    "libevolvingdomains_cut_cells" * (Sys.iswindows() ? ".dll" : Sys.isapple() ? ".dylib" : ".so"))
if !haskey(ENV, "EVOLVINGDOMAINS_CUTCELLS_LIBRARY") && isfile(LOCAL_CUT_LIBRARY)
    ENV["EVOLVINGDOMAINS_CUTCELLS_LIBRARY"] = LOCAL_CUT_LIBRARY
end

using EvolvingDomains
using StaticArrays: SVector
using Trixi
using Serialization
using Printf
using LinearAlgebra: dot

# ---------------------------------------------------------------------------------------
# Stage A — provenance. The resolved EvolvingDomains source must be this repository.
# ---------------------------------------------------------------------------------------

const EXPECTED_ROOT = normpath(joinpath(@__DIR__, "..", "..", ".."))
const OUTDIR = joinpath(@__DIR__, "output")

function provenance()
    edpath = normpath(pathof(EvolvingDomains))
    startswith(edpath, joinpath(EXPECTED_ROOT, "src")) || error(
        "Wrong EvolvingDomains checkout resolved: $edpath (expected $EXPECTED_ROOT). Aborting.")
    rev = strip(read(`git -C $EXPECTED_ROOT rev-parse HEAD`, String))
    dirty = split(strip(read(`git -C $EXPECTED_ROOT status --porcelain`, String)), '\n';
                  keepempty=false)
    native = get(ENV, "EVOLVINGDOMAINS_CUTCELLS_LIBRARY", "<none>")
    pin = match(r"GIT_TAG ([0-9a-f]+)",
                read(joinpath(EXPECTED_ROOT, "cxxBindings", "CMakeLists.txt"), String))
    cutcells_commit = pin === nothing ? "unknown" : pin.captures[1]
    info = (; julia=string(VERSION), root=EXPECTED_ROOT, source=edpath, revision=rev,
            dirty_files=length(dirty), version=string(pkgversion(EvolvingDomains)),
            cutcells_library=native, cutcells_bytes=isfile(native) ? filesize(native) : 0,
            cutcells_commit,
            trixi_version=string(pkgversion(Trixi)), trixi_source=normpath(pathof(Trixi)))
    @printf("Julia                 %s\n", info.julia)
    @printf("EvolvingDomains src   %s\n", info.source)
    @printf("EvolvingDomains git   %s (%d modified paths)\n", info.revision, info.dirty_files)
    @printf("EvolvingDomains ver   %s\n", info.version)
    @printf("CutCells native lib   %s (%d bytes)\n", info.cutcells_library, info.cutcells_bytes)
    @printf("CutCells source pin   %s\n", info.cutcells_commit)
    @printf("Trixi                 %s (%s)\n", info.trixi_version, info.trixi_source)
    return info
end

# ---------------------------------------------------------------------------------------
# §2 — analytic problem, with the periodic extension the Trixi setup requires.
# ---------------------------------------------------------------------------------------

const C0 = (0.35, 0.45)          # initial annulus centre
const RIN, ROUT = 0.12, 0.24     # inner / outer radius
const VEL = (0.24, 0.10)         # prescribed translating velocity
const TEND = 1.0
const BUMP = (0.41, 0.48)        # density bump centre in the initial frame
const WIDTH = 0.055              # density bump width

# Periodic extension of the Gaussian: the minimal image plus its eight neighbours. On the
# unit torus the nearest image is at distance >= 0.35, so the neighbour terms are <= 4e-18;
# they are summed explicitly anyway so that the initial datum, the analytic reference and
# the Trixi periodic solve all refer to the same torus-continuous function.
function gaussianPeriodic(dx, dy)
    s = 0.0
    for mx in -1:1, my in -1:1
        s += exp(-((dx + mx)^2 + (dy + my)^2) / WIDTH^2)
    end
    return s
end

rhoExact(x, y, t) = 1 + 0.5 * gaussianPeriodic(x - BUMP[1] - t * VEL[1],
                                              y - BUMP[2] - t * VEL[2])

radius(x, y, t) = hypot(x - C0[1] - t * VEL[1], y - C0[2] - t * VEL[2])

# φ < 0 is the annular material (§2).
phiExact(x, y, t) = max(RIN - radius(x, y, t), radius(x, y, t) - ROUT)

exactArea() = π * (ROUT^2 - RIN^2)

# Reference mass M(0) = ∫_{Ω(0)} ρ₀: 5-point Gauss–Legendre × 1024 midpoint sweep in θ.
const GL5 = ((0.0, 0.5688888888888889),
             (0.5384693101056831, 0.4786286704993665), (-0.5384693101056831, 0.4786286704993665),
             (0.9061798459386640, 0.2369268850561891), (-0.9061798459386640, 0.2369268850561891))

function exactMass(; nr=20, nθ=1024)
    total = 0.0
    for i in 1:nr
        a = RIN + (ROUT - RIN) * (i - 1) / nr
        b = RIN + (ROUT - RIN) * i / nr
        c, half = (a + b) / 2, (b - a) / 2
        for (ξ, w) in GL5
            r = c + half * ξ
            weight = half * w * r * (2π / nθ)
            for j in 0:(nθ - 1)
                θ = 2π * (j + 0.5) / nθ
                total += weight * rhoExact(C0[1] + r * cos(θ), C0[2] + r * sin(θ), 0.0)
            end
        end
    end
    return total
end

# ---------------------------------------------------------------------------------------
# §4 — one common integration rule: 7-point degree-5 Dunavant rule on each cut triangle.
# Quadrature = area × Σ w_q f(x_q). Exact for total degree <= 5, which covers the bilinear
# nodal reconstruction (degree 2) and the degree-3 DG polynomial exactly.
# ---------------------------------------------------------------------------------------

const DUN = (0.225,
             0.13239415278850617, 0.13239415278850617, 0.13239415278850617,
             0.12593918054482715, 0.12593918054482715, 0.12593918054482715)

const DUNBARY = ((1 / 3, 1 / 3, 1 / 3),
                 (0.0597158717897698, 0.4701420641051151, 0.4701420641051151),
                 (0.4701420641051151, 0.0597158717897698, 0.4701420641051151),
                 (0.4701420641051151, 0.4701420641051151, 0.0597158717897698),
                 (0.7974269853530873, 0.1012865073234564, 0.1012865073234564),
                 (0.1012865073234564, 0.7974269853530873, 0.1012865073234564),
                 (0.1012865073234564, 0.1012865073234564, 0.7974269853530873))

# ---------------------------------------------------------------------------------------
# Cut geometry — documented `Cut` arrays only.
# ---------------------------------------------------------------------------------------

"""Nodal support of a cut: a node is a corner of some in-phase or straddling background cell."""
function supportMask(cut, n::Int)
    nx = n - 1
    mask = falses(n * n)
    for j in 1:nx, i in 1:nx
        state = cut.bg_states[i + (j - 1) * nx]
        (state == CUT_IN || state == CUT_CUT) || continue
        for jj in (j, j + 1), ii in (i, i + 1)
            mask[ii + (jj - 1) * n] = true
        end
    end
    return mask
end

subcellVertices(cut, s::Int) = begin
    k = Int(cut.point_pointers[s])
    (cut.physical_points[cut.point_ids[k]], cut.physical_points[cut.point_ids[k + 1]],
     cut.physical_points[cut.point_ids[k + 2]])
end

triangleArea(p1, p2, p3) =
    0.5 * abs((p2[1] - p1[1]) * (p3[2] - p1[2]) - (p3[1] - p1[1]) * (p2[2] - p1[2]))

"""In-phase triangles as `(area, x1,y1,x2,y2,x3,y3, i, j)` records.

The parent background cell `(i, j)` is identified from the triangle's own bounding box: a
cut subcell never crosses a background cell boundary, so this needs no assumption about the
ordering of `cell_to_bgcell`. (`shortcut: bounding-box parent lookup, revisit only if
subcells can span background cells`.)"""
function inPhaseTriangles(cut, axis)
    h = axis[2] - axis[1]
    nx = length(axis) - 1
    lo = axis[1]
    out = NTuple{9,Float64}[]
    for s in eachindex(cut.subcell_states)
        cut.subcell_states[s] == CUT_IN || continue
        p1, p2, p3 = subcellVertices(cut, s)
        a = triangleArea(p1, p2, p3)
        a > 0 || continue
        xmin = min(p1[1], p2[1], p3[1])
        ymin = min(p1[2], p2[2], p3[2])
        i = clamp(floor(Int, (xmin - lo) / h) + 1, 1, nx)
        j = clamp(floor(Int, (ymin - lo) / h) + 1, 1, nx)
        push!(out, (a, p1[1], p1[2], p2[1], p2[2], p3[1], p3[2],
                    Float64(i), Float64(j)))
    end
    return out
end

"""Nodal integration weights w_i(t) = ∫_{Ω_h(t)} N_i(x) dx.

Each in-phase triangle is quadratured with the common degree-5 rule; on its parent cell the
four bilinear hat functions are evaluated in reference coordinates and accumulated into the
cell's corner nodes. Σ_i w_i = |Ω_h| because Σ_i N_i ≡ 1 — asserted by the caller."""
function nodalWeights(cut, axis)
    n = length(axis)
    h = axis[2] - axis[1]
    w = zeros(n * n)
    for tri in inPhaseTriangles(cut, axis)
        a, x1, y1, x2, y2, x3, y3, fi, fj = tri
        i, j = Int(fi), Int(fj)
        base = i + (j - 1) * n
        for q in 1:7
            λ1, λ2, λ3 = DUNBARY[q]
            x = λ1 * x1 + λ2 * x2 + λ3 * x3
            y = λ1 * y1 + λ2 * y2 + λ3 * y3
            ξ = (x - axis[i]) / h
            η = (y - axis[j]) / h
            c = a * DUN[q]
            w[base] += c * (1 - ξ) * (1 - η)
            w[base + 1] += c * ξ * (1 - η)
            w[base + n] += c * (1 - ξ) * η
            w[base + n + 1] += c * ξ * η
        end
    end
    return w
end

"""Bilinear value of a nodal array at a point inside background cell `(i, j)`."""
@inline function bilinear(rho, axis, h, i, j, x, y)
    n = length(axis)
    ξ = (x - axis[i]) / h
    η = (y - axis[j]) / h
    base = i + (j - 1) * n
    return (1 - ξ) * (1 - η) * rho[base] + ξ * (1 - η) * rho[base + 1] +
           (1 - ξ) * η * rho[base + n] + ξ * η * rho[base + n + 1]
end

# ---------------------------------------------------------------------------------------
# §3D — Trixi adapter. Only `get_node_coords`/`get_node_vars`/`wrap_array` are used, which is
# the documented analysis interface of the examples; no private storage is read.
#
# The DG solution *is* the tensor-product Lagrange interpolant of its LGL nodal values, so
# evaluating that interpolant reproduces the DG polynomial exactly. The adapter derives the
# LGL abscissae and the element size from the reported node coordinates.
# ---------------------------------------------------------------------------------------

struct DGView{U<:AbstractArray{Float64,4}}
    u::U                     # live (1, nn, nn, nelem) view of the Trixi state
    xi::Vector{Float64}      # LGL abscissae on [-1, 1]
    xl::Vector{Float64}      # element origins
    yl::Vector{Float64}
    ex::Matrix{Int}          # (ix, iy) -> element index
    xmin::Float64
    ymin::Float64
    len::Float64
    h::Float64
    ndiv::Int
end

function dgView(semi, uflat)
    mesh, equations, dg, cache = Trixi.mesh_equations_solver_cache(semi)
    nc = cache.elements.node_coordinates
    ne = Trixi.nelements(dg, cache)
    nn = Trixi.nnodes(dg)
    uw = Trixi.wrap_array(uflat, semi)          # live view, not a copy
    xl = zeros(ne)
    yl = zeros(ne)
    for e in 1:ne
        p = Trixi.get_node_coords(nc, equations, dg, 1, 1, e)
        xl[e] = p[1]
        yl[e] = p[2]
    end
    xlo = Trixi.get_node_coords(nc, equations, dg, 1, 1, 1)[1]
    xhi = Trixi.get_node_coords(nc, equations, dg, nn, 1, 1)[1]
    ylo = Trixi.get_node_coords(nc, equations, dg, 1, 1, 1)[2]
    elen = xhi - xlo
    xi = [2 * (Trixi.get_node_coords(nc, equations, dg, i, 1, 1)[1] - xlo) / elen - 1
          for i in 1:nn]
    ndiv = round(Int, 1 / elen)
    ex = zeros(Int, ndiv, ndiv)
    for e in 1:ne
        ex[round(Int, (xl[e] - xlo) / elen) + 1, round(Int, (yl[e] - ylo) / elen) + 1] = e
    end
    return DGView(uw, xi, xl, yl, ex, xlo, ylo, 1.0, elen, ndiv)
end

@inline function lagrange1d!(w, ξ, r)
    n = length(ξ)
    @inbounds for i in 1:n
        v = 1.0
        for j in 1:n
            j == i && continue
            v *= (r - ξ[j]) / (ξ[i] - ξ[j])
        end
        w[i] = v
    end
    return w
end

"""Evaluate the DG polynomial at a physical point of the periodic unit square."""
function evalDG(v::DGView, x, y; wx=zeros(length(v.xi)), wy=zeros(length(v.xi)))
    ax = mod(x - v.xmin, v.len)
    ay = mod(y - v.ymin, v.len)
    ix = clamp(floor(Int, ax / v.h) + 1, 1, v.ndiv)
    iy = clamp(floor(Int, ay / v.h) + 1, 1, v.ndiv)
    e = v.ex[ix, iy]
    lagrange1d!(wx, v.xi, 2 * (ax - v.xl[e]) / v.h - 1)
    lagrange1d!(wy, v.xi, 2 * (ay - v.yl[e]) / v.h - 1)
    s = 0.0
    @inbounds for j in eachindex(wy), i in eachindex(wx)
        s += wx[i] * wy[j] * v.u[1, i, j, e]
    end
    return s
end

# ---------------------------------------------------------------------------------------
# §4 — diagnostics on the common rule.
# ---------------------------------------------------------------------------------------

"""Per-method integrals over Ω_h: physical mass, relative L1 and L2 against the reference.

The three nodal arrays are read through their bilinear reconstruction on the parent cell and
the DG field through its polynomial; both are integrated on the same triangles with the same
degree-5 rule, so no method is favoured by its quadrature."""
function commonIntegrals(cut, axis, fields, dg, ref, t)
    n = length(axis)
    h = axis[2] - axis[1]
    # (mass, l1, l2) for the 3 nodal methods, then the DG method
    mass = zeros(4)
    l1 = zeros(4)
    l2 = zeros(4)
    ref_mass = 0.0
    ref_l2 = 0.0
    area = 0.0
    mass_centroid = 0.0
    wx = zeros(length(dg.xi))
    wy = zeros(length(dg.xi))
    for tri in inPhaseTriangles(cut, axis)
        a, x1, y1, x2, y2, x3, y3, fi, fj = tri
        i, j = Int(fi), Int(fj)
        area += a
        mass_centroid += a * ref((x1 + x2 + x3) / 3, (y1 + y2 + y3) / 3, t)
        for q in 1:7
            λ1, λ2, λ3 = DUNBARY[q]
            x = λ1 * x1 + λ2 * x2 + λ3 * x3
            y = λ1 * y1 + λ2 * y2 + λ3 * y3
            wq = a * DUN[q]
            r = ref(x, y, t)
            ref_mass += wq * r
            ref_l2 += wq * r^2
            for (k, rho) in enumerate(fields)
                v = bilinear(rho, axis, h, i, j, x, y)
                mass[k] += wq * v
                l1[k] += wq * abs(v - r)
                l2[k] += wq * (v - r)^2
            end
            v = evalDG(dg, x, y; wx=wx, wy=wy)
            mass[4] += wq * v
            l1[4] += wq * abs(v - r)
            l2[4] += wq * (v - r)^2
        end
    end
    return (; area,
            mass=NTuple{4,Float64}(mass),
            l1=NTuple{4,Float64}(l1 ./ ref_mass),
            l2=NTuple{4,Float64}(sqrt.(l2 ./ ref_l2)),
            mass_centroid, ref_mass)
end

"""Final-time nodal error binned by signed distance in node units, φ/h."""
function layerBins(qs, exactv, phin, h)
    edges = (-Inf, -4.0, -2.0, -1.0, 0.0, 1.0, 2.0, Inf)
    rows = Vector{NTuple{4,Float64}}[]     # (count, linf, l1mean, l2mean)
    for m in 1:4
        push!(rows, NTuple{4,Float64}[])
        for b in 1:(length(edges) - 1)
            sel = (phin .>= edges[b] * h) .& (phin .< edges[b + 1] * h)
            cnt = count(sel)
            if cnt == 0
                push!(rows[m], (0.0, 0.0, 0.0, 0.0))
            else
                e = abs.(qs[m][sel] .- exactv[sel])
                push!(rows[m], (Float64(cnt), maximum(e), sum(e) / cnt,
                                sqrt(sum(e .^ 2) / cnt)))
            end
        end
    end
    return rows
end

# ---------------------------------------------------------------------------------------
# §6 — statistics, advances, and the simulation loop.
# ---------------------------------------------------------------------------------------

struct Stats
    time::Base.RefValue{Float64}
    bytes::Base.RefValue{Int}
    calls::Base.RefValue{Int}
end
Stats() = Stats(Ref(0.0), Ref(0), Ref(0))
stuple(s::Stats) = (time=s.time[], bytes=s.bytes[], calls=s.calls[])

macro measure(stats, expr)
    quote
        local s = $(esc(stats))
        local r = @timed $(esc(expr))
        s.time[] += r.time
        s.bytes[] += r.bytes
        s.calls[] += 1
        r.value
    end
end

advanceGeometry!(levelset, velocity, dt, level, cache) =
    eulerianTransport!(levelset, velocity, dt; level_id=level, cache=cache)

advanceWithWENO!(rho, velocity, dt, level, cache) =
    eulerianTransport!(rho, velocity, dt; level_id=level, cache=cache)

function advanceWithCIP!(dst, src, velocity, dt, level, cache)
    lagrangianTransport!(dst, src, velocity, dt; level_id=level, cache=cache)
    return nothing
end

function advanceWithMap!(target, source, map)
    lagrangianTransport!(target, source, map)
    return nothing
end

# One SSP-RK3 step of Trixi's own integrator. The integrator is constructed with the timestep
# shared by all four methods (`Trixi.init(...; dt)`), so no integrator field is written here.
advanceWithTrixi!(integrator) = Trixi.step!(integrator)

# CIP refuses outer-boundary inflow: `_departureCell` throws when a departure point leaves the
# box, and the public operator offers no periodic or inflow boundary option (the whole grid is
# checked, whatever the values are). With a constant translation and vx, vy > 0 only the two
# lower faces can violate that; every point the diagnostics touch satisfies
# x - vx t >= 0.11 and y - vy t >= 0.21 for t <= 1, so holding the advection stationary in the
# single outermost node layer cannot reach the material region (the CIP stencil is 2 cells
# wide; that layer is >= 7 cells away even on the coarse grid). The WENO5 and TransportMap
# methods use the untampered velocity; this one deviation is CIP-specific and documented.
function cipVelocity(axis)
    n = length(axis)
    v = Vector{SVector{2,Float64}}(undef, n * n)
    for j in 1:n, i in 1:n
        v[i + (j - 1) * n] = SVector(i == 1 ? 0.0 : VEL[1], j == 1 ? 0.0 : VEL[2])
    end
    return v
end

const SNAPSHOT_TIMES = (0.0, 0.5, 1.0)

"""Trixi semidiscretization: degree-3 DGSEM, `nside²` elements, periodic unit square."""
function trixiSetup(nside, rho0; dt=1e-2)
    equations = LinearScalarAdvectionEquation2D(VEL)
    solver = DGSEM(polydeg=3, surface_flux=flux_lax_friedrichs)
    mesh = TreeMesh((0.0, 0.0), (1.0, 1.0); initial_refinement_level=round(Int, log2(nside)),
                    periodicity=true)
    ic(x, t, eq) = SVector(rho0(x[1], x[2], t))   # torus-continuous datum
    semi = SemidiscretizationHyperbolic(mesh, equations, ic, solver;
                                        boundary_conditions=boundary_condition_periodic)
    ode = semidiscretize(semi, (0.0, TEND))
    dtmax = Trixi.calculate_dt(ode.u0, 0.0, 1.0, 0.0, semi)
    integrator = Trixi.init(ode, Trixi.SimpleSSPRK33(); dt=dt)
    return (; semi, solver, equations, ode, integrator, dtmax,
            nelements=Trixi.nelements(solver, semi.cache),
            ndof=Trixi.nelements(solver, semi.cache) * Trixi.nnodes(solver)^2)
end

function simulate(d, level, nside, label; cfl=0.4, tfinal=TEND, constant=false)
    grid = d.levels[level]
    n = length(grid.axes[1])
    axis = grid.axes[1]
    nn = n * n
    h = axis[2] - axis[1]

    nsteps = max(1, ceil(Int, tfinal / (cfl * h / (abs(VEL[1]) + abs(VEL[2])))))
    dt = tfinal / nsteps
    ref = constant ? ((x, y, t) -> 1.0) : rhoExact

    levelset = Field(d, Dict(level => [phiExact(x, y, 0.0) for x in axis, y in axis]))
    rhoweno = Field(d, Dict(level => [ref(x, y, 0.0) for x in axis, y in axis]))
    rhoa = Field(d, Dict(level => [ref(x, y, 0.0) for x in axis, y in axis]))
    rhob = Field(d, Dict(level => [ref(x, y, 0.0) for x in axis, y in axis]))
    rhomap = [ref(x, y, 0.0) for x in axis, y in axis]

    velocity = vec([SVector(VEL[1], VEL[2]) for _ in axis, _ in axis])
    cipvel = cipVelocity(axis)
    phicache = WENO5Cache()
    rhocache = WENO5Cache()
    cipcache = CIPCache(vec(rhoa.data[level]), grid)

    cipdev = maximum(maximum(abs.(c - v)) for (c, v) in zip(cipvel, velocity))
    trixi = trixiSetup(nside, ref; dt)
    dt <= trixi.dtmax || error("shared dt $dt exceeds Trixi's stable dt $(trixi.dtmax)")
    dg = dgView(trixi.semi, trixi.integrator.u)
    mesh, equations, solver, cache = Trixi.mesh_equations_solver_cache(trixi.semi)
    uw = Trixi.wrap_array(trixi.integrator.u, trixi.semi)

    stGeo, stWeno, stCip, stMap, stTrixi, stAdapter = Stats(), Stats(), Stats(), Stats(),
                                                       Stats(), Stats()

    # Compile every operator on throwaway state, so the reported work measures the methods and
    # not JIT compilation (the first call of each path otherwise dominates its total).
    @timed begin
        wphi = Field(d, Dict(level => copy(levelset.data[level])))
        advanceGeometry!(wphi, velocity, dt, level, WENO5Cache())
        wsrc = Field(d, Dict(level => copy(rhoweno.data[level])))
        wdst = Field(d, Dict(level => copy(rhoweno.data[level])))
        advanceWithWENO!(wsrc, velocity, dt, level, WENO5Cache())
        advanceWithCIP!(wdst, wsrc, cipvel, dt, level, CIPCache(vec(wsrc.data[level]), grid))
        wmap = zeros(nn)
        icut = cut(levelset, level)
        advanceWithMap!(wmap, vec(rhoweno.data[level]), TransportMap(icut, icut, grid, velocity, dt))
        wtrixi = Trixi.init(trixi.ode, Trixi.SimpleSSPRK33(); dt=dt)
        advanceWithTrixi!(wtrixi)
        advanceWithTrixi!(wtrixi)
        nodalWeights(icut, axis)
        commonIntegrals(icut, axis, (vec(rhoweno.data[level]), vec(rhoa.data[level]),
                                     vec(rhomap)), dg, ref, 0.0)
        layerBins((vec(rhomap), vec(rhomap), vec(rhomap), vec(rhomap)),
                  vec(rhomap), vec(levelset.data[level]), h)
    end

    times = Float64[]
    areas = Float64[]
    masses = Vector{NTuple{4,Float64}}()
    nodalsums = Vector{NTuple{3,Float64}}()
    linfcore = Vector{NTuple{4,Float64}}()
    rmins = Vector{NTuple{4,Float64}}()
    rmaxs = Vector{NTuple{4,Float64}}()
    mapdrift = Float64[]
    weightdrift = Float64[]
    dgfullnorm = Float64[]
    mass_centroid = Float64[]
    mass_exactnodal = Float64[]
    snaps = NamedTuple[]

    function diagnostics(t, cutcurrent, rhoCip)
        fields = (vec(rhoweno.data[level]), vec(rhoCip), vec(rhomap))
        ci = @measure stAdapter commonIntegrals(cutcurrent, axis, fields, dg, ref, t)
        w = @measure stAdapter nodalWeights(cutcurrent, axis)
        # Σ w_i must reproduce the cut area: independent check of the weight construction.
        @assert abs(sum(w) / ci.area - 1) < 1e-12 "nodal weights do not integrate to the cut area"
        mass = (dot(w, fields[1]), dot(w, fields[2]), dot(w, fields[3]), ci.mass[4])
        nodal = (sum(fields[1]), sum(fields[2]), sum(fields[3]))
        phin = vec(levelset.data[level])
        inside = phin .< 0
        core = phin .< -h
        exactv = vec([ref(x, y, t) for x in axis, y in axis])
        tval = @measure stAdapter begin
            local out = zeros(nn)
            local wx = zeros(length(dg.xi))
            local wy = zeros(length(dg.xi))
            for j in 1:n, i in 1:n
                out[i + (j - 1) * n] = evalDG(dg, axis[i], axis[j]; wx, wy)
            end
            out
        end
        qs = (fields[1], fields[2], fields[3], tval)
        linf = ntuple(k -> maximum(abs.(qs[k][core] .- exactv[core])), 4)
        rmin = ntuple(k -> minimum(qs[k][inside]), 4)
        rmax = ntuple(k -> maximum(qs[k][inside]), 4)
        return (; area=ci.area, mass, l1=ci.l1, l2=ci.l2, nodal, linf, rmin, rmax, w, tval,
                mass_centroid=ci.mass_centroid, ref_mass=ci.ref_mass, qs, exactv, phin,
                mass_exactnodal=dot(w, exactv),
                exactrange=(minimum(exactv[inside]), maximum(exactv[inside])))
    end

    function record!(t, cutcurrent, rhoCip; snapshot=false)
        di = diagnostics(t, cutcurrent, rhoCip)
        push!(times, t)
        push!(areas, di.area)
        push!(masses, di.mass)
        push!(nodalsums, di.nodal)
        push!(linfcore, di.linf)
        push!(rmins, di.rmin)
        push!(rmaxs, di.rmax)
        push!(mass_centroid, di.mass_centroid)
        push!(mass_exactnodal, di.mass_exactnodal)
        push!(dgfullnorm, Trixi.integrate((u, eq) -> u[1], uw, mesh, equations, solver, cache;
                                          normalize=false))
        if snapshot
            push!(snaps, (; t, phi=copy(levelset.data[level]),
                          weno=copy(rhoweno.data[level]), cip=copy(rhoCip),
                          map=copy(rhomap), trixi=reshape(di.tval, n, n),
                          triangles=[(q[2:3], q[4:5], q[6:7])
                                     for q in inPhaseTriangles(cutcurrent, axis)],
                          l1=di.l1, l2=di.l2, mass=di.mass, area=di.area,
                          exactrange=di.exactrange,
                          bins=layerBins(di.qs, di.exactv, di.phin, h)))
        end
        return di
    end

    # t = 0: geometry is the analytic initial datum (initialisation only, §3A)
    cutcurrent = cut(levelset, level)
    record!(0.0, cutcurrent, rhoa.data[level]; snapshot=true)

    t = 0.0
    for step in 1:nsteps
        # 2. old cut and the map's source vector (zeroed outside the old support, which is the
        #    documented contract of the conservative map)
        old_cut = cut(levelset, level)
        old_support = supportMask(old_cut, n)
        srcmap = copy(vec(rhomap))
        srcmap[.!old_support] .= 0.0
        srcsum = sum(srcmap)

        # 3. advance the shared geometry numerically
        @measure stGeo advanceGeometry!(levelset, velocity, dt, level, phicache)
        t += dt

        # 4. advance each scalar method on its own state
        @measure stWeno advanceWithWENO!(rhoweno, velocity, dt, level, rhocache)
        @measure stCip advanceWithCIP!(rhob, rhoa, cipvel, dt, level, cipcache)
        rhoa.data[level], rhob.data[level] = rhob.data[level], rhoa.data[level]

        new_cut = cut(levelset, level)
        new_support = supportMask(new_cut, n)
        map = TransportMap(old_cut, new_cut, grid, velocity, dt)
        @assert map.source_mask == old_support "map source support differs from the cut mask"
        @assert map.target_mask == new_support "map target support differs from the cut mask"
        dstmap = zeros(nn)
        @measure stMap advanceWithMap!(dstmap, srcmap, map)
        rhomap .= reshape(dstmap, n, n)
        push!(mapdrift, (sum(dstmap) - srcsum) / max(abs(srcsum), eps()))

        # Is the map conservative for the *physical* measure? w_newᵀ P ?= w_oldᵀ.
        wold = @measure stAdapter nodalWeights(old_cut, axis)
        wnew = @measure stAdapter nodalWeights(new_cut, axis)
        wpush = zeros(nn)
        @measure stMap lagrangianTransport!(wpush, copy(wold), map)
        push!(weightdrift, dot(wnew, wpush) / dot(wold, wold) - 1)

        # 5. Trixi, on its own state
        @measure stTrixi advanceWithTrixi!(trixi.integrator)

        # 6. shared diagnostics on the reconstructed geometry
        cutcurrent = new_cut
        record!(t, cutcurrent, rhoa.data[level];
                snapshot=any(abs(t - ts) < 1e-9 for ts in SNAPSHOT_TIMES))
    end

    # Trixi must have reached the same physical time with the same number of steps.
    @assert trixi.integrator.iter == nsteps "Trixi took $(trixi.integrator.iter) steps, expected $nsteps"
    @assert abs(trixi.integrator.t - t) < 1e-12 "Trixi reached t=$(trixi.integrator.t), expected $t"

    return (; label, n, nside, h, dt, nsteps, ndof_ed=nn, ndof_trixi=trixi.ndof,
            trixi_t=trixi.integrator.t, trixi_steps=trixi.integrator.iter,
            nelements_trixi=trixi.nelements, dtmax_trixi=trixi.dtmax,
            times, areas, masses, nodalsums, linfcore, rmins, rmaxs, mapdrift, weightdrift,
            dgfullnorm,
            mass_centroid, mass_exactnodal, snaps, constant, cipdev,
            stats=(geo=stuple(stGeo), weno=stuple(stWeno), cip=stuple(stCip),
                   map=stuple(stMap), trixi=stuple(stTrixi), adapter=stuple(stAdapter)))
end

# ---------------------------------------------------------------------------------------
# Stage B — the DG adapter's polynomial reproduction, checked without any Trixi solve.
# ---------------------------------------------------------------------------------------

"""Feed Trixi a Q3 polynomial as its initial datum and check that `evalDG` reproduces it at
random points to machine precision. The DG solution of degree 3 *is* the tensor-product
Lagrange interpolant of its LGL values, so this is the exact statement the adapter claims."""
function verifyAdapter(; nside=8, npts=200, seed=7)
    poly(x, y) = 1.0 + 0.3x - 0.7y + 0.5x * y + 0.2x^3 - 0.4x^2 * y - 0.9x * y^2 +
                 0.1y^3 + 0.05x^3 * y^2 - 0.2x^2 * y^3 + 0.03x^3 * y^3
    tr = trixiSetup(nside, (x, y, t) -> 1.0)
    mesh, equations, dg, cache = Trixi.mesh_equations_solver_cache(tr.semi)
    u = Trixi.compute_coefficients((x, t, eq) -> SVector(poly(x[1], x[2])), 0.0, tr.semi)
    v = dgView(tr.semi, u)
    # deterministic pseudo-random points (no extra dependency)
    s = UInt64(seed)
    nextpoint() = begin
        s = s * 6364136223846793005 + 1442695040888963407
        Float64(s % UInt64(1 << 40)) / Float64(1 << 40)
    end
    worst = 0.0
    for _ in 1:npts
        x, y = nextpoint(), nextpoint()
        worst = max(worst, abs(evalDG(v, x, y) - poly(x, y)))
    end
    # sanity: the same evaluator on the true (non-polynomial) datum is only interpolation-accurate
    tr2 = trixiSetup(nside, rhoExact)
    mesh2, eq2, dg2, cache2 = Trixi.mesh_equations_solver_cache(tr2.semi)
    u2 = Trixi.compute_coefficients((x, t, eq) -> SVector(rhoExact(x[1], x[2], 0.0)), 0.0,
                                    tr2.semi)
    v2 = dgView(tr2.semi, u2)
    icerr = 0.0
    s = UInt64(seed)
    for _ in 1:npts
        x, y = nextpoint(), nextpoint()
        icerr = max(icerr, abs(evalDG(v2, x, y) - rhoExact(x, y, 0.0)))
    end
    return (; polydeg=3, nside, worst_polynomial_error=worst, ic_interpolation_error=icerr,
            npts)
end
# ---------------------------------------------------------------------------------------
# Stage D/E — the two resolutions, the constant-state test, and the report.
# ---------------------------------------------------------------------------------------

const METHOD_NAMES = ("WENO5", "CIP", "Map", "Trixi DG3")

function makeDiscretization(ncells::Int)
    axis = collect(range(0.0, 1.0; length=ncells + 1))
    return Discretization(Dict(:main => Grid((axis, copy(axis)))))
end

fmt(x) = @sprintf("%9.3e", x)

function reportCase(c, exact_area, exact_mass)
    mref = c.constant ? exact_area : exact_mass
    m0 = c.masses[1]
    n = length(c.times)
    @printf("\n=== %s: %d EvolvingDomains cells / %d Trixi elements (degree 3) ===\n",
            c.label, c.n - 1, c.nelements_trixi)
    @printf("nodes %d, ED dof %d, Trixi dof %d, steps %d, dt %.4e (Trixi dtmax %.4e)\n",
            c.n^2, c.ndof_ed, c.ndof_trixi, c.nsteps, c.dt, c.dtmax_trixi)
    @printf("Trixi advanced %d steps to t = %.15f (EvolvingDomains t = %.15f)\n",
            c.trixi_steps, c.trixi_t, c.times[end])
    @printf("exact area %.10e, measured at t=0 %.10e, t=1 %.10e\n",
            exact_area, c.areas[1], c.areas[end])
    worst_area = maximum(abs.(c.areas ./ exact_area .- 1))
    @printf("worst relative area error %.3e\n", worst_area)
    @printf("\n%-9s %13s %13s %13s %13s %13s %13s\n", "method", "M0/Mref-1", "M1/Mref-1",
            "max|drift|", "L1(1)", "L2(1)", "Linf core(1)")
    for k in 1:4
        drift = [c.masses[i][k] / mref - 1 for i in 1:n]
        @printf("%-9s %13.5e %13.5e %13.5e %13.5e %13.5e %13.5e\n", METHOD_NAMES[k],
                m0[k] / mref - 1, c.masses[end][k] / mref - 1, maximum(abs.(drift)),
                c.snaps[end].l1[k], c.snaps[end].l2[k], c.linfcore[end][k])
    end
    @printf("\ncontrols on the same cut: analytic ρ at centroids  M/Mref-1 = %.3e (t=1)\n",
            c.mass_centroid[end] / mref - 1)
    @printf("                         exact nodal values through w   M/Mref-1 = %.3e (t=1)\n",
            c.mass_exactnodal[end] / mref - 1)
    @printf("\nrelative L1 error by phi/h bin at t=1, printed as count/mean error:\n")
    edges = ("(-inf,-4)", "[-4,-2)", "[-2,-1)", "[-1,0)", "[0,1)", "[1,2)", "[2,inf)")
    @printf("%-9s", "method")
    for e in edges
        @printf(" %15s", e)
    end
    @printf("\n")
    for k in 1:4
        @printf("%-9s", METHOD_NAMES[k])
        for b in 1:7
            cnt, linf, l1m, l2m = c.snaps[end].bins[k][b]
            @printf(" %7.0f/%.1e", cnt, l1m)
        end
        @printf("\n")
    end
    @printf("\nnodal-sum drift (relative, worst over steps): WENO5 %.2e, CIP %.2e, Map %.2e\n",
            maximum(abs.([c.nodalsums[i][1] / c.nodalsums[1][1] - 1 for i in 1:n])),
            maximum(abs.([c.nodalsums[i][2] / c.nodalsums[1][2] - 1 for i in 1:n])),
            maximum(abs.([c.nodalsums[i][3] / c.nodalsums[1][3] - 1 for i in 1:n])))
    @printf("\ntransport-only mass change (M(1)-M(0))/M(0): ")
    for k in 1:4
        @printf("%s %.2e  ", METHOD_NAMES[k], c.masses[end][k] / c.masses[1][k] - 1)
    end
    @printf("\n")
    @printf("\ndensity range on material nodes at t=1 (exact %.4f .. %.4f):\n",
            minimum(c.snaps[end].exactrange), maximum(c.snaps[end].exactrange))
    for k in 1:4
        @printf("  %-9s [%.4f, %.4f]   min over steps %.4f\n", METHOD_NAMES[k],
                c.rmins[end][k], c.rmaxs[end][k], minimum(getindex.(c.rmins, k)))
    end
    if c.label != :constant
        @printf("\nTrixi full-domain integral: %.12e -> %.12e (rel %.3e)\n", c.dgfullnorm[1],
                c.dgfullnorm[end], c.dgfullnorm[end] / c.dgfullnorm[1] - 1)
        @printf("map unweighted nodal sum drift: worst %.3e (per step, rel)\n",
                maximum(abs.(c.mapdrift)))
        @printf("weighted identity w_new'P vs w_old': worst %.3e (per step, rel), t=1 %.3e\n",
                maximum(abs.(c.weightdrift)), c.weightdrift[end])
        @printf("unweighted nodal sums at t=1: WENO5 %.6e, CIP %.6e, Map %.6e\n",
                c.nodalsums[end][1], c.nodalsums[end][2], c.nodalsums[end][3])
    end
    @printf("\nwork (solver only, seconds / MiB):\n")
    for (nm, st) in (("geometry", c.stats.geo), ("WENO5", c.stats.weno),
                     ("CIP", c.stats.cip), ("Map", c.stats.map), ("Trixi DG3", c.stats.trixi),
                     ("adapter+diag", c.stats.adapter))
        @printf("  %-12s %8.3f s  %9.1f MiB  %d calls\n", nm, st.time, st.bytes / 2.0^20,
                st.calls)
    end
    return (; label=c.label, ncells=c.n - 1, nside=c.nside, ndof_ed=c.ndof_ed,
            ndof_trixi=c.ndof_trixi, nsteps=c.nsteps, dt=c.dt,
            area_error=worst_area,
            mass_drift=ntuple(k -> c.masses[end][k] / mref - 1, 4),
            mass_worst=ntuple(k -> maximum(abs.([c.masses[i][k] / mref - 1 for i in 1:n])), 4),
            mass0=ntuple(k -> m0[k] / mref - 1, 4),
            l1=Tuple(c.snaps[end].l1), l2=Tuple(c.snaps[end].l2),
            linf=Tuple(c.linfcore[end]),
            time=map(s -> getfield(c.stats, s).time, (:geo, :weno, :cip, :map, :trixi, :adapter)),
            bytes=map(s -> getfield(c.stats, s).bytes, (:geo, :weno, :cip, :map, :trixi, :adapter)),
            mass_centroid=c.mass_centroid[end] / mref - 1,
            mass_exactnodal=c.mass_exactnodal[end] / mref - 1,
            map_nodal_drift=maximum(abs.(c.mapdrift)),
            dgfull_drift=c.dgfullnorm[end] / c.dgfullnorm[1] - 1,
            rmin=ntuple(k -> minimum(getindex.(c.rmins, k)), 4),
            rmax=ntuple(k -> maximum(getindex.(c.rmaxs, k)), 4),
            linf_bins=c.snaps[end].bins, exact_range=c.snaps[end].exactrange)
end

function main(; pilot=false)
    prov = provenance()
    mkpath(OUTDIR)
    exact_area = exactArea()
    exact_mass = exactMass()
    @printf("\nexact area %.12e   exact mass M(0) %.12e\n", exact_area, exact_mass)

    ad = verifyAdapter()
    @printf("DG adapter polynomial reproduction: max |evalDG - poly| = %.3e over %d points\n",
            ad.worst_polynomial_error, ad.npts)
    ad.worst_polynomial_error < 1e-11 || error("DG adapter does not reproduce Q3 polynomials")

    cases = pilot ? ((:ncells64, 64, 16),) : ((:ncells64, 64, 16), (:ncells128, 128, 32))
    results = NamedTuple[]
    summary = NamedTuple[]
    for (label, ncells, nside) in cases
        d = makeDiscretization(ncells)
        c = simulate(d, :main, nside, label)
        push!(summary, reportCase(c, exact_area, exact_mass))
        push!(results, c)
        cc = simulate(d, :main, nside, Symbol(label, :_const); constant=true)
        push!(results, cc)
        push!(summary, reportCase(cc, exact_area, exact_mass))
    end

    serialize(joinpath(OUTDIR, "results.jls"),
              (; problem=(; C0, RIN, ROUT, VEL, TEND, BUMP, WIDTH),
                 exact_area, exact_mass, prov, adapter=ad,
                 summary, cases=results,
                 pilot))
    @printf("\nwrote %s\n", joinpath(OUTDIR, "results.jls"))
    return results
end

if abspath(PROGRAM_FILE) == @__FILE__
    main(pilot="pilot" in ARGS)
end
