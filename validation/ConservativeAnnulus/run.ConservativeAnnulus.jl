# ConservativeAnnulus — two-dimensional conservative cut-to-cut transport benchmark.
#
# Model (§2). A translating annular material domain Ω(t) = Ω(0) + t v under
# ∂_t ρ + ∇·(ρ v) = 0 in Ω(t). Because v is a constant translation the exact solution is
# ρ(x,t) = ρ₀(x - t v), the exact area is π(r_out² - r_in²) and the physical mass
# M(t) = ∫_{Ω(t)} ρ dx = M(0) is constant.
#
# The object under test is EvolvingDomains' conservative cut-to-cut map (Lentine,
# Grétarsson & Fedkiw 2011). Its claimed invariant is the *unweighted nodal sum* over the
# map's source support. This script keeps that apart from the physical integral on purpose:
#
#   q      nodal array of length prod(dims); Σ q is the operator's invariant.
#   M_h    ∫ρ over the *reconstructed* cut geometry, summed as (triangle area) × ρ(centroid)
#          over every CUT_IN subcell of the current Cut. This is the physical integral and it
#          does not know about Σ q.
#
# Nothing here reads a package internal: the cut ↔ array mapping uses only `bg_states`,
# `subcell_states`, `cell_to_bgcell`, `point_ids`/`point_pointers` and `physical_points`,
# all documented on `Cut`/`CutSurface`.
#
# Run:  julia validation/ConservativeAnnulus/run.ConservativeAnnulus.jl

using Pkg
Pkg.activate(joinpath(@__DIR__, "..", "env"); io = devnull)

# The native CutCells backend ships as a Julia artifact, which this (unreleased) checkout
# does not have. A locally built library is used unless the caller names one.
const LOCAL_CUT_LIBRARY = joinpath(@__DIR__, "lib",
    "libevolvingdomains_cut_cells" * (Sys.iswindows() ? ".dll" : Sys.isapple() ? ".dylib" : ".so"))
if !haskey(ENV, "EVOLVINGDOMAINS_CUTCELLS_LIBRARY") && isfile(LOCAL_CUT_LIBRARY)
    ENV["EVOLVINGDOMAINS_CUTCELLS_LIBRARY"] = LOCAL_CUT_LIBRARY
end

using EvolvingDomains
using StaticArrays: SVector
using Serialization
using Printf

# ---------------------------------------------------------------------------------------
# Analytic problem. C0, RIN, ROUT, VEL and the density parameters are the mission's.
# ---------------------------------------------------------------------------------------

const C0 = (0.35, 0.45)          # initial annulus centre
const RIN, ROUT = 0.12, 0.24     # inner / outer radius
const VEL = (0.24, 0.10)         # prescribed translating velocity
const TEND = 1.0
const BUMP = (0.41, 0.48)        # density bump centre, initial frame
const WIDTH = 0.055              # density bump width

radius(x, y, t) = hypot(x - C0[1] - t * VEL[1], y - C0[2] - t * VEL[2])

# φ < 0 is the annular material, as specified.
phiExact(x, y, t) = max(RIN - radius(x, y, t), radius(x, y, t) - ROUT)

# ρ(x,t) = ρ₀(x - t v).
function rhoExact(x, y, t)
    dx = x - BUMP[1] - t * VEL[1]
    dy = y - BUMP[2] - t * VEL[2]
    return 1 + 0.5 * exp(-(dx^2 + dy^2) / WIDTH^2)
end

exactArea() = π * (ROUT^2 - RIN^2)

# M(0) = ∫_{Ω(0)} ρ₀ dx in polar coordinates around C0: 5-point Gauss–Legendre on 20
# subintervals of [r_in, r_out] times a midpoint-rule θ sweep. Both rules converge spectrally
# on this smooth integrand, so the value is converged to ~1e-14 relative.
const GL5 = ((0.0, 0.5688888888888889),
             (0.5384693101056831, 0.4786286704993665), (-0.5384693101056831, 0.4786286704993665),
             (0.9061798459386640, 0.2369268850561891), (-0.9061798459386640, 0.2369268850561891))

function exactMass(; nr = 20, nθ = 1024)
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
# Cut → nodal support, and cut → physical integrals.
# ---------------------------------------------------------------------------------------

"""
    supportMask(cut, n)

Nodal support of a cut: a node belongs when it is a corner of some in-phase (`CUT_IN`) or
straddling (`CUT_CUT`) background cell. Linear index is `i + (j-1)*n`, first axis fastest,
which is exactly the flattening of the `(n, n)` Field array — so `vec(field.data[level])`
and this mask index the same entries. The cut's own point ordering is never involved.

`bg_states` is indexed by background cell in the same first-axis-fastest order.
"""
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

# The three physical vertices of subcell `s` (simplices have N+1 = 3 vertices in 2D).
function subcellVertices(cut, s::Int)
    k = Int(cut.point_pointers[s])
    return (cut.physical_points[cut.point_ids[k]],
            cut.physical_points[cut.point_ids[k + 1]],
            cut.physical_points[cut.point_ids[k + 2]])
end

triangleArea(p1, p2, p3) =
    0.5 * abs((p2[1] - p1[1]) * (p3[2] - p1[2]) - (p3[1] - p1[1]) * (p2[2] - p1[2]))

"""
    cutIntegrals(cut, value)

Physical integrals over the material region described by `cut`:

    area = Σ_K |K|,   centre = Σ_K |K| x_K / area,   mass = Σ_K |K| value(x_K)

summed over the CUT_IN triangular subcells, with `x_K` the triangle centroid. This is
centroid (one-point) quadrature of the density: exact for a linear density, second order for
a smooth one. It is a quadrature *weighted by physical cell areas*; the unweighted nodal sum
`Σ q` of the transport operator is a different quantity and is reported separately.
"""
function cutIntegrals(cut, value)
    area = 0.0
    mx = 0.0
    my = 0.0
    mass = 0.0
    for s in eachindex(cut.subcell_states)
        cut.subcell_states[s] == CUT_IN || continue
        p1, p2, p3 = subcellVertices(cut, s)
        a = triangleArea(p1, p2, p3)
        x = (p1[1] + p2[1] + p3[1]) / 3
        y = (p1[2] + p2[2] + p3[2]) / 3
        area += a
        mx += a * x
        my += a * y
        mass += a * value(x, y)
    end
    return (area = area, cx = mx / area, cy = my / area, mass = mass)
end

# Subcell triangles for the figures (vectorised so Makie draws them as a polygon list).
function cutTriangles(cut)
    tris = NTuple{3,NTuple{2,Float64}}[]
    for s in eachindex(cut.subcell_states)
        cut.subcell_states[s] == CUT_IN || continue
        push!(tris, subcellVertices(cut, s))
    end
    return tris
end

# Reconstructed interface: the cut's own boundary facets, as 2-point segments.
function cutSegments(cut)
    segs = NTuple{2,NTuple{2,Float64}}[]
    b = cut.boundary
    for f in eachindex(b.facet_to_bgcell)
        k = Int(b.point_pointers[f])
        push!(segs, (b.physical_points[b.point_ids[k]], b.physical_points[b.point_ids[k + 1]]))
    end
    return segs
end

# ---------------------------------------------------------------------------------------
# Multilevel study: hierarchy navigation and transfer! against a known analytic field.
# ---------------------------------------------------------------------------------------

smoothField(x, y) = exp(sin(2π * x) + cos(2π * y))

gridValues(grid, f) = [f(x, y) for x in grid.axes[1], y in grid.axes[2]]
gridSpacing(grid) = 1 / (length(grid.axes[1]) - 1)

"""
    transferStudy(d)

Prolongates the smooth analytic field `smoothField` along `:coarse → :fine → :finer` with the
public `transfer!` (default `linearInterpolator`) and measures the interpolation error
against the exact values at the target nodes. The two successive errors give the observed
order of the transfer. It also runs the restriction direction, which is plain node sampling
and therefore exact on nested nodes — evidence for the docstring's warning that it is *not* a
conservative averaging restriction.
"""
function transferStudy(d)
    stats = NamedTuple[]
    for (from, to) in ((:coarse, :fine), (:fine, :finer))
        gfrom, gto = d.levels[from], d.levels[to]
        source_values = gridValues(gfrom, smoothField)
        field = Field(d, Dict(from => source_values, to => zeros(size(gridValues(gto, smoothField)))))
        transfer!(field, from, to)
        err = field.data[to] .- gridValues(gto, smoothField)
        push!(stats, (; from, to,
            h_from = gridSpacing(gfrom), h_to = gridSpacing(gto),
            error_inf = maximum(abs, err), error_l1 = sum(abs, err) / length(err),
            error_fine = copy(err), values = copy(field.data[to])))
    end

    # opposite direction: fills :coarse from the (exact) fine nodal values
    gc = d.levels[:coarse]
    back = Field(d, Dict(:coarse => zeros(size(gridValues(gc, smoothField))),
                         :fine => gridValues(d.levels[:fine], smoothField)))
    transfer!(back, :fine, :coarse)
    restriction_inf = maximum(abs, back.data[:coarse] .- gridValues(gc, smoothField))

    return (;
        steps = stats,
        order = log2(stats[1].error_inf / stats[2].error_inf),
        restriction_inf,
        error_fine = stats[1].error_fine,
        fine_values = stats[1].values,
    )
end

# ---------------------------------------------------------------------------------------
# The simulation. Every mathematical step of §3 is visible in the loop.
# ---------------------------------------------------------------------------------------

"""
    simulate(d, level; cfl = 0.4, tfinal = TEND, snapshot_times = (0.0, 0.5, 1.0),
             label = string(level))

One level of the multilevel discretization `d`: cut-to-cut conservative transport of the
density on the translating annulus, with every diagnostic recorded per step.

Per step, following §3 exactly:

1. prescribed velocity (a frozen nodal field; constant here),
2. old cut from the level set,
3. `eulerianTransport!` on the level set (WENO5 + SSP-RK3),
4. new cut,
5. `TransportMap(old, new, grid, velocity, dt)`,
6. `lagrangianTransport!(dest_vector, source_vector, map)`,
7. diagnostics.

The transported vectors are the flattened `(n, n)` background node arrays, i.e.
`vec(field.data[level])`. The source vector is zeroed outside the old cut's nodal support:
the map transports only its own support and drops everything else (documented), and the
dropped entries are, by construction, nodes that lie outside the material region.
"""
function simulate(d, level; cfl = 0.4, tfinal = TEND,
                  snapshot_times = (0.0, 0.5, 1.0), label = string(level))
    grid = d.levels[level]
    n = length(grid.axes[1])
    nn = n * n
    x, y = grid.axes[1], grid.axes[2]

    # uniform grids: the CFL rate uses h and |v|_1
    h = 1 / (n - 1)
    rate = (abs(VEL[1]) + abs(VEL[2])) / h
    nsteps = max(1, ceil(Int, tfinal * rate / cfl))
    dt = tfinal / nsteps

    levelset = Field(d, Dict(level => [phiExact(xx, yy, 0.0) for xx in x, yy in y]))
    density = Field(d, Dict(level => [rhoExact(xx, yy, 0.0) for xx in x, yy in y]))
    phi = levelset.data[level]
    rho = density.data[level]

    velocity = vec([SVector(VEL[1], VEL[2]) for _ in x, _ in y])
    transport_cache = WENO5Cache()
    interior = (xx, yy) -> linearInterpolator((xx, yy), rho, grid)

    exact_area = exactArea()
    exact_mass = exactMass()
    exact_center = t -> (C0[1] + t * VEL[1], C0[2] + t * VEL[2])

    initial_cut = cut(levelset, level)
    geometry0 = cutIntegrals(initial_cut, interior)
    records = NamedTuple[]
    snapshots = NamedTuple[]
    push!(snapshots, (; t = 0.0, triangles = cutTriangles(initial_cut),
                      segments = cutSegments(initial_cut), phi = copy(phi), rho = copy(rho)))

    t = 0.0
    timing = @timed for step in 1:nsteps
        # 1. prescribed velocity: `velocity` above (frozen over the step)
        # 2. old cut, and the source vector that belongs to it
        old_cut = cut(levelset, level)
        old_support = supportMask(old_cut, n)
        source = copy(vec(rho))
        source[.!old_support] .= 0.0
        source_sum = sum(source)

        # 3. move the level set
        eulerianTransport!(levelset, velocity, dt; level_id = level, cache = transport_cache)
        t += dt

        # 4. new cut
        new_cut = cut(levelset, level)
        new_support = supportMask(new_cut, n)

        # 5. conservative map for this step
        map = TransportMap(old_cut, new_cut, grid, velocity, dt)
        @assert map.source_mask == old_support "map source support differs from the cut mask"
        @assert map.target_mask == new_support "map target support differs from the cut mask"

        # 6. remap the density old-cut → new-cut degrees of freedom
        target = zeros(nn)
        lagrangianTransport!(target, source, map)
        rho .= reshape(target, n, n)
        target_sum = sum(target)
        @assert isapprox(target_sum, source_sum; rtol = 1e-12, atol = 0.0) ||
                source_sum == 0.0 "conservative map lost the nodal sum"

        # 7. diagnostics
        geometry = cutIntegrals(new_cut, interior)
        # counterfactual: the same quadrature with the void-side (φ > 0) nodes replaced by
        # their exact values. Separates material-side error from pollution of the shared cut
        # cells, whose corners are in-phase nodes on one side and void-side nodes on the other.
        rho_material = copy(rho)
        rho_material[phi .> 0] .= [rhoExact(xx, yy, t) for xx in x, yy in y][phi .> 0]
        material_mass = cutIntegrals(new_cut, (xx, yy) ->
            linearInterpolator((xx, yy), rho_material, grid)).mass
        geometry_exact_field = cutIntegrals(new_cut, (xx, yy) -> rhoExact(xx, yy, t))
        phin = vec(phi)
        core = new_support .&& (phin .< -h)             # one lattice step inside the phase
        inside = new_support .&& (phin .< 0)
        exact_nodal = [rhoExact(xx, yy, t) for xx in x, yy in y]
        exactv = vec(exact_nodal)
        nodal_error = target .- exactv

        # constant-state preservation, measured on the same map
        probe_dest = zeros(nn)
        lagrangianTransport!(probe_dest, Float64.(old_support), map)

        cexact = exact_center(t)
        push!(records, (;
            t = t,
            source_sum = source_sum,
            target_sum = target_sum,
            nodal_drift = source_sum == 0 ? 0.0 : target_sum / source_sum - 1,
            mass = geometry.mass,
            mass_from_exact_field = geometry_exact_field.mass,
            mass_from_material_nodes = material_mass,
            area = geometry.area,
            cx = geometry.cx,
            cy = geometry.cy,
            centre_error = hypot(geometry.cx - cexact[1], geometry.cy - cexact[2]),
            area_error = geometry.area / exact_area - 1,
            mass_error = geometry.mass / exact_mass - 1,
            mass_exact_field_error = geometry_exact_field.mass / exact_mass - 1,
            mass_material_error = material_mass / exact_mass - 1,
            constant_state_band = maximum(abs, probe_dest[new_support] .- 1),
            constant_state_core = maximum(abs, probe_dest[core] .- 1),
            qmin_band = minimum(target[new_support]),
            qmax_band = maximum(target[new_support]),
            qmin_inside = minimum(target[inside]),
            qmax_inside = maximum(target[inside]),
            exact_min_inside = minimum(exactv[inside]),
            exact_max_inside = maximum(exactv[inside]),
            error_inf_core = maximum(abs, nodal_error[core]),
            error_l1_core = sum(abs, nodal_error[core]) / count(core),
            error_inf_all = maximum(abs, nodal_error),
            nodes_core = count(core),
            nodes_support = count(new_support),
        ))
        if any(ts -> abs(t - ts) <= dt / 2, snapshot_times) && length(snapshots) < length(snapshot_times)
            push!(snapshots, (; t, triangles = cutTriangles(new_cut),
                              segments = cutSegments(new_cut), phi = copy(phi), rho = copy(rho)))
        end
    end

    return (;
        label = label,
        level = level,
        n = n,
        h = h,
        dt = dt,
        nsteps = nsteps,
        cfl = cfl,
        elapsed = timing.time,
        allocated = timing.bytes,
        gctime = timing.gctime,
        exact_area = exact_area,
        exact_mass = exact_mass,
        mass0_quadrature = geometry0.mass,
        records = records,
        snapshots = snapshots,
        final_rho = copy(rho),
        final_phi = copy(phi),
    )
end

# ---------------------------------------------------------------------------------------
# Reporting
# ---------------------------------------------------------------------------------------

lastrecord(run) = run.records[end]

function reportRun(run, transfer, prolongation)
    r0, rl = run.records[1], lastrecord(run)
    @printf("\n=== %s (level %s): %d² nodes, %d steps, dt = %.6g, CFL = %.2f ===\n",
            run.label, run.level, run.n, run.nsteps, run.dt, run.cfl)
    @printf("time %.3f s, allocated %.1f MiB, gc %.3f s\n",
            run.elapsed, run.allocated / 2^20, run.gctime)
    @printf("\n-- conservation: the two quantities are different --\n")
    @printf("nodal sum  Σq: %.15e → %.15e   worst |drift| %.2e\n",
            r0.source_sum, rl.target_sum,
            maximum(abs, (r.source_sum == 0 ? 0.0 : r.target_sum / r.source_sum - 1)
                          for r in run.records))
    @printf("h²·Σq (uniform-weight 'mass'): %.6f vs physical mass %.6f  → %.1f%% off\n",
            run.h^2 * r0.source_sum, run.mass0_quadrature,
            100 * (run.h^2 * r0.source_sum / run.mass0_quadrature - 1))
    @printf("physical mass M_h: %.6f → %.6f   drift %+.3e (exact %.9f)\n",
            run.mass0_quadrature, rl.mass, rl.mass_error, run.exact_mass)
    @printf("  same quadrature, exact density:  drift %+.3e  (upper bound from cut geometry)\n",
            rl.mass_exact_field_error)
    @printf("  same quadrature, void-side nodes replaced by exact values: drift %+.3e\n",
            rl.mass_material_error)
    @printf("\n-- geometry --\n")
    @printf("area:  %.9f (exact %.9f)  worst |rel err| %.3e\n",
            rl.area, run.exact_area, maximum(abs, r.area_error for r in run.records))
    @printf("centre: (%.6f, %.6f) exact (%.6f, %.6f)  worst error %.3e\n",
            rl.cx, rl.cy, C0[1] + TEND * VEL[1], C0[2] + TEND * VEL[2],
            maximum(r.centre_error for r in run.records))
    @printf("\n-- density --\n")
    @printf("L∞ error on core nodes (φ < -h): %.6e   L¹: %.6e\n",
            rl.error_inf_core, rl.error_l1_core)
    @printf("L∞ error over all nodes:         %.6e\n", rl.error_inf_all)
    @printf("final range on in-phase nodes: [%.4f, %.4f], exact [%.4f, %.4f] → overshoot %+.2e, undershoot %+.2e\n",
            rl.qmin_inside, rl.qmax_inside, rl.exact_min_inside, rl.exact_max_inside,
            rl.qmax_inside - rl.exact_max_inside, rl.exact_min_inside - rl.qmin_inside)
    @printf("final range on the whole transported support (includes halo nodes): [%.4f, %.4f]\n",
            rl.qmin_band, rl.qmax_band)
    @printf("positivity: min over steps, in-phase nodes %.6f; support nodes %.6f\n",
            minimum(r.qmin_inside for r in run.records),
            minimum(r.qmin_band for r in run.records))
    @printf("\n-- constant-state preservation (source ≡ 1 on the old support) --\n")
    @printf("max |target - 1|: band %.6e, core %.6e\n",
            maximum(r.constant_state_band for r in run.records),
            maximum(r.constant_state_core for r in run.records))
    @printf("\n-- where the density error lives at t = %.1f (%s) --\n", TEND, run.label)
    @printf("        φ/h bin        nodes       L∞ err       L¹ err\n")
    for row in layerProfile(run)
        @printf("  [%+5.1f, %+5.1f)  %10d  %12.3e  %12.3e\n",
                row.low, row.high, row.nodes, row.error_inf, row.error_l1)
    end

    @printf("\n-- coarse → fine → finer transfer of the smooth field exp(sin2πx+cos2πy) --\n")
    for s in transfer.steps
        @printf("%s → %s: h %.4f → %.4f, L∞ %.6e, L¹ %.6e\n",
                string(s.from), string(s.to), s.h_from, s.h_to, s.error_inf, s.error_l1)
    end
    @printf("observed order %.2f (second order for a bilinear prolongation of a smooth field)\n",
            transfer.order)
    @printf("fine → coarse resampling: L∞ %.3e (node sampling, exactly nodal on nested grids)\n",
            transfer.restriction_inf)
    @printf("\n-- coarse/fine density comparison at t = %.1f --\n", TEND)
    @printf("prolongated coarse density vs fine density: L∞ %.6e, L¹ %.6e on fine core nodes\n",
            prolongation.error_inf, prolongation.error_l1)
    @printf("\n-- step-by-step trend (%s) --\n", run.label)
    @printf("      t        Σq drift     M_h drift     area err      L∞ core      L¹ core\n")
    for r in run.records[max(1, run.nsteps ÷ 4):max(1, run.nsteps ÷ 4):end]
        @printf("%7.3f  %11.2e  %11.2e  %11.2e  %11.2e  %11.2e\n",
                r.t, r.nodal_drift, r.mass_error, r.area_error,
                r.error_inf_core, r.error_l1_core)
    end
end

"""
    layerProfile(run, bins)

Final density error binned by signed distance in node units, `φ/h`. Shows where the error
lives: the transported band carries a boundary layer a couple of nodes thick whose
partition-of-unity defect (measured separately) decays about one decade per node.
"""
function layerProfile(run; bins = (-6.0, -4.0, -2.0, -1.0, 0.0, 1.0, 2.0))
    n, h = run.n, run.h
    x = collect(range(0.0, 1.0; length = n))
    exact = [rhoExact(a, b, TEND) for a in x, b in x]
    err = run.final_rho .- exact
    phi = run.final_phi
    rows = NamedTuple[]
    for k in 1:(length(bins) - 1)
        mask = (phi .>= bins[k] * h) .& (phi .< bins[k + 1] * h)
        count(mask) == 0 && continue
        push!(rows, (; low = bins[k], high = bins[k + 1], nodes = count(mask),
                     error_inf = maximum(abs, err[mask]),
                     error_l1 = sum(abs, err[mask]) / count(mask)))
    end
    return rows
end

# ---------------------------------------------------------------------------------------
# Entry point
# ---------------------------------------------------------------------------------------

function main()
    outdir = joinpath(@__DIR__, "output")
    mkpath(outdir)

    # 65², 129² and 257² nodes: 64², 128² and 256² cells, exactly nested (factor 2).
    axis = collect(range(0.0, 1.0; length = 65))
    d = Discretization(Dict(:coarse => Grid((axis, copy(axis)))))
    refine!(d, :coarse, :fine; factor = 2)
    refine!(d, :fine, :finer; factor = 2)
    println("levels: ", sort!(collect(keys(d.levels))))
    println("coarsest = ", coarsest(d), ", finest = ", finest(d))
    for id in (:coarse, :fine, :finer)
        println("  ", id, ": ", length(d.levels[id].axes[1]), " nodes/axis, ",
                "immediateCoarser = ", repr(immediateCoarser(d, id)),
                ", immediateFiner = ", repr(immediateFiner(d, id)))
    end

    transfer = transferStudy(d)

    # warm-up: compile the cut/transport path outside the reported timings
    simulate(d, :coarse; cfl = 0.4, tfinal = 0.06, snapshot_times = (), label = "warm-up")

    quick = simulate(d, :coarse; cfl = 0.4, tfinal = TEND,
                     snapshot_times = (0.0, 0.5, 1.0), label = "quick")
    fine = simulate(d, :fine; cfl = 0.4, tfinal = TEND,
                    snapshot_times = (0.0,), label = "fine")

    # per-step error: one step from exact initial data, for the accumulation argument
    onestep = [(level, simulate(d, level; cfl = 0.4, tfinal = quick.dt, snapshot_times = (),
                                label = "one step $level")) for level in (:coarse, :fine)]

    # the timestep is the caller's choice; it must not decide the qualitative answer
    sensitivity = [simulate(d, :coarse; cfl = c, tfinal = TEND, snapshot_times = (),
                            label = "cfl $c") for c in (0.2, 0.4, 0.8)]

    # coarse → fine prolongation of the transported density itself
    prolongated = zeros(fine.n, fine.n)
    density_field = Field(d, Dict(:coarse => quick.final_rho, :fine => prolongated))
    transfer!(density_field, :coarse, :fine)
    core = vec(fine.final_phi) .< -fine.h
    diff = vec(density_field.data[:fine]) .- vec(fine.final_rho)
    prolongation = (;
        error_inf = maximum(abs, diff[core]),
        error_l1 = sum(abs, diff[core]) / count(core),
        values = copy(density_field.data[:fine]),
    )

    reportRun(quick, transfer, prolongation)
    reportRun(fine, transfer, prolongation)

    @printf("\n-- accumulation of the interior error (%d steps in the quick run) --\n", quick.nsteps)
    for (level, run) in onestep
        base = level === :coarse ? quick : fine
        r = lastrecord(run)
        @printf("%s: one-step L¹ core %.3e (L∞ %.3e); after %d steps %.3e → ratio %.1f\n",
                string(level), r.error_l1_core, r.error_inf_core, base.nsteps,
                lastrecord(base).error_l1_core, lastrecord(base).error_l1_core / r.error_l1_core)
    end

    @printf("\n-- timestep sensitivity at %d² nodes (the nodal sum is an invariant, the mass is not) --\n",
            quick.n)
    @printf("  CFL   steps      Σq drift     M_h drift      L¹ core      area err\n")
    for run in sensitivity
        r = lastrecord(run)
        @printf("%5.1f  %6d  %11.2e  %11.2e  %11.2e  %11.2e\n", run.cfl, run.nsteps,
                maximum(abs, (x.source_sum == 0 ? 0.0 : x.target_sum / x.source_sum - 1)
                        for x in run.records),
                r.mass_error, r.error_l1_core, r.area_error)
    end

    results = (; problem = (; C0, RIN, ROUT, VEL, TEND, BUMP, WIDTH),
               transfer, prolongation, quick, fine, sensitivity,
               onestep = Dict(level => run for (level, run) in onestep))
    path = joinpath(outdir, "results.jls")
    serialize(path, results)
    @printf("\nwrote %s\n", path)
    @printf("run the figures with: julia %s\n",
            joinpath(@__DIR__, "visualize.ConservativeAnnulus.jl"))
    return results
end

if abspath(PROGRAM_FILE) == @__FILE__
    main()
end
