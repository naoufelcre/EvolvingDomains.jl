#!/usr/bin/env julia
#
# NutrientDrivenGrowth
# ====================
# Two-dimensional, free-boundary, nutrient-driven tissue growth on a Cartesian background
# box D = [0,1]^2, orchestrated by EvolvingDomains and discretised by Gridap +
# GridapEmbedded (CutFEM: Nitsche Dirichlet + ghost-penalty stabilisation).
#
# Mathematical model (all fields live on the current physical domain Omega(t)):
#
#   nutrient   -Dc Δc + λ c = 0      in Omega,      c = 1        on Gamma = ∂Omega
#   growth      g(c) = g0 max(c-c*,0)
#   mechanics  -Δp = g(c)            in Omega,      p = 0        on Gamma   (Darcy: v = -∇p)
#   interface  V_n = v⋅n             <=>           ∂_t φ + v⋅∇φ = 0            (φ<0 inside)
#   marker     ∂_t q + v⋅∇q = 0      in Omega       (passive material field)
#
# Data flow per time step (see `dynamics` and `simulate`), all through public interfaces:
#
#   EvolvingDomains Field(φ) --> Cut snapshot (diagnostics: physical area from sub-cell
#        triangles, cell centroids, connected components, aggregateCutCells probe).
#   φ (nodal Array, first axis fastest) --> Gridap Q1 CellField on the same background model
#        --> DiscreteGeometry --> GE.cut --> PHYSICAL_IN / ACTIVE_IN triangulations,
#        EmbeddedBoundary, GhostSkeleton, measures.
#   solveNutrient --> g(c) as an FE field --> solvePressure --> ∇p sampled per active cell
#        --> nodal velocity on the background grid --> eulerianTransport!(φ) (WENO5+SSP-RK3)
#        --> reinitialize!(φ) --> lagrangianTransport!(q) (CIP, separate source/target buffers).
#
# Verified public-API facts (checked at runtime, never assumed silently):
#   * Gridap's CartesianDiscreteModel node order equals `vec` of an (nx, ny) node array, i.e.
#     first axis fastest (`checkNodeOrdering`).
#   * EvolvingDomains Cut 2D subcells are triangles (N+1 vertices) and `subcell_states` marks
#     the physical phase, so the physical part of each cell is available from public arrays.
#   * The active-cell FE space has one dof per node of an active cell, in the node order of
#     `get_node_coordinates(Ωact)`; `backgroundNodal` reuses that correspondence.
#
# EvolvingDomains identity is *enforced*, not assumed: `checkEvolvingDomains` refuses to run
# unless the loaded package is version 1.0 of this repository (the default environment stacked
# under `validation/env` pins a different pre-1.0 checkout) and records the root, version and
# source hashes in run_metadata.txt.
#
# Usage:  julia --project=../env --startup-file=no run.NutrientDrivenGrowth.jl \
#             [ncells] [tfinal] [output-name] [dtmax]
# Output: ./output/<tag>/ (metrics.csv, fields.jls, run_metadata.txt)

import EvolvingDomains as ED
import GridapEmbedded as GE

using Gridap
using Gridap.Geometry: get_node_coordinates
using GridapEmbedded.LevelSetCutters: DiscreteGeometry
using StaticArrays: SVector
using LinearAlgebra: norm
using Printf: @printf, @sprintf
using Serialization: serialize
using SHA: sha256

const MODEL_ORDER = 1     # P1/Q1 spaces: DiscreteGeometry requires one dof per node.

# `normpath` keeps a trailing separator when it resolves to a directory boundary; comparisons
# between package roots must not depend on that.
_cleanpath(path) = rstrip(normpath(path), ('/', '\\'))

# The repository this case belongs to; the run refuses any other EvolvingDomains.
const REPOSITORY_ROOT = _cleanpath(joinpath(@__DIR__, "..", ".."))

"""
    checkEvolvingDomains() -> Vector{String}

Hard guard on the *identity* of the EvolvingDomains that is loaded, because the Julia load path
here is `validation/env` **stacked on** the user default environment
`~/.julia/environments/v1.11`, which pins a different EvolvingDomains
(`path = "…/EvolvingDomains.jl"`, version 0.5.0, the pre-redesign checkout).  The path source
`[sources] EvolvingDomains = {path = "../.."}` in `validation/env` wins over that copy; this
function proves it instead of assuming it:

* package root, loaded source file and version must be **this repository's** 1.0,
* the v1.0 public surface must be present,
* names that only exist in the pre-1.0 checkout must be absent.

Returns the provenance lines that go into `run_metadata.txt`.
"""
function checkEvolvingDomains()
    root = _cleanpath(pkgdir(ED))
    version = pkgversion(ED)
    source = _cleanpath(pathof(ED))
    root == REPOSITORY_ROOT || error(
        "EvolvingDomains was loaded from \"$root\" but this case requires the current " *
        "repository \"$REPOSITORY_ROOT\". Run with --project=" *
        joinpath(REPOSITORY_ROOT, "validation", "env"))
    version >= v"1.0.0" || error(
        "EvolvingDomains $version was loaded; this case requires version 1.0 from the " *
        "current repository (the 0.5.0 checkout next to it exports a different API)")
    startswith(source, root * "/") || error(
        "pathof(EvolvingDomains) = $source is outside the repository root $root")

    for name in (:Grid, :Discretization, :Field, :cut, :Cut, :CutSurface, :CUT_IN, :CUT_OUT,
        :CUT_CUT, :aggregateCutCells, :reinitialize!, :eulerianTransport!,
        :lagrangianTransport!, :WENO5Cache, :CIPCache, :TransportMap)
        isdefined(ED, name) || error("EvolvingDomains 1.0 must define $name, it does not")
    end
    for name in (:EvolvingDiscreteGeometry, :CartesianGridInfo,
        :AbstractVelocitySource, :aggregate_cut_cells)
        isdefined(ED, name) && error(
            "EvolvingDomains defines $name, which belongs to the pre-1.0 API: this is not " *
            "the redesigned package")
    end

    lines = String["evolvingdomains_root = $root",
        "evolvingdomains_version = $version",
        "evolvingdomains_source = $source"]
    for file in ("src/EvolvingDomains.jl", "src/DataModel.jl", "src/CutData.jl",
        "src/Kinematic/Kinematic.jl", "src/Geometric/Geometric.jl")
        path = joinpath(root, file)
        push!(lines, "evolvingdomains_sha256[$file] = $(bytes2hex(sha256(read(path))))")
    end
    push!(lines, "load_path = $(join(Base.load_path(), " : "))")
    return lines
end

# =============================================================================
# 1. Configuration
# =============================================================================

"""
    defaultConfig(; kwargs...) -> NamedTuple

Every physical, numerical and bookkeeping parameter of the experiment. Keyword defaults are
the baseline configuration of the validation protocol.
"""
function defaultConfig(; ncells::Int = 64,
        tfinal::Float64 = 0.5,
        dtmax::Float64 = 5.0e-3,
        cfl::Float64 = 0.4,
        Dc::Float64 = 0.02,
        lambda::Float64 = 0.05,
        g0::Float64 = 1.5,
        cstar::Float64 = 0.2,
        R0::Float64 = 0.15,
        center1::Tuple{Float64,Float64} = (0.33, 0.50),
        center2::Tuple{Float64,Float64} = (0.67, 0.50),
        nitsche::Float64 = 10.0,
        ghost::Float64 = 0.1,
        qdegree::Int = 3,
        outdir::String = joinpath(@__DIR__, "output"))
    ncells >= 8 || throw(ArgumentError("need at least 8 cells per axis"))
    tfinal > 0 || throw(ArgumentError("tfinal must be positive"))
    (; ncells, tfinal, dtmax, cfl, Dc, lambda, g0, cstar, R0, center1, center2,
        nitsche, ghost, qdegree, outdir)
end

# =============================================================================
# 2. Background grid and initial data
# =============================================================================

"""
    backgroundModel(cfg) -> NamedTuple

The fixed Cartesian backdrop: a Gridap `CartesianDiscreteModel` with `ncells^2` cells and
`(ncells+1)^2` nodes, the nodal Q1 space used only to hand the level set to GridapEmbedded,
and the axis vector shared with EvolvingDomains.
"""
function backgroundModel(cfg)
    nx = cfg.ncells + 1
    axes = collect(range(0.0, 1.0; length = nx))
    model = CartesianDiscreteModel((0.0, 1.0, 0.0, 1.0), (cfg.ncells, cfg.ncells))
    Vnodal = FESpace(Triangulation(model), ReferenceFE(lagrangian, Float64, 1))
    (; model, axes, Vnodal, nx, h = 1.0 / cfg.ncells, ncx = cfg.ncells)
end

"""
    checkNodeOrdering(backdrop) -> Float64

Verify that `vec` of an `(nx, nx)` EvolvingDomains node array (first axis fastest) enumerates
the same points in the same order as Gridap's model nodes. Throws if the two representations
of the grid ever disagree; the conversion would then have to be an explicit permutation.
"""
function checkNodeOrdering(backdrop)
    gridap = [Tuple(p) for p in collect(get_node_coordinates(backdrop.model))]
    evolving = [(a, b) for b in backdrop.axes for a in backdrop.axes]
    length(gridap) == length(evolving) || error(
        "node count mismatch: Gridap $(length(gridap)) vs EvolvingDomains $(length(evolving))")
    worst = 0.0
    for k in eachindex(gridap)
        worst = max(worst, abs(gridap[k][1] - evolving[k][1]),
            abs(gridap[k][2] - evolving[k][2]))
    end
    worst <= 100eps(Float64) || error(
        "Gridap and EvolvingDomains node orders differ (max deviation $worst); " *
        "an explicit permutation would be required")
    return worst
end

"""Negative-inside level set of two discs: min of the two signed distances."""
function initialLevelSet(axes, cfg)
    (c1, c2, R0) = (cfg.center1, cfg.center2, cfg.R0)
    return [min(hypot(a - c1[1], b - c1[2]) - R0, hypot(a - c2[1], b - c2[2]) - R0)
            for a in axes, b in axes]
end

"""Passive material marker: 0.8 near the left disc, 0.2 near the right one."""
function initialMarker(axes, cfg)
    (c1, c2) = (cfg.center1, cfg.center2)
    return [0.8exp(-((a - c1[1])^2 + (b - c1[2])^2) / 0.06^2) +
            0.2exp(-((a - c2[1])^2 + (b - c2[2])^2) / 0.06^2) for a in axes, b in axes]
end

# =============================================================================
# 3. EvolvingDomains cut snapshot -> geometric quantities
# =============================================================================

"""
    cutGeometryInfo(snapshot, grid, h) -> NamedTuple

Physical-part geometry of every background cell, read from the public `Cut` arrays only.
In 2D each subcell is a triangle; `subcell_states == CUT_IN` marks the part of a cell that
belongs to the physical phase, so `inside_area[c]` is the exact area of `Omega ∩ cell` for the
piecewise-linear cut model and `cx[c], cy[c]` its centroid. `min_fraction` is the smallest
non-empty cell area in units of `h^2` (cut-cell conditioning indicator); `aggregate_roots` is
the number of distinct roots returned by the public `aggregateCutCells` on this topology.
"""
function cutGeometryInfo(snapshot::ED.Cut, grid::ED.Grid, h)
    ncells = length(snapshot.bg_states)
    inside_area = zeros(ncells)
    cx = zeros(ncells)
    cy = zeros(ncells)
    pointers = snapshot.point_pointers
    ids = snapshot.point_ids
    points = snapshot.physical_points
    for s in eachindex(snapshot.subcell_states)
        snapshot.subcell_states[s] == ED.CUT_IN || continue
        lo = Int(pointers[s])
        hi = Int(pointers[s + 1]) - 1
        hi - lo + 1 == 3 || error("expected triangular 2D subcells, got $(hi - lo + 1) vertices")
        a = points[Int(ids[lo])]
        b = points[Int(ids[lo + 1])]
        c = points[Int(ids[lo + 2])]
        area = abs((b[1] - a[1]) * (c[2] - a[2]) - (b[2] - a[2]) * (c[1] - a[1])) / 2
        cell = Int(snapshot.cell_to_bgcell[s])
        inside_area[cell] += area
        cx[cell] += area * (a[1] + b[1] + c[1]) / 3
        cy[cell] += area * (a[2] + b[2] + c[2]) / 3
    end
    for cell in eachindex(inside_area)
        inside_area[cell] > 0 || continue
        cx[cell] /= inside_area[cell]
        cy[cell] /= inside_area[cell]
    end
    occupied = findall(>(0.0), inside_area)
    min_fraction = isempty(occupied) ? 0.0 : minimum(inside_area[occupied]) / h^2
    roots = ED.aggregateCutCells(snapshot, grid)
    (; inside_area, cx, cy, occupied, min_fraction, area = sum(inside_area),
        n_cut = count(==(ED.CUT_CUT), snapshot.bg_states),
        n_in = count(==(ED.CUT_IN), snapshot.bg_states),
        aggregate_roots = length(unique(r for r in roots if r > 0)))
end

"""
    connectedComponents(active, n) -> Int

Number of 4-connected components of a boolean mask over an `n × n` lattice, in the node order
(first axis fastest) used everywhere in this file.
"""
function connectedComponents(active::AbstractVector{Bool}, n::Int)
    labels = zeros(Int, n * n)
    stack = Int[]
    components = 0
    for seed in 1:(n * n)
        (active[seed] && labels[seed] == 0) || continue
        components += 1
        labels[seed] = components
        push!(stack, seed)
        while !isempty(stack)
            item = pop!(stack)
            i = (item - 1) % n + 1
            j = (item - 1) ÷ n + 1
            for (di, dj) in ((1, 0), (-1, 0), (0, 1), (0, -1))
                ii = i + di
                jj = j + dj
                (1 <= ii <= n && 1 <= jj <= n) || continue
                other = ii + (jj - 1) * n
                (active[other] && labels[other] == 0) || continue
                labels[other] = components
                push!(stack, other)
            end
        end
    end
    return components
end

"""
    componentCounts(snapshot, levelset, backdrop) -> NamedTuple

Two resolution-dependent topology diagnostics, both with 4-connectivity and both measured on
this Cartesian grid:

  * `nodes`: components of the set of grid nodes where `φ < 0`;
  * `cells`: components of the set of active background cells (`bg_states != CUT_OUT`), i.e.
    cells that are fully inside or crossed by the interface.

The nodal count flips when a node first becomes negative, the cell count one cell earlier
(a cut band already connects the two regions while no node is negative yet). Neither is a
sub-cell-exact topology test.
"""
componentCounts(snapshot, levelset, backdrop) = (;
    nodes = connectedComponents(vec(levelset .< 0.0), backdrop.nx),
    cells = connectedComponents(snapshot.bg_states .!= ED.CUT_OUT, backdrop.ncx))

# =============================================================================
# 4. Embedded geometry handed to GridapEmbedded
# =============================================================================

"""
    buildGeometry(levelset, backdrop, cfg) -> NamedTuple

Convert the current EvolvingDomains level-set array into the embedded discretisation used by
both elliptic solvers. Documented conversion: the nodal array (first axis fastest) becomes the
free-dof vector of the Q1 background space, which `DiscreteGeometry` reads as a nodal level
set on the background model; `GE.cut` then produces the sub-cell geometry.
"""
function buildGeometry(levelset, backdrop, cfg)
    values = vec(levelset)
    if any(v -> abs(v) <= 1.0e-12, values)          # a dof exactly on the interface is degenerate
        values = copy(values)
        values[findall(v -> abs(v) <= 1.0e-12, values)] .+= 100eps(Float64)
    end
    φh = FEFunction(backdrop.Vnodal, values)
    geo = DiscreteGeometry(φh, backdrop.model)
    cutgeo = GE.cut(backdrop.model, geo)
    Ω = Triangulation(cutgeo, GE.PHYSICAL_IN)         # physical part (sub-cell triangles)
    Ωact = Triangulation(cutgeo, GE.ACTIVE_IN)        # active cells (inside + cut)
    Γ = GE.EmbeddedBoundary(cutgeo)                   # material boundary ∂Omega
    Γg = GE.GhostSkeleton(cutgeo, GE.ACTIVE_IN)       # ghost-penalty skeleton
    dΩ = Measure(Ω, cfg.qdegree)
    dΓ = Measure(Γ, cfg.qdegree)
    dΓg = Measure(Γg, cfg.qdegree)
    (; cutgeo, Ω, Ωact, Γ, Γg, dΩ, dΓ, dΓg,
        nΓ = get_normal_vector(Γ), nΓg = get_normal_vector(Γg))
end

"""Conforming P1 test/trial spaces on the active cells (one dof per node of an active cell)."""
function buildSpaces(geo)
    reffe = ReferenceFE(lagrangian, Float64, MODEL_ORDER)
    V = TestFESpace(geo.Ωact, reffe, conformity = :H1)
    return (; V, U = TrialFESpace(V))
end

"""Physical area of the current domain, `∫_Ω 1 dΩ`."""
domainArea(geo) = sum(∫(1.0) * geo.dΩ)

"""Length of the material boundary, `∫_Γ 1 dΓ`."""
interfaceLength(geo) = sum(∫(1.0) * geo.dΓ)

# =============================================================================
# 5. Elliptic solves (CutFEM)
# =============================================================================

"""
    solveElliptic(geo, spaces, cfg, h; diffusion, reaction, source, boundaryValue, label)

Assemble and solve the CutFEM Dirichlet problem

    -∇⋅(diffusion ∇u) + reaction u = source   in Ω,      u = boundaryValue   on Γ

with symmetric Nitsche imposition of the Dirichlet data (penalty `diffusion*cfg.nitsche/h`)
and the facet ghost penalty of Burman & Hansbo (`diffusion*cfg.ghost*h`). Both penalty terms
carry the diffusion coefficient so that the stabilisation is independent of the physical
coefficient of the operator. Returns the FE function, its dofs, the algebraic residual
`‖b - A u‖/‖b‖` and the assembly/solve timings.
"""
function solveElliptic(geo, spaces, cfg, h; diffusion, reaction, source, boundaryValue, label)
    γd = cfg.nitsche / h
    γg = cfg.ghost * h
    nΓ = geo.nΓ
    nΓg = geo.nΓg
    V = spaces.V
    U = spaces.U

    a(u, v) = ∫(diffusion * (∇(v) ⋅ ∇(u)) + reaction * v * u) * geo.dΩ +
              ∫(diffusion * (γd * v * u - v * (nΓ ⋅ ∇(u)) - (nΓ ⋅ ∇(v)) * u)) * geo.dΓ +
              ∫(diffusion * γg * jump(nΓg ⋅ ∇(v)) * jump(nΓg ⋅ ∇(u))) * geo.dΓg
    if source === nothing
        l = v -> ∫(diffusion * (γd * v - (nΓ ⋅ ∇(v))) * boundaryValue) * geo.dΓ
    else
        l = v -> ∫(v * source) * geo.dΩ +
                  ∫(diffusion * (γd * v - (nΓ ⋅ ∇(v))) * boundaryValue) * geo.dΓ
    end

    t0 = time()
    op = AffineFEOperator(a, l, U, V)
    matrix = get_matrix(op)
    rhs = get_vector(op)
    t_assembled = time()
    uh = solve(op)
    t_solved = time()

    dofs = get_free_dof_values(uh)
    all(isfinite, dofs) || error("$label: solver returned non-finite values")
    residual = norm(rhs - matrix * dofs) / max(norm(rhs), eps(Float64))
    (; uh, dofs, residual, t_assembly = t_assembled - t0, t_solve = t_solved - t_assembled)
end

"""
    solveNutrient(geo, spaces, cfg, h) -> NamedTuple

Quasi-steady nutrient problem `-Dc Δc + λc = 0` in Ω with `c = 1` on Γ. Besides the FE
solution it returns the boundary flux `∫_Γ n⋅∇c`, the mass `∫_Ω c`, the Dirichlet defect
`∫_Γ (c-1)` (the Nitsche boundary error) and the concentration extrema.
"""
function solveNutrient(geo, spaces, cfg, h)
    state = solveElliptic(geo, spaces, cfg, h; diffusion = cfg.Dc, reaction = cfg.lambda,
        source = nothing, boundaryValue = 1.0, label = "nutrient")
    flux = sum(∫(geo.nΓ ⋅ ∇(state.uh)) * geo.dΓ)
    (; state..., cmin = minimum(state.dofs), cmax = maximum(state.dofs),
        flux, mass = sum(∫(state.uh) * geo.dΩ),
        boundary_defect = sum(∫(state.uh - 1.0) * geo.dΓ))
end

"""
    computeGrowthRate(nutrient, geo, spaces, cfg) -> NamedTuple

Volumetric growth rate `g(c) = g0 max(c - c*, 0)` as an FE field on the same space, its dof
vector, and its integral over the current physical domain `∫_Ω g(c) dΩ`.
"""
function computeGrowthRate(nutrient, geo, spaces, cfg)
    values = cfg.g0 .* max.(nutrient.dofs .- cfg.cstar, 0.0)
    gh = FEFunction(spaces.V, values)
    (; gh, values, integral = sum(∫(gh) * geo.dΩ))
end

"""
    solvePressure(geo, growth, spaces, cfg, h) -> NamedTuple

Pressure Poisson problem `-Δp = g(c)` in Ω with `p = 0` on Γ, plus the boundary flux
`∫_Γ n⋅∇p` and the assembly-consistency check `-∫_Γ n⋅∇p + (γ/h)∫_Γ p = ∫_Ω g` obtained by
testing the discrete weak form with the constant function.
"""
function solvePressure(geo, growth, spaces, cfg, h)
    state = solveElliptic(geo, spaces, cfg, h; diffusion = 1.0, reaction = 0.0,
        source = growth.gh, boundaryValue = 0.0, label = "pressure")
    flux = sum(∫(geo.nΓ ⋅ ∇(state.uh)) * geo.dΓ)
    flux_identity = -flux + (cfg.nitsche / h) * sum(∫(state.uh) * geo.dΓ)
    (; state..., flux, flux_identity, pmin = minimum(state.dofs),
        pmax = maximum(state.dofs))
end

# =============================================================================
# 6. Velocity recovery and extension to the background grid
# =============================================================================

"""
    backgroundNodal(dofs, geo, backdrop) -> Matrix

Place the dofs of the active-cell space (one per node of an active cell, in the space's node
order) on the background `(nx, nx)` grid; inactive nodes become `NaN`. The dof/node
correspondence is taken from `get_node_coordinates(Ωact)` and is checked to be a bijection.
"""
function backgroundNodal(dofs, geo, backdrop)
    nx = backdrop.nx
    h = backdrop.h
    coordinates = collect(get_node_coordinates(geo.Ωact))
    length(coordinates) == length(dofs) || error(
        "active-space dof count $(length(dofs)) differs from active node count $(length(coordinates))")
    out = fill(NaN, nx, nx)
    for (k, point) in enumerate(coordinates)
        i = round(Int, point[1] / h) + 1
        j = round(Int, point[2] / h) + 1
        (1 <= i <= nx && 1 <= j <= nx) || error("active node outside the background grid")
        isfinite(out[i, j]) && error("two active dofs map to node ($i, $j)")
        out[i, j] = dofs[k]
    end
    return out
end

"""
    velocityFromPressure(pressure, geo, cutinfo, snapshot, backdrop) -> NamedTuple

Recover `v = -∇p` on the background grid. Two different objects are produced, because they
serve two different purposes.

*Cell-sampled velocity* `cell_vx, cell_vy` (used for the discrete divergence diagnostic and the
plots): one sample per active background cell, the bilinear derivative of the Q1 pressure field
evaluated at the area-weighted centroid of the cell's physical part (inside Ω, as close to the
interface as the cell allows).

*Nodal transport velocity* `velocity` (consumed by `eulerianTransport!` and CIP): a zero-order
normal velocity extension in the sense of the **closest interface point**. Every interface facet
of the EvolvingDomains cut snapshot supplies one sample: the bilinear derivative of `p`
evaluated at the facet midpoint in the facet's own background cell, i.e. `-∇p` *at* the
interface. Each background node then takes the velocity of the nearest facet midpoint, which
gives every node on either side of a given interface the velocity of that interface.

Known limitation, documented rather than hidden: the transport velocity is a *single-valued*
field, so two opposing interfaces closer than a cell cannot both be represented. The nearest
facet wins and the topology diagnostics (`components_nodes`, `components_cells`) and the
area-growth residual show the consequences.
"""
function velocityFromPressure(pressure, geo, cutinfo, snapshot, backdrop)
    nx = backdrop.nx
    ncx = backdrop.ncx
    h = backdrop.h
    pn = backgroundNodal(pressure.dofs, geo, backdrop)

    # 1. Cell samples at the centroid of each cell's physical part.
    cell_vx = zeros(ncx, ncx)
    cell_vy = zeros(ncx, ncx)
    for cell in cutinfo.occupied
        i = (cell - 1) % ncx + 1
        j = (cell - 1) ÷ ncx + 1
        gradient = _bilinearGradient(pn, i, j, cutinfo.cx[cell], cutinfo.cy[cell], h)
        gradient === nothing && continue
        cell_vx[i, j] = -gradient[1]
        cell_vy[i, j] = -gradient[2]
    end

    # 2. Interface samples: -∇p at each material-boundary facet midpoint.
    surface = snapshot.boundary
    nfacets = length(surface.facet_to_bgcell)
    nfacets > 0 || error("the cut snapshot has no interface facet; cannot extend the velocity")
    midpoint_x = zeros(nfacets)
    midpoint_y = zeros(nfacets)
    facet_vx = zeros(nfacets)
    facet_vy = zeros(nfacets)
    interface_flux = 0.0
    facet_speed_max = 0.0
    for facet in 1:nfacets
        lo = Int(surface.point_pointers[facet])
        p1 = surface.physical_points[Int(surface.point_ids[lo])]
        p2 = surface.physical_points[Int(surface.point_ids[lo + 1])]
        mx = (p1[1] + p2[1]) / 2
        my = (p1[2] + p2[2]) / 2
        cell = Int(surface.facet_to_bgcell[facet])
        i = (cell - 1) % ncx + 1
        j = (cell - 1) ÷ ncx + 1
        midpoint_x[facet] = mx
        midpoint_y[facet] = my
        # The facet inherits the velocity of its own cell, sampled at the centroid of that
        # cell's physical part: the Q1 gradient at the facet midpoint itself is visibly more
        # sensitive to the ghost nodes of cells the interface only clips.
        facet_vx[facet] = cell_vx[i, j]
        facet_vy[facet] = cell_vy[i, j]
        normal = surface.normals[facet]
        interface_flux += (facet_vx[facet] * normal[1] + facet_vy[facet] * normal[2]) *
                          hypot(p2[1] - p1[1], p2[2] - p1[2])
        facet_speed_max = max(facet_speed_max, hypot(facet_vx[facet], facet_vy[facet]))
    end
    all(isfinite, facet_vx) || error("non-finite interface velocity sample")

    # 3. Nearest-facet extension of the interface velocity to every background node.
    vx = zeros(nx, nx)
    vy = zeros(nx, nx)
    for j in 1:nx, i in 1:nx
        x = (i - 1) * h
        y = (j - 1) * h
        best = 0
        best_distance = Inf
        for facet in 1:nfacets
            distance = (x - midpoint_x[facet])^2 + (y - midpoint_y[facet])^2
            if distance < best_distance
                best_distance = distance
                best = facet
            end
        end
        vx[i, j] = facet_vx[best]
        vy[i, j] = facet_vy[best]
    end

    velocity = SVector{2,Float64}[SVector{2,Float64}(vx[k], vy[k]) for k in eachindex(vx)]
    (; velocity, vx, vy, cell_vx, cell_vy, interface_flux, facet_speed_max,
        maximum_speed = maximum(norm, velocity))
end

"""
    _bilinearGradient(nodal, i, j, x, y, h) -> SVector or nothing

Gradient of the Q1 interpolant of `nodal` on background cell `(i, j)` at the physical point
`(x, y)`, using the standard bilinear derivative formula; `nothing` if any of the four nodal
values is undefined (a node outside the active region).
"""
function _bilinearGradient(nodal, i, j, x, y, h)
    nx = size(nodal, 1)
    (1 <= i <= nx - 1 && 1 <= j <= nx - 1) || return nothing
    u00 = nodal[i, j]
    u10 = nodal[i + 1, j]
    u01 = nodal[i, j + 1]
    u11 = nodal[i + 1, j + 1]
    all(isfinite, (u00, u10, u01, u11)) || return nothing
    ξ = clamp((x - (i - 1) * h) / h, 0.0, 1.0)
    η = clamp((y - (j - 1) * h) / h, 0.0, 1.0)
    return SVector(((u10 - u00) * (1 - η) + (u11 - u01) * η) / h,
        ((u01 - u00) * (1 - ξ) + (u11 - u10) * ξ) / h)
end

# =============================================================================
# 7. Diagnostics
# =============================================================================

"""
    divergenceMismatch(cell_vx, cell_vy, growth_nodes, snapshot, backdrop) -> NamedTuple

Mismatch between the discrete divergence of the recovered velocity and `g(c)`.

Discretisation, stated explicitly: the divergence is taken with central differences on the
*cell-sampled* velocity field created by `velocityFromPressure` (one sample per active cell,
`(vx[i+1,j]-vx[i-1,j])/(2h) + (vy[i,j+1]-vy[i,j-1])/(2h)`), and `g(c)` is bilinear-interpolated to
the same cell centres. Only cells whose full 3×3 neighbourhood is `CUT_IN` are used, so no
interface cell, ghost node or extended velocity enters the norm. The samples are Q1 derivatives
of the pressure field, so the mismatch measures the accuracy of the discrete Darcy law, not an
algebraic residual of the pressure solve.
"""
function divergenceMismatch(cell_vx, cell_vy, growth_nodes, snapshot, backdrop)
    ncx = backdrop.ncx
    h = backdrop.h
    states = snapshot.bg_states
    inside(i, j) = 1 <= i <= ncx && 1 <= j <= ncx && states[i + (j - 1) * ncx] == ED.CUT_IN
    num = 0.0
    den = 0.0
    worst = 0.0
    worst_scale = 0.0
    used = 0
    for j in 1:ncx, i in 1:ncx
        (inside(i - 1, j - 1) && inside(i, j - 1) && inside(i + 1, j - 1) &&
         inside(i - 1, j) && inside(i, j) && inside(i + 1, j) &&
         inside(i - 1, j + 1) && inside(i, j + 1) && inside(i + 1, j + 1)) || continue
        g = (growth_nodes[i, j] + growth_nodes[i + 1, j] +
             growth_nodes[i, j + 1] + growth_nodes[i + 1, j + 1]) / 4
        isfinite(g) || continue
        div = (cell_vx[i + 1, j] - cell_vx[i - 1, j]) / (2h) +
              (cell_vy[i, j + 1] - cell_vy[i, j - 1]) / (2h)
        isfinite(div) || continue
        num += (div - g)^2
        den += g^2
        used += 1
        if abs(div - g) > worst
            worst = abs(div - g)
            worst_scale = abs(g)
        end
    end
    used == 0 && return (; l2 = NaN, linf = NaN, nodes = 0)
    (; l2 = sqrt(num / max(den, eps(Float64))),
        linf = worst / max(worst_scale, eps(Float64)), nodes = used)
end

# =============================================================================
# 8. Geometry and marker advance
# =============================================================================

"""
    advanceGeometry!(levelset, velocity, dt, cfg, grid, h) -> (Cut, Float64)

Advance the level set with the public `eulerianTransport!` (WENO5 + SSP-RK3, frozen nodal
velocity) and re-distance it with `reinitialize!` at the new interface. Returns the fresh cut
snapshot and the physical area right after the transport, so that the area change caused by
re-distancing is measured rather than assumed away.
"""
function advanceGeometry!(levelset, velocity, dt, cfg, grid, h)
    ED.eulerianTransport!(levelset, velocity, dt; level_id = :init, cache = cfg.weno_cache)
    area_transported = cutGeometryInfo(ED.cut(levelset, :init), grid, h).area
    snapshot = ED.cut(levelset, :init)
    ED.reinitialize!(levelset, snapshot)
    return snapshot, area_transported
end

"""Transport the passive marker with intensive CIP; source and target buffers must not alias."""
function advanceMarker!(target, source, velocity, dt, cfg)
    ED.lagrangianTransport!(target, source, velocity, dt; level_id = :init, cache = cfg.cip_cache)
    return target, source
end

"""Explicit geometry CFL restriction from the computed nodal velocity."""
function timestep(velocity, backdrop, cfg, t)
    vmax = maximum(norm, velocity)
    cfl_dt = vmax > 0 ? cfg.cfl * backdrop.h / vmax : cfg.dtmax
    return min(cfl_dt, cfg.dtmax, cfg.tfinal - t), vmax
end

# =============================================================================
# 9. Independent verification against analytic solutions
# =============================================================================

"""Modified Bessel function I0 by its power series (avoids an extra dependency)."""
function besselI0(z)
    term = 1.0
    total = 1.0
    for m in 1:40
        term *= (z / 2)^2 / (m * m)
        total += term
        abs(term) < 1.0e-18 * abs(total) && break
    end
    return total
end

"""
    analyticVerification(backdrop, cfg) -> String

Independent check of both elliptic solvers on the *initial* geometry, where closed-form
solutions exist on the two discs of radius R0:

    pressure   -Δp = 1,             p = 0 on each disc  =>  p = (R0² - r²)/4
    nutrient   -Dc Δc + λc = 0,     c = 1 on each disc  =>  c = I0(k r)/I0(k R0), k=√(λ/Dc)

The domain is measured against `2πR0²` (two discs), and both errors are measured against the
analytic fields on the *discrete* (polygonal) domain: nodal maximum over nodes of fully-inside
cells, and the relative L2 error against the interpolant of the analytic solution.
"""
function analyticVerification(backdrop, cfg)
    grid = ED.Grid((backdrop.axes, backdrop.axes))
    levelset = ED.Field(ED.Discretization(Dict(:init => grid)),
        Dict(:init => initialLevelSet(backdrop.axes, cfg)))
    snapshot = ED.cut(levelset, :init)
    cutinfo = cutGeometryInfo(snapshot, grid, backdrop.h)
    geo = buildGeometry(levelset.data[:init], backdrop, cfg)
    spaces = buildSpaces(geo)

    R0 = cfg.R0
    c1 = cfg.center1
    c2 = cfg.center2
    radius(x) = min(hypot(x[1] - c1[1], x[2] - c1[2]),
        hypot(x[1] - c2[1], x[2] - c2[2]))
    pfun(x) = (R0^2 - radius(x)^2) / 4
    kappa = sqrt(cfg.lambda / cfg.Dc)
    cfun(x) = besselI0(kappa * radius(x)) / besselI0(kappa * R0)

    inside_nodes = falses(backdrop.nx, backdrop.nx)
    for j in 1:backdrop.ncx, i in 1:backdrop.ncx
        snapshot.bg_states[i + (j - 1) * backdrop.ncx] == ED.CUT_IN || continue
        for (ii, jj) in ((i, j), (i + 1, j), (i, j + 1), (i + 1, j + 1))
            inside_nodes[ii, jj] = true
        end
    end

    p_exact = [pfun((a, b)) for a in backdrop.axes, b in backdrop.axes]
    c_exact = [cfun((a, b)) for a in backdrop.axes, b in backdrop.axes]

    pressure = solveElliptic(geo, spaces, cfg, backdrop.h; diffusion = 1.0, reaction = 0.0,
        source = 1.0, boundaryValue = 0.0, label = "verify-pressure")
    nutrient = solveElliptic(geo, spaces, cfg, backdrop.h; diffusion = cfg.Dc,
        reaction = cfg.lambda, source = nothing, boundaryValue = 1.0, label = "verify-nutrient")
    p_nodes = backgroundNodal(pressure.dofs, geo, backdrop)
    c_nodes = backgroundNodal(nutrient.dofs, geo, backdrop)

    p_nodal = maximum(abs(p_nodes[i, j] - p_exact[i, j]) for i in 1:backdrop.nx,
        j in 1:backdrop.nx if inside_nodes[i, j])
    c_nodal = maximum(abs(c_nodes[i, j] - c_exact[i, j]) for i in 1:backdrop.nx,
        j in 1:backdrop.nx if inside_nodes[i, j])
    p_interp = interpolate(pfun, spaces.V)
    c_interp = interpolate(cfun, spaces.V)
    l2(value, reference) = sqrt(sum(∫((value - reference) * (value - reference)) * geo.dΩ) /
                                sum(∫(reference * reference) * geo.dΩ))
    exact_area = 2π * R0^2
    return join((
        "verification_area_relative_error = $(abs(cutinfo.area - exact_area) / exact_area)",
        "verification_pressure_nodal_max_error = $p_nodal",
        "verification_pressure_l2_relative_error = $(l2(pressure.uh, p_interp))",
        "verification_nutrient_nodal_max_error = $c_nodal",
        "verification_nutrient_l2_relative_error = $(l2(nutrient.uh, c_interp))",
        "verification_pressure_max = $(maximum(pressure.dofs))",
        "verification_pressure_max_exact = $(R0^2 / 4)"), "\n")
end

# =============================================================================
# 10. Dynamics and orchestration
# =============================================================================

"""
    dynamics(levelset, backdrop, cfg) -> NamedTuple

One evaluation of the dynamics: cut the current geometry, solve the nutrient problem, form
`g(c)`, solve the pressure problem with that source, and recover the nodal tissue velocity.
"""
function dynamics(levelset, backdrop, cfg)
    grid = levelset.discretization.levels[:init]
    t_start = time()
    snapshot = ED.cut(levelset, :init)
    cutinfo = cutGeometryInfo(snapshot, grid, backdrop.h)
    geo = buildGeometry(levelset.data[:init], backdrop, cfg)
    spaces = buildSpaces(geo)
    t_cut = time() - t_start

    nutrient = solveNutrient(geo, spaces, cfg, backdrop.h)
    growth = computeGrowthRate(nutrient, geo, spaces, cfg)
    pressure = solvePressure(geo, growth, spaces, cfg, backdrop.h)
    velocity = velocityFromPressure(pressure, geo, cutinfo, snapshot, backdrop)
    divergence = divergenceMismatch(velocity.cell_vx, velocity.cell_vy,
        backgroundNodal(growth.values, geo, backdrop), snapshot, backdrop)

    (; snapshot, cutinfo, geo, nutrient, growth, pressure, velocity, divergence, t_cut,
        c_nodes = backgroundNodal(nutrient.dofs, geo, backdrop),
        p_nodes = backgroundNodal(pressure.dofs, geo, backdrop),
        area = domainArea(geo), perimeter = interfaceLength(geo))
end

"""
    simulate(cfg) -> NamedTuple

The explicit orchestration loop. Everything physical arrives through `dynamics`; this
function integrates, book-keeps and measures the area-growth identity

    d|Omega|/dt = ∫_Ω g(c) dΩ

from the independently estimated reconstructed-area change and the source integral.
"""
function simulate(cfg)
    backdrop = backgroundModel(cfg)
    ordering_error = checkNodeOrdering(backdrop)

    grid = ED.Grid((backdrop.axes, backdrop.axes))
    discretization = ED.Discretization(Dict(:init => grid))
    levelset = ED.Field(discretization, Dict(:init => initialLevelSet(backdrop.axes, cfg)))
    marker_source = ED.Field(discretization, Dict(:init => initialMarker(backdrop.axes, cfg)))
    marker_target = ED.Field(discretization, Dict(:init => zeros(backdrop.nx, backdrop.nx)))
    cfg = merge(cfg, (; weno_cache = ED.WENO5Cache(),
                        cip_cache = ED.CIPCache(vec(marker_source.data[:init]), grid)))

    mkpath(cfg.outdir)
    metrics_path = joinpath(cfg.outdir, "metrics.csv")
    metrics = open(metrics_path, "w")
    println(metrics, join(("step", "t", "dt", "area_cut", "area_fem", "perimeter",
        "d_area_dt", "integral_growth", "balance_residual", "nutrient_residual",
        "pressure_residual", "flux_identity_residual", "interface_flux", "facet_speed_max",
        "divergence_l2", "divergence_linf",
        "divergence_nodes", "c_min", "c_max", "p_min", "p_max", "nitsche_boundary_defect",
        "components_nodes", "components_cells", "min_cut_fraction", "aggregate_roots",
        "vmax", "marker_min",
        "marker_max", "reinit_area_change", "t_cut", "t_nutrient_solve", "t_pressure_solve",
        "t_step"), ","))

    frames = Vector{Any}()
    t = 0.0
    step = 0
    merged_at = NaN
    min_fraction = Inf
    max_balance = 0.0
    sum_balance = 0.0
    max_divergence = 0.0
    max_nutrient = 0.0
    max_pressure = 0.0
    max_flux_identity = 0.0
    sum_flux_ratio = 0.0
    max_facet_speed = 0.0
    max_boundary_defect = 0.0
    sqrt_balance = 0.0
    min_marker = Inf
    max_marker = -Inf
    sum_reinit = 0.0
    max_reinit = 0.0
    t_cut_total = 0.0
    t_nutrient_total = 0.0
    t_pressure_total = 0.0
    dt_min = Inf
    dt_max = 0.0
    components_initial = 0
    area_initial = NaN
    run_start = time()

    while t < cfg.tfinal - 1.0e-12
        step += 1
        step_start = time()
        state = dynamics(levelset, backdrop, cfg)
        step == 1 && (components_initial = componentCounts(state.snapshot,
                                                           levelset.data[:init], backdrop);
                      area_initial = state.cutinfo.area)

        dt, vmax = timestep(state.velocity.velocity, backdrop, cfg, t)
        dt <= 0 && break
        dt_min = min(dt_min, dt)
        dt_max = max(dt_max, dt)

        advanceMarker!(marker_target, marker_source, state.velocity.velocity, dt, cfg)
        marker_source, marker_target = marker_target, marker_source

        new_snapshot, area_transported = advanceGeometry!(levelset, state.velocity.velocity,
            dt, cfg, grid, backdrop.h)

        # Area right after re-distancing, and the change that re-distancing alone produced.
        area_after = cutGeometryInfo(new_snapshot, grid, backdrop.h).area
        reinit_change = abs(area_after - area_transported)

        # Both ends of the area difference use the same (cut) reconstruction.
        observed_rate = (area_after - state.cutinfo.area) / dt
        balance = abs(observed_rate - state.growth.integral) /
                  max(abs(state.growth.integral), 1.0e-12)
        flux_residual = abs(state.pressure.flux_identity - state.growth.integral) /
                        max(abs(state.growth.integral), 1.0e-12)
        components = componentCounts(new_snapshot, levelset.data[:init], backdrop)
        isnan(merged_at) && components.nodes == 1 && (merged_at = t + dt)

        t += dt
        marker_values = marker_source.data[:init]
        step_time = time() - step_start

        min_fraction = min(min_fraction, state.cutinfo.min_fraction)
        max_balance = max(max_balance, balance)
        sum_balance += balance
        sqrt_balance += balance^2
        max_divergence = max(max_divergence, state.divergence.l2)
        max_nutrient = max(max_nutrient, state.nutrient.residual)
        max_pressure = max(max_pressure, state.pressure.residual)
        max_flux_identity = max(max_flux_identity, flux_residual)
        sum_flux_ratio += state.velocity.interface_flux / max(abs(state.growth.integral), eps(Float64))
        max_facet_speed = max(max_facet_speed, state.velocity.facet_speed_max)
        max_boundary_defect = max(max_boundary_defect, abs(state.nutrient.boundary_defect))
        min_marker = min(min_marker, minimum(marker_values))
        max_marker = max(max_marker, maximum(marker_values))
        sum_reinit += reinit_change
        max_reinit = max(max_reinit, reinit_change)
        t_cut_total += state.t_cut
        t_nutrient_total += state.nutrient.t_solve + state.nutrient.t_assembly
        t_pressure_total += state.pressure.t_solve + state.pressure.t_assembly

        @printf("%3d  t=%.4f dt=%.4f  A=%.6f P=%.4f  ∫g=%.4f dA/dt=%.4f bal=%.3f  comps=%d/%d  c∈[%.3f,%.3f] p∈[%.4f,%.4f] vmax=%.4f (%.2fs)\n",
            step, t, dt, area_after, state.perimeter, state.growth.integral, observed_rate,
            balance, components.nodes, components.cells,
            state.nutrient.cmin, state.nutrient.cmax,
            state.pressure.pmin, state.pressure.pmax, vmax, step_time)
        flush(stdout)

        println(metrics, join((step, t, dt, area_after, state.area, state.perimeter,
            observed_rate, state.growth.integral, balance, state.nutrient.residual,
            state.pressure.residual, flux_residual, state.velocity.interface_flux,
            state.velocity.facet_speed_max, state.divergence.l2,
            state.divergence.linf, state.divergence.nodes, state.nutrient.cmin,
            state.nutrient.cmax, state.pressure.pmin, state.pressure.pmax,
            state.nutrient.boundary_defect, components.nodes, components.cells,
            state.cutinfo.min_fraction,
            state.cutinfo.aggregate_roots, vmax, minimum(marker_values),
            maximum(marker_values), reinit_change, state.t_cut, state.nutrient.t_solve,
            state.pressure.t_solve, step_time), ","))
        flush(metrics)

        push!(frames, (; t, phi = copy(levelset.data[:init]), c = state.c_nodes,
            p = state.p_nodes, growth = backgroundNodal(state.growth.values, state.geo,
                backdrop), vx = copy(state.velocity.vx), vy = copy(state.velocity.vy),
            cell_vx = copy(state.velocity.cell_vx), cell_vy = copy(state.velocity.cell_vy),
            q = copy(marker_values), area = area_after, perimeter = state.perimeter,
            components_nodes = components.nodes, components_cells = components.cells,
            integral_growth = state.growth.integral, observed_rate,
            balance))
    end
    close(metrics)
    total_runtime = time() - run_start
    serialize(joinpath(cfg.outdir, "fields.jls"),
        (; axes = backdrop.axes, ncells = cfg.ncells, frames))
    (; metrics_path, total_runtime, steps = step, t_final = t,
        tfinal_reached = t >= cfg.tfinal - 1.0e-9, merged_at, ordering_error,
        components_initial, area_initial, final_area = frames[end].area,
        min_fraction, balance_mean = sum_balance / max(step, 1),
        balance_rms = sqrt(sqrt_balance / max(step, 1)), balance_max = max_balance,
        divergence_max = max_divergence, nutrient_residual_max = max_nutrient,
        pressure_residual_max = max_pressure, flux_identity_max = max_flux_identity,
        boundary_defect_max = max_boundary_defect,
        interface_flux_ratio_mean = sum_flux_ratio / max(step, 1),
        facet_speed_max = max_facet_speed, marker_min = min_marker,
        marker_max = max_marker, reinit_area_sum = sum_reinit, reinit_area_max = max_reinit,
        t_cut_total, t_nutrient_total, t_pressure_total, dt_min, dt_max,
        final_components = frames[end].components_nodes, frames = length(frames))
end

# =============================================================================
# 11. Entry point
# =============================================================================

function main(args = ARGS)
    ncells = length(args) >= 1 ? parse(Int, args[1]) : 64
    tfinal = length(args) >= 2 ? parse(Float64, args[2]) : 0.5
    name = length(args) >= 3 ? args[3] : @sprintf("n%d", ncells)
    dtmax = length(args) >= 4 ? parse(Float64, args[4]) : 5.0e-3
    cfg = defaultConfig(; ncells, tfinal, dtmax, outdir = joinpath(@__DIR__, "output", name))
    @printf("NutrientDrivenGrowth: %d^2 cells (h=%.5f), tfinal=%.3f, output=%s\n",
        ncells, 1.0 / ncells, tfinal, cfg.outdir)

    provenance = checkEvolvingDomains()
    for line in provenance
        println(line)
    end

    verification = analyticVerification(backgroundModel(cfg), cfg)
    result = simulate(cfg)

    report = open(joinpath(cfg.outdir, "run_metadata.txt"), "w")
    for line in provenance
        println(report, line)
    end
    println(report, verification)
    for (key, value) in pairs(result)
        key === :frames && continue
        println(report, key, " = ", value)
    end
    println(report, "ncells = ", cfg.ncells)
    println(report, "h = ", 1.0 / cfg.ncells)
    println(report, "tfinal_requested = ", cfg.tfinal)
    println(report, "dtmax = ", cfg.dtmax)
    println(report, "cfl = ", cfg.cfl)
    println(report, "Dc = ", cfg.Dc, " lambda = ", cfg.lambda)
    println(report, "g0 = ", cfg.g0, " cstar = ", cfg.cstar)
    println(report, "nitsche = ", cfg.nitsche, " ghost = ", cfg.ghost,
        " quadrature_degree = ", cfg.qdegree)
    close(report)
    @printf("done in %.1f s (%d steps, t=%.4f, merged_at=%s)\n", result.total_runtime,
        result.steps, result.t_final, result.merged_at)
    @printf("report: %s\n", joinpath(cfg.outdir, "run_metadata.txt"))
    return result
end

if abspath(PROGRAM_FILE) == @__FILE__
    main()
end
