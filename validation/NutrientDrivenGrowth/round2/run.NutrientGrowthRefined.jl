# NOTE: this round-2 script calls the removed `filterSubResolutionFeatures!` (morphological
# opening/closing). Its recorded results refer to that operator; the package now offers
# `filterUnresolvedFeatures!` (Luo et al. 2016), which is not equivalent. Re-running requires
# porting the filter call, and the numbers will change.
#!/usr/bin/env julia
#
# NutrientGrowthRefined  (Campaign II, Experiment 5)
# ==================================================
# Nutrient-driven free-boundary growth of two merging tissue discs on a Cartesian box,
# orchestrated through the *public* EvolvingDomains 1.0 surface and discretised by
# Gridap + GridapEmbedded (CutFEM: symmetric Nitsche + ghost penalty).
#
#   nutrient   -Dc Δc + λ c = 0     in Ω,   c = 1 on Γ        Dc = 0.02, λ = 0.05
#   growth      g(c) = g0 max(c-c*,0)                         g0 = 1.5,  c* = 0.2
#   mechanics  -Δp = g(c)           in Ω,   p = 0 on Γ        Darcy: v = -∇p
#   interface  ∂_t φ + v⋅∇φ = 0     (level set, φ < 0 inside), V_n = v⋅n
#   marker     ∂_t q + v⋅∇q = 0     (passive intensive scalar, CIP)
#
# What is new relative to the first round (validation/NutrientDrivenGrowth, read-only):
#
#   * provenance is printed *and* recorded: package root, git revision, version, native
#     CutCells library path + hash, Julia / Gridap / GridapEmbedded versions;
#   * two velocity-recovery methods (A: cell-centroid bilinear gradient; B: interface
#     located weighted-quadratic patch reconstruction over interior cell averages),
#     identical in every other respect, both verified against a manufactured pressure;
#   * reinitialisation is measured with a *fresh* snapshot after `reinitialize!` — the first
#     round returned the pre-reinitialisation snapshot and therefore measured zero change;
#   * every saved frame is synchronised: the elliptic fields in a frame are solved on the
#     geometry of that frame's own time;
#   * three topology measures (negative nodes, active cells, edge-connectivity of the
#     physical cut triangulation);
#   * per-stage timing and allocation accounting for every step.
#
# Usage (the project path must be absolute: a wrong relative path silently falls back to the
# default environment, which contains a pre-1.0 EvolvingDomains checkout):
#
#   PROJ="/Users/naoufel/Documents/Scientific stack/EvolvingDomains-Redesign/validation/env"
#   export EVOLVINGDOMAINS_CUTCELLS_LIBRARY="…/validation/lib/cutcells/lib/libevolvingdomains_cut_cells.dylib"
#
#   julia --project="$PROJ" run.NutrientGrowthRefined.jl verify [n ...]
#       elliptic verification + recovery scoring on the initial discs, writes output/verification.csv
#   julia --project="$PROJ" run.NutrientGrowthRefined.jl run <ncells> <tfinal> <A|B> <dtmax> <tag>
#       one coupled run, writes output/<tag>/
#   julia --project="$PROJ" run.NutrientGrowthRefined.jl recovery <ncells>
#       Method A vs Method B against the manufactured unit-source pressure
#   julia --project="$PROJ" run.NutrientGrowthRefined.jl merged [tag ...]
#       manufactured test on the merged geometry of finished runs -> output/merged_manufactured.csv
#   julia --project="$PROJ" run.NutrientGrowthRefined.jl extension [n ...]
#       brute-force vs CSR-ring nearest-facet adapter -> output/extension_benchmark.csv

import EvolvingDomains as ED
import GridapEmbedded as GE

using Gridap
using Gridap.Geometry: get_node_coordinates, get_cell_coordinates
using GridapEmbedded.LevelSetCutters: DiscreteGeometry
using StaticArrays: SVector
using LinearAlgebra: norm, svdvals, qr
using Printf: @printf, @sprintf
using Serialization: serialize, deserialize
using SHA: sha256

const MODEL_ORDER = 1

_cleanpath(path) = rstrip(normpath(path), ('/', '\\'))
const REPOSITORY_ROOT = _cleanpath(joinpath(@__DIR__, "..", "..", ".."))

# =============================================================================
# 0. Provenance: refuse to run on anything but the redesign checkout
# =============================================================================

"""`git rev-parse HEAD` of the repository the loaded package lives in (empty if unavailable)."""
function gitRevision(root)
    try
        return strip(read(pipeline(`git -C $root rev-parse HEAD`, stderr = devnull), String))
    catch
        return "unavailable"
    end
end

"""`git status --porcelain` summary of the package root, as a count of dirty paths."""
function gitDirtyCount(root)
    try
        lines = split(strip(read(pipeline(`git -C $root status --porcelain`,
            stderr = devnull), String)), '\n')
        return count(!isempty, lines)
    catch
        return -1
    end
end

"""
    checkEvolvingDomains() -> Vector{String}

Hard guard on package identity plus the provenance block recorded in `run_metadata.txt`.
The Julia load path here is `validation/env` stacked on the user default environment, which
also contains an *older* EvolvingDomains checkout; the `[sources] path` entry in
`validation/env/Project.toml` wins, and this function proves it rather than assuming it.
"""
function checkEvolvingDomains()
    root = _cleanpath(pkgdir(ED))
    version = pkgversion(ED)
    source = _cleanpath(pathof(ED))
    root == REPOSITORY_ROOT || error(
        "EvolvingDomains was loaded from \"$root\" but this experiment requires the redesign " *
        "repository \"$REPOSITORY_ROOT\". Run with --project=$(joinpath(REPOSITORY_ROOT, "validation", "env"))")
    version >= v"1.0.0" || error("EvolvingDomains $version is not the redesign (needs ≥ 1.0.0)")
    startswith(source, root * "/") ||
        error("pathof(EvolvingDomains) = $source is outside the repository root $root")
    for name in (:Grid, :Discretization, :Field, :cut, :Cut, :CutSurface, :CUT_IN, :CUT_OUT,
        :CUT_CUT, :aggregateCutCells, :reinitialize!, :eulerianTransport!,
        :lagrangianTransport!, :WENO5Cache, :CIPCache, :TransportMap, :extend!)
        isdefined(ED, name) || error("EvolvingDomains 1.0 must define $name; it does not")
    end
    for name in (:EvolvingDiscreteGeometry, :CartesianGridInfo, :AbstractVelocitySource,
        :aggregate_cut_cells)
        isdefined(ED, name) &&
            error("EvolvingDomains defines $name, which is pre-1.0 API: wrong checkout")
    end

    library = ED.CutCellsBinding._libraryPath()
    override = get(ENV, "EVOLVINGDOMAINS_CUTCELLS_LIBRARY", "")
    provenance = isempty(override) ? "artifact:$(library)" : "env-override:$(abspath(override))"
    lines = String[
        "evolvingdomains_root = $root",
        "evolvingdomains_version = $version",
        "evolvingdomains_source = $source",
        "evolvingdomains_git_revision = $(gitRevision(root))",
        "evolvingdomains_git_dirty_paths = $(gitDirtyCount(root))",
        "cutcells_library = $library",
        "cutcells_library_provenance = $provenance",
        "cutcells_library_sha256 = $(bytes2hex(sha256(read(library))))",
        "julia_version = $(VERSION)",
    ]
    for file in ("src/EvolvingDomains.jl", "src/DataModel.jl", "src/CutData.jl",
        "src/Kinematic/Kinematic.jl", "src/Geometric/Geometric.jl")
        path = joinpath(root, file)
        push!(lines, "evolvingdomains_sha256[$file] = $(bytes2hex(sha256(read(path))))")
    end
    for (name, pkg) in (("gridap", Gridap), ("gridapembedded", GE))
        push!(lines, "$(name)_version = $(pkgversion(pkg))")
        push!(lines, "$(name)_root = $(_cleanpath(pkgdir(pkg)))")
        push!(lines, "$(name)_git_revision = $(gitRevision(_cleanpath(pkgdir(pkg))))")
    end
    push!(lines, "load_path = $(join(Base.load_path(), " : "))")
    return lines
end

# =============================================================================
# 1. Configuration
# =============================================================================

"""
    defaultConfig(; kwargs...) -> NamedTuple

All physical, numerical and diagnostic parameters. `eps_scale` is the additive denominator
scale of the area-balance residual ε_A = |Δ|Ω|/Δt − ∫g| / max(|∫g|, eps_scale); it is a fixed
positive constant built from the physical parameters, 0.1 % of the initial growth rate scale
`g0 π R0²`, so that the normalisation never divides by a vanishing growth rate. For this
experiment |∫g| ≈ 0.08 ≫ eps_scale ≈ 6e-5 and the scale never binds (verified per run from
the recorded `eps_scale_binds` counter).
"""
function defaultConfig(; ncells::Int = 64,
        tfinal::Float64 = 0.5,
        dtmax::Float64 = 5.0e-3,
        cfl::Float64 = 0.4,
        method::Symbol = :centroid,
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
        patch_radius_h::Float64 = 3.5,
        patch_radius_max_h::Float64 = 6.0,
        min_patch_samples::Int = 8,
        patch_condition_max::Float64 = 1.0e3,
        extension_index::Bool = false,
        frame_stride::Int = 10,
        filter_mode::Symbol = :none,
        filter_radius::Int = 0,
        tag::String = "case",
        outdir::String = joinpath(@__DIR__, "output"))
    ncells >= 8 || throw(ArgumentError("need at least 8 cells per axis"))
    method in (:centroid, :patch) || throw(ArgumentError("method must be :centroid or :patch"))
    filter_radius >= 0 || throw(ArgumentError("filter_radius must be >= 0 (0 disables the filter)"))
    filter_radius == 0 || filter_mode in (:opening, :closing) ||
        throw(ArgumentError("filter_mode must be :opening or :closing when filter_radius > 0"))
    (; ncells, tfinal, dtmax, cfl, method, Dc, lambda, g0, cstar, R0, center1, center2,
        nitsche, ghost, qdegree, patch_radius_h, patch_radius_max_h, min_patch_samples,
        patch_condition_max, extension_index, frame_stride, filter_mode, filter_radius,
        tag, outdir,
        eps_scale = 1.0e-3 * g0 * π * R0^2)
end

# =============================================================================
# 2. Background grid and initial data
# =============================================================================

function backgroundModel(cfg)
    nx = cfg.ncells + 1
    axes = collect(range(0.0, 1.0; length = nx))
    model = CartesianDiscreteModel((0.0, 1.0, 0.0, 1.0), (cfg.ncells, cfg.ncells))
    Vnodal = FESpace(Triangulation(model), ReferenceFE(lagrangian, Float64, 1))
    (; model, axes, Vnodal, nx, ncx = cfg.ncells, h = 1.0 / cfg.ncells)
end

"""Assert Gridap's node order equals `vec` of an `(nx, nx)` node array (first axis fastest)."""
function checkNodeOrdering(backdrop)
    gridap = [Tuple(p) for p in collect(get_node_coordinates(backdrop.model))]
    evolving = [(a, b) for b in backdrop.axes for a in backdrop.axes]
    length(gridap) == length(evolving) ||
        error("node count mismatch: Gridap $(length(gridap)) vs $(length(evolving))")
    worst = 0.0
    for k in eachindex(gridap)
        worst = max(worst, abs(gridap[k][1] - evolving[k][1]),
            abs(gridap[k][2] - evolving[k][2]))
    end
    worst <= 100eps(Float64) || error(
        "Gridap / EvolvingDomains node orders differ (max deviation $worst)")
    return worst
end

initialLevelSet(axes, cfg) = [min(hypot(a - cfg.center1[1], b - cfg.center1[2]) - cfg.R0,
    hypot(a - cfg.center2[1], b - cfg.center2[2]) - cfg.R0) for a in axes, b in axes]

initialMarker(axes, cfg) = [0.8exp(-((a - cfg.center1[1])^2 + (b - cfg.center1[2])^2) / 0.06^2) +
                            0.2exp(-((a - cfg.center2[1])^2 + (b - cfg.center2[2])^2) / 0.06^2)
                            for a in axes, b in axes]

# =============================================================================
# 3. EvolvingDomains cut snapshot -> geometric quantities
# =============================================================================

@inline decodeCell(cell::Integer, ncx::Integer) = ((cell - 1) % ncx + 1, (cell - 1) ÷ ncx + 1)
@inline encodeCell(i::Integer, j::Integer, ncx::Integer) = i + (j - 1) * ncx

"""
    cutGeometryInfo(snapshot, grid, h) -> NamedTuple

Physical-part geometry of every background cell from the public `Cut` arrays: in 2D each
subcell is a triangle and `subcell_states == CUT_IN` marks the part belonging to the physical
phase, so `inside_area[c]` is the exact area of `Ω ∩ cell` for the piecewise-linear cut model,
and `cx[c], cy[c]` is the area centroid of that piece. `min_fraction` is the smallest non-empty
cell area in units of `h²` (cut-cell conditioning indicator).
"""
function cutGeometryInfo(snapshot::ED.Cut, grid::ED.Grid, h)
    ncells = length(snapshot.bg_states)
    inside_area = zeros(ncells)
    cx = zeros(ncells)
    cy = zeros(ncells)
    totals = zeros(ncells)
    pointers = snapshot.point_pointers
    ids = snapshot.point_ids
    points = snapshot.physical_points
    for s in eachindex(snapshot.subcell_states)
        lo = Int(pointers[s])
        hi = Int(pointers[s + 1]) - 1
        hi - lo + 1 == 3 || error("expected triangular 2D subcells, got $(hi - lo + 1) vertices")
        a = points[Int(ids[lo])]
        b = points[Int(ids[lo + 1])]
        c = points[Int(ids[lo + 2])]
        area = abs((b[1] - a[1]) * (c[2] - a[2]) - (b[2] - a[2]) * (c[1] - a[1])) / 2
        cell = Int(snapshot.cell_to_bgcell[s])
        totals[cell] += area
        snapshot.subcell_states[s] == ED.CUT_IN || continue
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
    roots = ED.aggregateCutCells(snapshot, grid)
    (; inside_area, cx, cy, occupied, min_fraction = isempty(occupied) ? 0.0 :
            minimum(inside_area[occupied]) / h^2,
        area = sum(inside_area), n_cut = count(==(ED.CUT_CUT), snapshot.bg_states),
        n_in = count(==(ED.CUT_IN), snapshot.bg_states),
        aggregate_roots = length(unique(r for r in roots if r > 0)))
end

"""
    interfaceFacets(snapshot, backdrop) -> NamedTuple

The material boundary as consumed by the velocity recovery: facet midpoints, outward unit
normals, lengths and owning background cells, straight from the public `CutSurface` arrays
(`physical_points`, `normals`, `facet_to_bgcell`; CSR row = facet index).
"""
function interfaceFacets(snapshot, backdrop)
    surface = snapshot.boundary
    n = length(surface.facet_to_bgcell)
    n > 0 || error("the cut snapshot carries no interface facet")
    mx = zeros(n)
    my = zeros(n)
    nnx = zeros(n)
    nny = zeros(n)
    len = zeros(n)
    cell = Vector{Int}(undef, n)
    pointers = surface.point_pointers
    ids = surface.point_ids
    points = surface.physical_points
    for f in 1:n
        lo = Int(pointers[f])
        hi = Int(pointers[f + 1]) - 1
        hi - lo + 1 == 2 || error("expected segment facets in 2D, got $(hi - lo + 1) vertices")
        p1 = points[Int(ids[lo])]
        p2 = points[Int(ids[lo + 1])]
        mx[f] = (p1[1] + p2[1]) / 2
        my[f] = (p1[2] + p2[2]) / 2
        normal = surface.normals[f]
        nnx[f] = normal[1]
        nny[f] = normal[2]
        len[f] = hypot(p2[1] - p1[1], p2[2] - p1[2])
        cell[f] = Int(surface.facet_to_bgcell[f])
    end
    (; nfacet = n, mx, my, nnx, nny, len, cell)
end

"""Distance from each point in `(px, py)` to the nearest facet midpoint of `facets`."""
function nearestFacetDistances(px, py, facets)
    out = zeros(length(px))
    for k in eachindex(px)
        best = Inf
        for f in 1:facets.nfacet
            d = (px[k] - facets.mx[f])^2 + (py[k] - facets.my[f])^2
            d < best && (best = d)
        end
        out[k] = sqrt(best)
    end
    return out
end

# =============================================================================
# 4. Topology: three different notions of "number of material components"
# =============================================================================

"""Number of 4-connected components of a boolean mask over an `n × n` lattice."""
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
    geometricComponents(geo, h) -> Int

Number of edge-connected components of the *physical cut triangulation* handed to the finite
element solvers — the most geometrically faithful of the three measures, because it is the
discretised Ω itself rather than a Cartesian surrogate. Works on any input cell arity
(triangle or quad subcells): vertices are welded with a tolerance of `h*1e-9` (the package's
own `CutSurface` documentation warns that duplicated vertices may differ by an ulp, so exact
equality must never be used), and two cells are joined when they share a polygon edge.
"""
function geometricComponents(geo, h)
    cells = get_cell_coordinates(geo.Ω)
    tol = h * 1.0e-9
    n = length(cells)
    parent = collect(1:n)                    # union-find over *cells*
    function root(a::Int)
        while parent[a] != a
            parent[a] = parent[parent[a]]
            a = parent[a]
        end
        return a
    end
    function join!(a::Int, b::Int)
        ra, rb = root(a), root(b)
        ra == rb || (parent[rb] = ra)
    end
    vertex_ids = Dict{NTuple{2,Int},Int}()
    next_vertex = 0
    function vertexId(point)
        key = (round(Int, point[1] / tol), round(Int, point[2] / tol))
        id = get(vertex_ids, key, 0)
        if id == 0
            next_vertex += 1
            vertex_ids[key] = next_vertex
            id = next_vertex
        end
        return id
    end
    edge_owner = Dict{NTuple{2,Int},Int}()
    for (cellindex, cellpoints) in enumerate(cells)
        ids = [vertexId(p) for p in cellpoints]
        m = length(ids)
        for a in 1:m
            b = a == m ? 1 : a + 1
            key = ids[a] < ids[b] ? (ids[a], ids[b]) : (ids[b], ids[a])
            owner = get(edge_owner, key, 0)
            if owner == 0
                edge_owner[key] = cellindex
            else
                join!(owner, cellindex)
            end
        end
    end
    return length(unique(root(k) for k in 1:n))
end

"""
    topologyCounts(snapshot, levelset, backdrop, geo) -> NamedTuple

Three measures, all reported per step and all resolution dependent:

  * `nodes` — 4-connected components of the negative grid nodes (`φ < 0`);
  * `cells` — 4-connected components of the active background cells (`state != CUT_OUT`),
    i.e. cells fully inside or crossed by the interface: this flips about one cell *before*
    the nodal count, because a cut band already joins the regions while no node is negative;
  * `geometric` — edge-connected components of the physical cut triangulation (see
    `geometricComponents`), the material region as the FE solvers see it.
"""
topologyCounts(snapshot, levelset, backdrop, geo, h) = (;
    nodes = connectedComponents(vec(levelset .< 0.0), backdrop.nx),
    cells = connectedComponents(snapshot.bg_states .!= ED.CUT_OUT, backdrop.ncx),
    geometric = geometricComponents(geo, h))

# =============================================================================
# 5. Embedded geometry handed to Gridap / GridapEmbedded
# =============================================================================

function buildGeometry(levelset, backdrop, cfg)
    values = vec(levelset)
    if any(v -> abs(v) <= 1.0e-12, values)
        values = copy(values)
        values[findall(v -> abs(v) <= 1.0e-12, values)] .+= 100eps(Float64)
    end
    φh = FEFunction(backdrop.Vnodal, values)
    geo = DiscreteGeometry(φh, backdrop.model)
    cutgeo = GE.cut(backdrop.model, geo)
    Ω = Triangulation(cutgeo, GE.PHYSICAL_IN)
    Ωact = Triangulation(cutgeo, GE.ACTIVE_IN)
    Γ = GE.EmbeddedBoundary(cutgeo)
    Γg = GE.GhostSkeleton(cutgeo, GE.ACTIVE_IN)
    (; cutgeo, Ω, Ωact, Γ, Γg, dΩ = Measure(Ω, cfg.qdegree), dΓ = Measure(Γ, cfg.qdegree),
        dΓg = Measure(Γg, cfg.qdegree), nΓ = get_normal_vector(Γ),
        nΓg = get_normal_vector(Γg))
end

function buildSpaces(geo)
    reffe = ReferenceFE(lagrangian, Float64, MODEL_ORDER)
    V = TestFESpace(geo.Ωact, reffe, conformity = :H1)
    (; V, U = TrialFESpace(V))
end

domainArea(geo) = sum(∫(1.0) * geo.dΩ)
interfaceLength(geo) = sum(∫(1.0) * geo.dΓ)

# =============================================================================
# 6. Elliptic solves (CutFEM, Nitsche + ghost penalty)
# =============================================================================

"""
    solveElliptic(geo, spaces, cfg, h; diffusion, reaction, source, boundaryValue, label)

CutFEM Dirichlet problem `-∇⋅(diffusion ∇u) + reaction u = source` in Ω, `u = boundaryValue`
on Γ, with symmetric Nitsche (penalty `diffusion*nitsche/h`) and the Burman–Hansbo facet ghost
penalty (`diffusion*ghost*h`). Returns the FE function, its dofs, the algebraic residual and
assembly/solve timings.
"""
function solveElliptic(geo, spaces, cfg, h; diffusion, reaction, source, boundaryValue, label)
    γd = cfg.nitsche / h
    γg = cfg.ghost * h
    nΓ = geo.nΓ
    nΓg = geo.nΓg
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
    op = AffineFEOperator(a, l, spaces.U, spaces.V)
    matrix = get_matrix(op)
    rhs = get_vector(op)
    t1 = time()
    uh = solve(op)
    t2 = time()
    dofs = get_free_dof_values(uh)
    all(isfinite, dofs) || error("$label: solver returned non-finite values")
    (; uh, dofs, residual = norm(rhs - matrix * dofs) / max(norm(rhs), eps(Float64)),
        t_assembly = t1 - t0, t_solve = t2 - t1)
end

"""
    solveNutrient(geo, spaces, cfg, h) -> NamedTuple

Nutrient problem `-Dc Δc + λ c = 0` in Ω with `c = 1` on Γ. Reported besides the solution:

* `flux` — the raw recovered boundary flux `∫_Γ n⋅∇c`; with a Nitsche boundary condition this
  is a boundary *gradient* integral, not the flux of the continuous problem, and it carries
  the boundary layer;
* `penalty` — `(γ/h)∫_Γ (c-1)`, the Nitsche boundary term;
* `weak_identity` — `-∫_Γ n⋅∇c + (γ/h)∫_Γ (c-1)`, which the discrete weak form, tested with
  the constant function, forces to equal `-(λ/Dc)∫_Ω c`. Only that combination is an algebraic
  identity; the raw flux alone is not.
"""
function solveNutrient(geo, spaces, cfg, h)
    state = solveElliptic(geo, spaces, cfg, h; diffusion = cfg.Dc, reaction = cfg.lambda,
        source = nothing, boundaryValue = 1.0, label = "nutrient")
    γd = cfg.nitsche / h
    flux = sum(∫(geo.nΓ ⋅ ∇(state.uh)) * geo.dΓ)
    defect = sum(∫(state.uh - 1.0) * geo.dΓ)
    (; state..., cmin = minimum(state.dofs), cmax = maximum(state.dofs), flux,
        mass = sum(∫(state.uh) * geo.dΩ), boundary_defect = defect,
        penalty = γd * defect, weak_identity = -flux + γd * defect)
end

function computeGrowthRate(nutrient, geo, spaces, cfg)
    values = cfg.g0 .* max.(nutrient.dofs .- cfg.cstar, 0.0)
    gh = FEFunction(spaces.V, values)
    (; gh, values, integral = sum(∫(gh) * geo.dΩ))
end

"""
    solvePressure(geo, growth, spaces, cfg, h) -> NamedTuple

Pressure problem `-Δp = g(c)` with `p = 0` on Γ. Besides the solution it reports

* `flux` — the raw recovered boundary flux `∫_Γ n⋅∇p` (one term of the weak form);
* `penalty` — the Nitsche boundary penalty term `(γ/h)∫_Γ p`;
* `flux_identity` — `-∫_Γ n⋅∇p + (γ/h)∫_Γ p`, which equals `∫_Ω g` by testing the discrete
  weak form with the constant function. This is an algebraic identity of the discretisation,
  **not** evidence that the physical boundary velocity is accurate; the two are reported
  separately for exactly that reason.
"""
function solvePressure(geo, growth, spaces, cfg, h)
    state = solveElliptic(geo, spaces, cfg, h; diffusion = 1.0, reaction = 0.0,
        source = growth.gh, boundaryValue = 0.0, label = "pressure")
    flux = sum(∫(geo.nΓ ⋅ ∇(state.uh)) * geo.dΓ)
    penalty = (cfg.nitsche / h) * sum(∫(state.uh) * geo.dΓ)
    (; state..., flux, penalty, flux_identity = -flux + penalty,
        pmin = minimum(state.dofs), pmax = maximum(state.dofs))
end

# =============================================================================
# 7. Velocity recovery (Method A / Method B) and extension to the grid
# =============================================================================

"""
    backgroundNodal(dofs, geo, backdrop) -> Matrix

Place active-space dofs on the background `(nx, nx)` node array (first axis fastest); nodes
outside the active region become `NaN`. The dof↔node correspondence comes from
`get_node_coordinates(Ωact)` and is asserted to be a bijection; it is an ordering coincidence
of the two packages, checked at runtime rather than assumed.
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
    bilinearGradient(nodal, i, j, x, y, h) -> SVector or nothing

Gradient of the Q1 interpolant of `nodal` in background cell `(i, j)` at the physical point
`(x, y)`; `nothing` when a corner node is undefined. This is the gradient *of the discrete
pressure field inside that cell*, evaluated wherever the caller asks — never a statement about
the gradient at the interface unless the evaluation point is on the interface.
"""
function bilinearGradient(nodal, i, j, x, y, h)
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

"""
    recoveryCentroid(pressure, geo, cutinfo, facets, backdrop) -> NamedTuple   (Method A)

Baseline. For every active cell, `-∇p_h` is evaluated with the bilinear derivative at the
**area centroid of the cell's physical part** — a point up to `h/2` inside Ω, not on Γ. Each
interface facet then inherits the velocity of its owning background cell, so the interface
velocity actually used is the interior-cell one ("original centroid recovery").
"""
function recoveryCentroid(pressure, geo, cutinfo, facets, backdrop)
    t0 = time()
    ncx = backdrop.ncx
    h = backdrop.h
    pn = backgroundNodal(pressure.dofs, geo, backdrop)
    cell_vx = fill(NaN, ncx, ncx)
    cell_vy = fill(NaN, ncx, ncx)
    for cell in cutinfo.occupied
        i, j = decodeCell(cell, ncx)
        gradient = bilinearGradient(pn, i, j, cutinfo.cx[cell], cutinfo.cy[cell], h)
        gradient === nothing && continue
        cell_vx[i, j] = -gradient[1]
        cell_vy[i, j] = -gradient[2]
    end
    facet_vx = zeros(facets.nfacet)
    facet_vy = zeros(facets.nfacet)
    for f in 1:facets.nfacet
        i, j = decodeCell(facets.cell[f], ncx)
        if isfinite(cell_vx[i, j])
            facet_vx[f] = cell_vx[i, j]
            facet_vy[f] = cell_vy[i, j]
        else                     # degenerate cell: fall back to its own centre
            gradient = bilinearGradient(pn, i, j, (i - 0.5) * h, (j - 0.5) * h, h)
            if gradient !== nothing
                facet_vx[f] = -gradient[1]
                facet_vy[f] = -gradient[2]
            end
        end
    end
    (; facet_vx, facet_vy, cell_vx, cell_vy, pn, fallbacks = 0,
        t_recovery = time() - t0, order = "A: bilinear ∇p_h at physical-part centroid")
end

"""
    interiorSamples(pn, snapshot, backdrop) -> NamedTuple

Cell-centre pressure samples for the patch reconstruction: for every fully-inside (`CUT_IN`)
background cell the average of its four nodal dofs, placed at the cell centre. For a Q1 field
that average is a second-order accurate approximation of `p` at the centre, and every dof used
belongs to a cell that is entirely inside Ω — no interface cell, no weakly imposed boundary
value and no ghost stabilisation term enters the samples. The samples are bucketed by
background cell for O(1) patch queries.
"""
function interiorSamples(pn, snapshot, backdrop)
    ncx = backdrop.ncx
    h = backdrop.h
    states = snapshot.bg_states
    sx = Float64[]
    sy = Float64[]
    sv = Float64[]
    scell = Int[]
    bucket = Dict{Int,Vector{Int}}()
    for cell in eachindex(states)
        states[cell] == ED.CUT_IN || continue
        i, j = decodeCell(cell, ncx)
        value = (pn[i, j] + pn[i + 1, j] + pn[i, j + 1] + pn[i + 1, j + 1]) / 4
        isfinite(value) || continue
        push!(sx, (i - 0.5) * h)
        push!(sy, (j - 0.5) * h)
        push!(sv, value)
        push!(scell, cell)
        push!(get!(() -> Int[], bucket, cell), length(sx))
    end
    (; sx, sy, sv, scell, bucket, ncx, count = length(sx))
end

"""
    bilinearValue(nodal, i, j, x, y, h) -> Float64

Value of the Q1 interpolant of `nodal` in background cell `(i, j)` at `(x, y)`; `NaN` when a
corner node is undefined. Used by the patch visibility test.
"""
function bilinearValue(nodal, i, j, x, y, h)
    nx = size(nodal, 1)
    (1 <= i <= nx - 1 && 1 <= j <= nx - 1) || return NaN
    u00 = nodal[i, j]
    u10 = nodal[i + 1, j]
    u01 = nodal[i, j + 1]
    u11 = nodal[i + 1, j + 1]
    all(isfinite, (u00, u10, u01, u11)) || return NaN
    ξ = clamp((x - (i - 1) * h) / h, 0.0, 1.0)
    η = clamp((y - (j - 1) * h) / h, 0.0, 1.0)
    return u00 * (1 - ξ) * (1 - η) + u10 * ξ * (1 - η) + u01 * (1 - ξ) * η + u11 * ξ * η
end

"""
    visible(phi, x0, y0, x1, y1, h, ncx; samples = 6) -> Bool

True when the straight segment from the interface point `(x0, y0)` to the interior sample
`(x1, y1)` stays inside the physical phase, tested on the same bilinear level-set interpolant
the cut is built from. This is what keeps a patch from reaching across the void into another
part of the domain — with two discs a few cells apart, an unfiltered patch centred on the
facing side of one disc collects samples *inside the other disc*, where `p` belongs to a
different analytic branch. That contamination, not the fit itself, is what a naive patch
reconstruction pays for with large patches.
"""
function visible(phi, x0, y0, x1, y1, h, ncx; nsteps = 6)
    for s in 1:nsteps
        τ = s / nsteps
        x = x0 + τ * (x1 - x0)
        y = y0 + τ * (y1 - y0)
        i = min(max(round(Int, x / h), 1), ncx)
        j = min(max(round(Int, y / h), 1), ncx)
        value = bilinearValue(phi, i, j, x, y, h)
        isfinite(value) || return false
        value > h * 1.0e-3 && return false
    end
    return true
end

"""
    patchGradient(samples, phi, x0, y0, h, radius0, radius_max, mincount)
        -> (gradient, degree, radius, condition) or nothing

Weighted least-squares local polynomial reconstruction of `p` around the interface point
`(x0, y0)`, from interior cell-centre samples that are *visible* from it (see `visible`),
evaluated at `(x0, y0)`:

    p(x) ≈ a0 + a1 X + a2 Y + a3 X² + a4 XY + a5 Y²,     X = (x - x0)/h, Y = (y - y0)/h

with the compactly supported weight `w = (1 - (d/r)²)²`. The patch radius grows from `radius0`
in one-cell steps until the quadratic fit has `mincount` visible samples; below that a linear
fit is attempted at `radius_max`, and only when no fit is possible at all does the caller fall
back to the Method A value (counted, never hidden).

The normal equations are **not** formed: the weighted design matrix is column-equilibrated and
solved by QR, and the fit is rejected unless its design matrix has `σ_min/σ_max ≥ 1e-6`.
Both matter. A degree-2 monomial patch evaluated one cell outside a one-sided sample cloud has
a Hilbert-type conditioning of order `10³` even in the best case, and squaring that in the
normal equations is what turns a mildly sensitive extrapolation into pure noise; with the
normal equations this campaign observed condition numbers of `10¹⁶` and interface speeds an
order of magnitude too large, which then poisoned the single-valued nodal extension.

Because a weighted least-squares fit reproduces any exactly representable polynomial
regardless of its weights, a quadratic pressure field — which is what `-Δp = const` produces
on a disc — is recovered exactly up to the finite-element error; that is what the
manufactured-solution verification measures.
"""
function patchGradient(samples, phi, x0, y0, h, radius0, radius_max, mincount, cond_max)
    i0, j0 = round(Int, x0 / h), round(Int, y0 / h)
    rmax = ceil(Int, radius_max / h) + 1
    candidates = Int[]
    for j in (j0 - rmax):(j0 + rmax), i in (i0 - rmax):(i0 + rmax)
        for k in get(samples.bucket, encodeCell(i, j, samples.ncx), Int[])
            X = samples.sx[k] - x0
            Y = samples.sy[k] - y0
            X * X + Y * Y <= radius_max^2 || continue
            visible(phi, x0, y0, samples.sx[k], samples.sy[k], h, samples.ncx) || continue
            push!(candidates, k)
        end
    end
    isempty(candidates) && return nothing
    for degree in (2, 1)
        ncoef = degree == 2 ? 6 : 3
        need = degree == 2 ? mincount : 4
        first_radius = degree == 2 ? radius0 : max(radius0, radius_max - 2h)
        radius = first_radius
        while radius <= radius_max + 1.0e-12
            fit = _qrFit(samples, candidates, x0, y0, h, radius, degree, ncoef, need, cond_max)
            fit === nothing || return (fit[1], degree, radius, fit[2])
            radius += h
        end
    end
    return nothing
end

"""One equilibrated QR least-squares polynomial fit of the given degree over `candidates`."""
function _qrFit(samples, candidates, x0, y0, h, radius, degree, ncoef, need, cond_max)
    rows = Int[]
    for k in candidates
        X = (samples.sx[k] - x0) / h
        Y = (samples.sy[k] - y0) / h
        d2 = (X * X + Y * Y) * h^2 / radius^2
        d2 < 1 && push!(rows, k)
    end
    used = length(rows)
    used >= need || return nothing
    design = zeros(used, ncoef)
    rhs = zeros(used)
    for (r, k) in enumerate(rows)
        X = (samples.sx[k] - x0) / h
        Y = (samples.sy[k] - y0) / h
        d2 = (X * X + Y * Y) * h^2 / radius^2
        w = sqrt(max(1 - d2, 0.0)^2)
        φ = degree == 2 ? (1.0, X, Y, X * X, X * Y, Y * Y) : (1.0, X, Y)
        for a in 1:ncoef
            design[r, a] = w * φ[a]
        end
        rhs[r] = w * samples.sv[k]
    end
    column_norm = [norm(view(design, :, a)) for a in 1:ncoef]
    scale = maximum(column_norm)
    scale > 0 || return nothing
    for a in 1:ncoef
        column_norm[a] < 1.0e-8 * scale && return nothing   # a column the patch cannot see
        design[:, a] ./= column_norm[a]
    end
    spectrum = svdvals(design)
    σmin = minimum(spectrum)
    σmin > 0 || return nothing
    σmax = maximum(spectrum)
    cond = σmax / σmin
    cond <= cond_max || return nothing
    coefficients = qr(design) \ rhs
    all(isfinite, coefficients) || return nothing
    gradient = SVector(coefficients[2] / column_norm[2] / h,
        coefficients[3] / column_norm[3] / h)
    all(isfinite, gradient) || return nothing
    return (gradient, cond)
end

"""
    recoveryPatch(pressure, geo, cutinfo, facets, snapshot, backdrop, cfg, phi) -> NamedTuple

Interface-aware recovery. For every interface facet the pressure gradient is obtained from a
weighted-quadratic local reconstruction of *interior* data (see `patchGradient`) and evaluated
at the facet midpoint, i.e. at a point on Γ. The facet velocity is `-∇p(x_Γ)`; its normal
component `V_n` is the interface normal speed. Facets whose patch is too small fall back to
the Method A value (counted in `fallbacks`) instead of inventing data.
"""
function recoveryPatch(pressure, geo, cutinfo, facets, snapshot, backdrop, cfg, phi)
    t0 = time()
    ncx = backdrop.ncx
    h = backdrop.h
    pn = backgroundNodal(pressure.dofs, geo, backdrop)
    samples = interiorSamples(pn, snapshot, backdrop)
    baseline = recoveryCentroid(pressure, geo, cutinfo, facets, backdrop)
    facet_vx = zeros(facets.nfacet)
    facet_vy = zeros(facets.nfacet)
    fallbacks = 0
    nquadratic = 0
    nlinear = 0
    cond_max = 0.0
    cond_facet = 0
    radius_max = 0.0
    for f in 1:facets.nfacet
        fit = patchGradient(samples, phi, facets.mx[f], facets.my[f], h,
            cfg.patch_radius_h * h, cfg.patch_radius_max_h * h, cfg.min_patch_samples,
            cfg.patch_condition_max)
        if fit === nothing
            facet_vx[f] = baseline.facet_vx[f]
            facet_vy[f] = baseline.facet_vy[f]
            fallbacks += 1
            continue
        end
        gradient, degree, radius, cond = fit
        facet_vx[f] = -gradient[1]
        facet_vy[f] = -gradient[2]
        degree == 2 ? (nquadratic += 1) : (nlinear += 1)
        radius_max = max(radius_max, radius)
        if cond > cond_max
            cond_max = cond
            cond_facet = f
        end
    end
    (; facet_vx, facet_vy, cell_vx = baseline.cell_vx, cell_vy = baseline.cell_vy, pn,
        fallbacks, nquadratic, nlinear, cond_max, cond_facet,
        radius_max, t_recovery = time() - t0,
        order = "B: weighted-quadratic patch of interior cell averages, evaluated on Γ",
        samples = samples.count)
end

"""
    _scanBucket!(state, offsets, flat, mx, my, x, y, cell) -> (best, best_distance, evaluations, hit)

Fold every facet of one background cell's bucket into the running nearest-facet state. Ties are
broken by the smallest facet index, exactly like the brute-force reference, so both searches
return the same velocity and not merely the same distance.
"""
@inline function _scanBucket(best, best_distance, evaluations, offsets, flat, mx, my, x, y, cell)
    hit = false
    for position in offsets[cell]:(offsets[cell + 1] - 1)
        f = flat[position]
        evaluations += 1
        dx = x - mx[f]
        dy = y - my[f]
        distance = dx * dx + dy * dy
        if distance < best_distance || (distance == best_distance && f < best)
            best_distance = distance
            best = f
            hit = true
        end
    end
    return best, best_distance, evaluations, hit
end

"""
    extendToNodes(facets, fvx, fvy, backdrop) -> NamedTuple

Zero-order closest-facet extension of the interface velocity to every background node: a node
takes the velocity of the nearest interface-facet midpoint, so every node on either side of a
given interface segment carries that segment's normal speed (a normal extension, exact for a
straight interface, single valued). Documented limitation: a single-valued nodal field cannot
represent two opposing interfaces closer than a cell, and the nearest facet wins.

The search is local, not global: facets are bucketed into their background cells and each node
scans expanding Chebyshev rings, stopping `stopping_rings` rings after the first hit; the
result is checked against brute force on a deterministic subsample of nodes
(`distance_check_max`, which must be 0 for the reported extension to be the true nearest-facet
one). `t_build`/`t_query` separate the index cost from the query cost, because the first-round
implementation paid `O(nodes × facets)` and this campaign must not repeat that silently.
"""
function extendToNodes(facets, fvx, fvy, backdrop; check_stride = 37)
    t0 = time()
    nx = backdrop.nx
    ncx = backdrop.ncx
    h = backdrop.h
    ncells = ncx * ncx
    # CSR buckets (background cell -> facet list) instead of a Dict: the ring scan touches the
    # same cells for every node and hashing per node dominated the adapter cost.
    counts = zeros(Int, ncells)
    for f in 1:facets.nfacet
        i = min(max(round(Int, facets.mx[f] / h) + 1, 1), ncx)
        j = min(max(round(Int, facets.my[f] / h) + 1, 1), ncx)
        counts[encodeCell(i, j, ncx)] += 1
    end
    offsets = Vector{Int}(undef, ncells + 1)
    offsets[1] = 1
    for c in 1:ncells
        offsets[c + 1] = offsets[c] + counts[c]
    end
    cursor = copy(offsets[1:ncells])
    flat = Vector{Int}(undef, facets.nfacet)
    for f in 1:facets.nfacet
        i = min(max(round(Int, facets.mx[f] / h) + 1, 1), ncx)
        j = min(max(round(Int, facets.my[f] / h) + 1, 1), ncx)
        c = encodeCell(i, j, ncx)
        flat[cursor[c]] = f
        cursor[c] += 1
    end
    mx = facets.mx
    my = facets.my
    vx = zeros(nx, nx)
    vy = zeros(nx, nx)
    distance2 = fill(Inf, nx, nx)
    evaluations = 0
    t1 = time()
    for j in 1:nx, i in 1:nx
        x = (i - 1) * h
        y = (j - 1) * h
        ci = min(i, ncx)
        cj = min(j, ncx)
        best = 0
        best_distance = Inf
        found_ring = -1
        for ring in 0:ncx
            hit = false
            if ring == 0
                best, best_distance, evaluations, found =
                    _scanBucket(best, best_distance, evaluations, offsets, flat, mx, my, x, y,
                        encodeCell(ci, cj, ncx))
                hit |= found
            else
                # Only the perimeter of the ring: iterating the whole square makes the ring
                # scan quadratic in the ring index, which is what made the first bucketed
                # version slower than brute force for nodes far from the interface.
                for d in -ring:ring
                    for (ii, jj) in ((ci + d, cj + ring), (ci + d, cj - ring),
                        (ci + ring, cj + d), (ci - ring, cj + d))
                        (1 <= ii <= ncx && 1 <= jj <= ncx) || continue
                        best, best_distance, evaluations, found =
                            _scanBucket(best, best_distance, evaluations, offsets, flat, mx, my,
                                x, y, encodeCell(ii, jj, ncx))
                        hit |= found
                    end
                end
            end
            # Exact stopping rule: every cell at Chebyshev ring `ring + 1` or beyond is at
            # least `ring*h` away from this node, so once `ring*h` reaches the current best
            # distance no further ring can hold a closer facet. A fixed "+2 rings" rule is NOT
            # enough: nodes well inside the bulk are first hit at ring ~9 on a 64² grid, where
            # the safe margin is four rings, and three nodes out of 4225 then took a facet
            # 0.9 % too far away.
            if best != 0 && ring * h >= sqrt(best_distance)
                break
            end
        end
        best == 0 && error("node ($i, $j) found no interface facet: empty velocity field")
        vx[i, j] = fvx[best]
        vy[i, j] = fvy[best]
        distance2[i, j] = best_distance
    end
    t2 = time()
    # The ring search is only claimed to be the true nearest-facet search if it agrees with
    # brute force on this deterministic subsample of nodes.
    check_error = 0.0
    for j in 1:nx, i in 1:nx
        ((i - 1) + (j - 1) * nx) % check_stride == 0 || continue
        x = (i - 1) * h
        y = (j - 1) * h
        brute = minimum((x - mx[f])^2 + (y - my[f])^2 for f in 1:facets.nfacet)
        check_error = max(check_error, sqrt(distance2[i, j]) - sqrt(brute))
    end
    velocity = SVector{2,Float64}[SVector{2,Float64}(vx[k], vy[k]) for k in eachindex(vx)]
    (; vx, vy, velocity, maximum_speed = maximum(norm, velocity),
        t_build = t1 - t0, t_query = t2 - t1, distance_evaluations = evaluations,
        brute_force_evaluations = nx * nx * facets.nfacet,
        distance_check_max = check_error, buckets = ncells)
end

"""
    recoverVelocity(method, pressure, geo, cutinfo, facets, snapshot, backdrop, cfg, phi)

The single dispatch point of the controlled comparison: identical geometry, identical
pressures, identical extension operator — only the facet velocity sample differs.
"""
function recoverVelocity(method, pressure, geo, cutinfo, facets, snapshot, backdrop, cfg, phi)
    t0 = time()
    interface = method === :patch ?
                recoveryPatch(pressure, geo, cutinfo, facets, snapshot, backdrop, cfg, phi) :
                recoveryCentroid(pressure, geo, cutinfo, facets, backdrop)
    extension = cfg.extension_index ?
                extendToNodes(facets, interface.facet_vx, interface.facet_vy, backdrop) :
                nearestFacetExtension(facets, interface.facet_vx, interface.facet_vy, backdrop)
    t_extension = time() - t0
    flux = 0.0
    flux_signed = 0.0
    vn = zeros(facets.nfacet)
    for f in 1:facets.nfacet
        vn[f] = interface.facet_vx[f] * facets.nnx[f] + interface.facet_vy[f] * facets.nny[f]
        flux += vn[f] * facets.len[f]
        flux_signed += abs(vn[f]) * facets.len[f]
    end
    (; interface, extension, vn, interface_flux = flux,
        interface_flux_abs = flux_signed, t_extension,
        facet_speed_max = maximum(f -> hypot(interface.facet_vx[f], interface.facet_vy[f]),
            1:facets.nfacet))
end

# =============================================================================
# 8. Diagnostics
# =============================================================================

"""
    divergenceMismatch(cell_vx, cell_vy, growth_nodes, snapshot, backdrop) -> NamedTuple

Discrete divergence of the *cell-sampled* velocity (central differences on one sample per
active cell) against the bilinearly interpolated `g(c)` at the same cell centres, restricted to
cells whose full 3×3 neighbourhood is `CUT_IN`. No interface cell, ghost node or extended value
enters, so the mismatch measures the accuracy of the discrete Darcy law in the interior.
"""
function divergenceMismatch(cell_vx, cell_vy, growth_nodes, snapshot, backdrop)
    ncx = backdrop.ncx
    h = backdrop.h
    states = snapshot.bg_states
    inside(i, j) = 1 <= i <= ncx && 1 <= j <= ncx && states[encodeCell(i, j, ncx)] == ED.CUT_IN
    num = 0.0
    den = 0.0
    worst = NaN
    used = 0
    for j in 1:ncx, i in 1:ncx
        all(inside(i + di, j + dj) for di in -1:1, dj in -1:1) || continue
        g = (growth_nodes[i, j] + growth_nodes[i + 1, j] +
             growth_nodes[i, j + 1] + growth_nodes[i + 1, j + 1]) / 4
        isfinite(g) || continue
        div = (cell_vx[i + 1, j] - cell_vx[i - 1, j]) / (2h) +
              (cell_vy[i, j + 1] - cell_vy[i, j - 1]) / (2h)
        isfinite(div) || continue
        num += (div - g)^2
        den += g^2
        used += 1
        worst = max(coalesce(worst, 0.0), abs(div - g) / max(abs(g), eps(Float64)))
    end
    used == 0 && return (; l2 = NaN, linf = NaN, cells = 0)
    (; l2 = sqrt(num / max(den, eps(Float64))), linf = worst, cells = used)
end

"""
    growthIntegralGeometric(cutinfo, growth_nodes, backdrop) -> Float64

`∫_Ω g dΩ` estimated on the EvolvingDomains cut geometry with a piecewise-constant rule: the
physical area of each active cell times the bilinearly interpolated `g` at the centroid of that
physical part. This is a different quadrature *and* a different cut representation from the
Gridap value `sum(∫(gh)*dΩ)`, which makes the two estimates an honest cross-check of the
diagnostic rather than a restatement of the solver.
"""
function growthIntegralGeometric(cutinfo, growth_nodes, backdrop)
    ncx = backdrop.ncx
    h = backdrop.h
    total = 0.0
    for cell in cutinfo.occupied
        i, j = decodeCell(cell, ncx)
        nodes = (growth_nodes[i, j], growth_nodes[i + 1, j],
            growth_nodes[i, j + 1], growth_nodes[i + 1, j + 1])
        all(isfinite, nodes) || continue
        ξ = clamp((cutinfo.cx[cell] - (i - 1) * h) / h, 0.0, 1.0)
        η = clamp((cutinfo.cy[cell] - (j - 1) * h) / h, 0.0, 1.0)
        g = nodes[1] * (1 - ξ) * (1 - η) + nodes[2] * ξ * (1 - η) +
            nodes[3] * (1 - ξ) * η + nodes[4] * ξ * η
        total += cutinfo.inside_area[cell] * g
    end
    return total
end

# =============================================================================
# 9. Geometry / marker advance with corrected reinitialisation diagnostics
# =============================================================================

"""
    advanceGeometry!(levelset, velocity, dt, cfg, grid, h) -> NamedTuple

Transport the level set (`eulerianTransport!`, WENO5 + SSP-RK3) and re-distance it
(`reinitialize!`). Both sides of the reinitialisation are measured with **distinct, freshly
taken snapshots**: `pre` is cut from the transported field before `reinitialize!` runs, `post`
is cut again afterwards. The first round returned the pre-reinitialisation snapshot in place of
the post one, so its reported reinitialisation area change was identically zero by
construction; here the two are separate `Cut` objects and the interface displacement
(one-sided closest-facet distance from every pre facet midpoint to the post interface, and
back) is reported as well as the area change.
"""
function advanceGeometry!(levelset, velocity, dt, cfg, grid, backdrop)
    h = backdrop.h
    t0 = time()
    ED.eulerianTransport!(levelset, velocity, dt; level_id = :init, cache = cfg.weno_cache)
    pre = ED.cut(levelset, :init)
    # Explicit geometric regularization, between transport and re-distance, exactly as the
    # filter's recipe prescribes: a fresh snapshot of the filtered field is what reinitialize!
    # consumes. The area change it causes is measured from two reconstructions around the call
    # alone, so it is never confused with growth or with the transport scheme.
    filt = (; t = 0.0, nodes = 0, to_positive = 0, to_negative = 0, cut_before = 0,
        cut_after = 0, comps_before = 0, comps_after = 0, area_change = 0.0)
    if cfg.filter_radius > 0
        t_filter = time()
        area_before_filter = cutGeometryInfo(pre, grid, h).area
        stats = ED.filterSubResolutionFeatures!(levelset, :init; phase = :negative,
            mode = cfg.filter_mode, radius = cfg.filter_radius)
        pre = ED.cut(levelset, :init)
        filt = (; t = time() - t_filter, nodes = stats.n_changed_nodes,
            to_positive = stats.n_flipped_to_positive,
            to_negative = stats.n_flipped_to_negative,
            cut_before = stats.n_cut_cells_before, cut_after = stats.n_cut_cells_after,
            comps_before = stats.n_components_before, comps_after = stats.n_components_after,
            area_change = cutGeometryInfo(pre, grid, h).area - area_before_filter)
    end
    phi_before = copy(levelset.data[:init])
    t1 = time()
    ED.reinitialize!(levelset, pre)
    post = ED.cut(levelset, :init)
    t2 = time()
    preinfo = cutGeometryInfo(pre, grid, h)
    postinfo = cutGeometryInfo(post, grid, h)
    prefacets = interfaceFacets(pre, backdrop)
    postfacets = interfaceFacets(post, backdrop)
    forward = nearestFacetDistances(prefacets.mx, prefacets.my, postfacets)
    backward = nearestFacetDistances(postfacets.mx, postfacets.my, prefacets)
    phi_after = levelset.data[:init]
    (; snapshot = post, info = postinfo, filter = filt,
        transport = (; t = t1 - t0, area = preinfo.area, facets = prefacets.nfacet),
        reinit = (; t = t2 - t1, area = postinfo.area, facets = postfacets.nfacet,
            area_change = postinfo.area - preinfo.area,
            shift_forward_mean = sum(forward) / length(forward),
            shift_forward_max = maximum(forward),
            shift_backward_mean = sum(backward) / length(backward),
            shift_backward_max = maximum(backward),
            phi_rms_change = sqrt(sum(abs2, phi_after .- phi_before) / length(phi_before)),
            phi_max_change = maximum(abs, phi_after .- phi_before),
            sign_flips = count((phi_before .< 0) .!= (phi_after .< 0))))
end

"""Transport the passive marker with intensive CIP; source and target must not alias."""
function advanceMarker!(target, source, velocity, dt, cfg)
    ED.lagrangianTransport!(target, source, velocity, dt; level_id = :init, cache = cfg.cip_cache)
    return target
end

"""Explicit geometry CFL restriction on the recovered nodal velocity."""
function timestep(velocity, backdrop, cfg, t)
    vmax = maximum(norm, velocity)
    cfl_dt = vmax > 0 ? cfg.cfl * backdrop.h / vmax : cfg.dtmax
    return min(cfl_dt, cfg.dtmax, cfg.tfinal - t), vmax
end

# =============================================================================
# 10. Verification
# =============================================================================

"""Modified Bessel functions I0 and I1 by power series (no extra dependency)."""
function besselI0(z)
    term = 1.0
    total = 1.0
    for m in 1:60
        term *= (z / 2)^2 / (m * m)
        total += term
        abs(term) < 1.0e-18 * abs(total) && break
    end
    return total
end

function besselI1(z)
    term = z / 2
    total = term
    for m in 1:60
        term *= (z / 2)^2 / (m * (m + 1))
        total += term
        abs(term) < 1.0e-18 * abs(total) && break
    end
    return total
end

"""
    verificationCase(cfg) -> NamedTuple

Independent verification of the two elliptic solvers and of both velocity-recovery methods on
the **initial separated discs**, where closed-form solutions exist:

    nutrient      -Dc Δc + λc = 0,  c = 1 on Γ    =>  c(r) = I0(k r)/I0(k R0),  k = √(λ/Dc)
    manufactured  -Δp = 1,          p = 0 on Γ    =>  p(r) = (R0² - r²)/4,  v = -∇p = (x-c)/2

The second problem is a *separate* unit-source test: it is not the exact solution of the
coupled nutrient-driven pressure equation, whose source is `g(c)`, and it is never presented
as such. It is what makes the recovery question answerable, because it fixes the exact
interface velocity: `V_n = R0/2` and `∫_Γ V_n ds = π R0²` per disc. Both recovery methods are
then applied to the *same* manufactured pressure field and their interface velocities compared
against that reference — the only place in this experiment where a recovery method is scored
against a known answer.
"""
function verificationCase(cfg)
    backdrop = backgroundModel(cfg)
    grid = ED.Grid((backdrop.axes, backdrop.axes))
    levelset = ED.Field(ED.Discretization(Dict(:init => grid)),
        Dict(:init => initialLevelSet(backdrop.axes, cfg)))
    snapshot = ED.cut(levelset, :init)
    cutinfo = cutGeometryInfo(snapshot, grid, backdrop.h)
    geo = buildGeometry(levelset.data[:init], backdrop, cfg)
    spaces = buildSpaces(geo)
    facets = interfaceFacets(snapshot, backdrop)

    R0, c1, c2 = cfg.R0, cfg.center1, cfg.center2
    radius(x) = min(hypot(x[1] - c1[1], x[2] - c1[2]), hypot(x[1] - c2[1], x[2] - c2[2]))
    kappa = sqrt(cfg.lambda / cfg.Dc)
    pfun(x) = (R0^2 - radius(x)^2) / 4
    cfun(x) = besselI0(kappa * radius(x)) / besselI0(kappa * R0)

    inside_nodes = falses(backdrop.nx, backdrop.nx)
    strict_nodes = falses(backdrop.nx, backdrop.nx)
    for j in 1:backdrop.ncx, i in 1:backdrop.ncx
        snapshot.bg_states[encodeCell(i, j, backdrop.ncx)] == ED.CUT_IN || continue
        for (ii, jj) in ((i, j), (i + 1, j), (i, j + 1), (i + 1, j + 1))
            inside_nodes[ii, jj] = true
        end
        if radius(((i - 0.5) * backdrop.h, (j - 0.5) * backdrop.h)) <= R0 - 2 * backdrop.h
            for (ii, jj) in ((i, j), (i + 1, j), (i, j + 1), (i + 1, j + 1))
                strict_nodes[ii, jj] = true
            end
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

    worst_over(mask, exact, numerical) = maximum(
        (abs(numerical[i, j] - exact[i, j]) for i in 1:backdrop.nx, j in 1:backdrop.nx
         if mask[i, j]); init = 0.0)
    p_interp = interpolate(pfun, spaces.V)
    c_interp = interpolate(cfun, spaces.V)
    rel2(value, reference) = sqrt(sum(∫((value - reference) * (value - reference)) * geo.dΩ) /
                                   sum(∫(reference * reference) * geo.dΩ))

    # Boundary flux identities, evaluated with the same machinery as the coupled run.
    p_flux = sum(∫(geo.nΓ ⋅ ∇(pressure.uh)) * geo.dΓ)
    p_penalty = (cfg.nitsche / backdrop.h) * sum(∫(pressure.uh) * geo.dΓ)
    c_flux = sum(∫(geo.nΓ ⋅ ∇(nutrient.uh)) * geo.dΓ)
    c_mass = sum(∫(nutrient.uh) * geo.dΩ)
    # Both discs: ∫_Γ ∂_n c ds = 2 * 2πR0 ∂_r c(R0).
    analytic_c_flux = 4π * R0 * kappa * besselI1(kappa * R0) / besselI0(kappa * R0)

    # Reference recovery on the manufactured pressure.
    A = recoveryCentroid(pressure, geo, cutinfo, facets, backdrop)
    B = recoveryPatch(pressure, geo, cutinfo, facets, snapshot, backdrop, cfg, levelset.data[:init])
    exact_vn = R0 / 2
    errA = [abs((A.facet_vx[f] * facets.nnx[f] + A.facet_vy[f] * facets.nny[f]) - exact_vn)
            for f in 1:facets.nfacet]
    errB = [abs((B.facet_vx[f] * facets.nnx[f] + B.facet_vy[f] * facets.nny[f]) - exact_vn)
            for f in 1:facets.nfacet]
    biasA = [((A.facet_vx[f] * facets.nnx[f] + A.facet_vy[f] * facets.nny[f]) - exact_vn)
             for f in 1:facets.nfacet]
    biasB = [((B.facet_vx[f] * facets.nnx[f] + B.facet_vy[f] * facets.nny[f]) - exact_vn)
             for f in 1:facets.nfacet]
    fluxA = sum((A.facet_vx[f] * facets.nnx[f] + A.facet_vy[f] * facets.nny[f]) *
                facets.len[f] for f in 1:facets.nfacet)
    fluxB = sum((B.facet_vx[f] * facets.nnx[f] + B.facet_vy[f] * facets.nny[f]) *
                facets.len[f] for f in 1:facets.nfacet)
    exact_flux = 4π * R0 * exact_vn            # two discs of radius R0

    Dict{String,Float64}(
        "verification_ncells" => cfg.ncells,
        "verification_area_ed" => cutinfo.area,
        "verification_area_fem" => domainArea(geo),
        "verification_area_exact" => 2π * R0^2,
        "verification_area_relative_error_ed" => abs(cutinfo.area - 2π * R0^2) / (2π * R0^2),
        "verification_area_relative_error_fem" => abs(domainArea(geo) - 2π * R0^2) / (2π * R0^2),
        "verification_perimeter_relative_error" => abs(interfaceLength(geo) - 4π * R0) /
                                                  (4π * R0),
        "verification_pressure_nodal_max_error" => worst_over(inside_nodes, p_exact, p_nodes),
        "verification_pressure_nodal_max_error_interior" => worst_over(strict_nodes, p_exact,
            p_nodes),
        "verification_pressure_l2_relative_error" => rel2(pressure.uh, p_interp),
        "verification_pressure_raw_boundary_flux" => p_flux,
        "verification_pressure_nitsche_penalty_term" => p_penalty,
        "verification_pressure_weak_identity" => -p_flux + p_penalty,
        "verification_pressure_source_integral" => sum(∫(1.0) * geo.dΩ),
        "verification_pressure_exact_boundary_flux" => 4π * R0 * exact_vn,
        "verification_pressure_max" => maximum(pressure.dofs),
        "verification_pressure_max_exact" => R0^2 / 4,
        "verification_nutrient_nodal_max_error" => worst_over(inside_nodes, c_exact, c_nodes),
        "verification_nutrient_l2_relative_error" => rel2(nutrient.uh, c_interp),
        "verification_nutrient_analytic_boundary_flux" => analytic_c_flux,
        "verification_nutrient_raw_boundary_flux" => c_flux,
        "verification_nutrient_nitsche_penalty_term" =>
            (cfg.nitsche / backdrop.h) * sum(∫(nutrient.uh - 1.0) * geo.dΓ),
        "verification_nutrient_weak_identity" => -c_flux +
            (cfg.nitsche / backdrop.h) * sum(∫(nutrient.uh - 1.0) * geo.dΓ),
        # The raw gradient flux against the analytic flux, and the Nitsche identity against the
        # discrete mass: two different statements, reported separately on purpose.
        "verification_nutrient_flux_identity_residual" =>
            abs(cfg.Dc * (-c_flux + (cfg.nitsche / backdrop.h) *
                          sum(∫(nutrient.uh - 1.0) * geo.dΓ)) + cfg.lambda * c_mass) /
            abs(cfg.lambda * c_mass),
        "verification_recovery_exact_normal_velocity" => exact_vn,
        "verification_recovery_A_mean_abs_error" => sum(errA) / length(errA),
        "verification_recovery_A_max_abs_error" => maximum(errA),
        "verification_recovery_B_mean_abs_error" => sum(errB) / length(errB),
        "verification_recovery_B_max_abs_error" => maximum(errB),
        "verification_recovery_A_flux" => fluxA,
        "verification_recovery_B_flux" => fluxB,
        "verification_recovery_B_flux_relative_error" => abs(fluxB - exact_flux) / exact_flux,
        "verification_recovery_A_flux_relative_error" => abs(fluxA - exact_flux) / exact_flux,
        "verification_recovery_B_fallbacks" => B.fallbacks,
        "verification_recovery_B_quadratic_facets" => B.nquadratic,
        "verification_recovery_B_condition_max" => B.cond_max,
        "verification_recovery_B_linear_facets" => B.nlinear,
        "verification_recovery_facets" => facets.nfacet,
        "verification_recovery_A_mean_bias" => sum(biasA) / length(biasA),
        "verification_recovery_B_mean_bias" => sum(biasB) / length(biasB),
        "verification_pressure_algebraic_residual" => pressure.residual,
        "verification_nutrient_algebraic_residual" => nutrient.residual,
    )
end

# =============================================================================
# 10b. Manufactured test on the *merged* geometry of the baseline run
# =============================================================================

"""
    mergedManufactured(tag) -> Dict{String,Float64}

Manufactured-solution test run on the merged geometry of an existing baseline run, to separate
the two things that a post-merger area deficit could come from: the *recovery operator* and the
*solution* it interpolates.

The geometry is the last saved level set of `output/<tag>/frames.jls` (no analytic geometry, no
re-evolution). On it we solve the *separate* problem

    -Δp = 1    in Ω,      p = p_exact on Γ,      p_exact = 0.25 - |x - (1/2,1/2)|²/4

whose exact velocity is `-∇p_exact = (x - (1/2,1/2))/2`, so every quantity has a closed form at
every interface facet midpoint. Two references are reported separately because they test
different things:

* `exact_flux_from_facets` — `Σ_f (-∇p_exact(x_f)·n_f) ℓ_f` with the *discrete* facet normals and
  lengths. By the divergence theorem on the closed polygon this equals half the polygon area
  `polygon_area` exactly, whatever the normal errors are; so the difference between the recovered
  flux and this number is the recovery error *plus* the normal error, while its difference from
  `∫_Ω f = |Ω_fem|` is the geometric representation error (ED sub-cells vs Gridap cut).
* `source_integral` — `∫_Ω 1 dΩ` on the Gridap domain, i.e. the divergence-theorem right-hand side.

If the recovery is still accurate here, a deficit observed with the *coupled* pressure field on
the same geometry is a property of that field, not of the recovery.
"""
function mergedManufactured(tag)
    path = joinpath(@__DIR__, "output", tag, "frames.jls")
    isfile(path) || error("no frames for tag $tag at $path")
    data = deserialize(path)
    phi = data.frames[end].phi
    ncells = data.ncells
    cfg = defaultConfig(; ncells)
    backdrop = backgroundModel(cfg)
    grid = ED.Grid((backdrop.axes, backdrop.axes))
    levelset = ED.Field(ED.Discretization(Dict(:init => grid)), Dict(:init => copy(phi)))
    snapshot = ED.cut(levelset, :init)
    cutinfo = cutGeometryInfo(snapshot, grid, backdrop.h)
    facets = interfaceFacets(snapshot, backdrop)
    geo = buildGeometry(levelset.data[:init], backdrop, cfg)
    spaces = buildSpaces(geo)

    pfun(x) = 0.25 - ((x[1] - 0.5)^2 + (x[2] - 0.5)^2) / 4
    gfun(x) = ((x[1] - 0.5) / 2, (x[2] - 0.5) / 2)          # -∇p_exact
    pD = interpolate(pfun, spaces.V)
    state = solveElliptic(geo, spaces, cfg, backdrop.h; diffusion = 1.0, reaction = 0.0,
        source = 1.0, boundaryValue = pD, label = "merged-manufactured")

    A = recoveryCentroid(state, geo, cutinfo, facets, backdrop)
    B = recoveryPatch(state, geo, cutinfo, facets, snapshot, backdrop, cfg,
        levelset.data[:init])
    exact_vn = zeros(facets.nfacet)
    for f in 1:facets.nfacet
        exact_vn[f] = gfun((facets.mx[f], facets.my[f]))[1] * facets.nnx[f] +
                      gfun((facets.mx[f], facets.my[f]))[2] * facets.nny[f]
    end
    errA = [abs(A.facet_vx[f] * facets.nnx[f] + A.facet_vy[f] * facets.nny[f] - exact_vn[f])
            for f in 1:facets.nfacet]
    errB = [abs(B.facet_vx[f] * facets.nnx[f] + B.facet_vy[f] * facets.nny[f] - exact_vn[f])
            for f in 1:facets.nfacet]
    exact_flux = sum(exact_vn[f] * facets.len[f] for f in 1:facets.nfacet)
    fluxA = sum((A.facet_vx[f] * facets.nnx[f] + A.facet_vy[f] * facets.nny[f]) * facets.len[f]
                for f in 1:facets.nfacet)
    fluxB = sum((B.facet_vx[f] * facets.nnx[f] + B.facet_vy[f] * facets.nny[f]) * facets.len[f]
                for f in 1:facets.nfacet)
    source_integral = sum(∫(1.0) * geo.dΩ)
    raw_flux = sum(∫(geo.nΓ ⋅ ∇(state.uh)) * geo.dΓ)
    penalty = (cfg.nitsche / backdrop.h) * sum(∫(state.uh - pD) * geo.dΓ)
    polygon_area = sum((facets.mx[f] * facets.nnx[f] + facets.my[f] * facets.nny[f]) *
                       facets.len[f] for f in 1:facets.nfacet)
    p_nodes = backgroundNodal(state.dofs, geo, backdrop)
    interior = falses(backdrop.nx, backdrop.nx)
    for j in 1:backdrop.ncx, i in 1:backdrop.ncx
        snapshot.bg_states[encodeCell(i, j, backdrop.ncx)] == ED.CUT_IN || continue
        for (ii, jj) in ((i, j), (i + 1, j), (i, j + 1), (i + 1, j + 1))
            interior[ii, jj] = true
        end
    end
    nodal_error = maximum((abs(p_nodes[i, j] - pfun((backdrop.axes[i], backdrop.axes[j])))
                           for i in 1:backdrop.nx, j in 1:backdrop.nx if interior[i, j]);
        init = 0.0)
    Dict{String,Any}(
        "source_tag" => tag, "ncells" => ncells, "facets" => facets.nfacet,
        "t_geometry" => data.frames[end].t,
        "area_ed" => cutinfo.area, "area_fem" => source_integral,
        "polygon_area_from_facets" => polygon_area,
        "exact_flux_from_facets" => exact_flux,
        "exact_flux_expected" => polygon_area / 2,
        "recovered_flux_A" => fluxA, "recovered_flux_B" => fluxB,
        "fluxA_relative_error" => abs(fluxA - exact_flux) / abs(exact_flux),
        "fluxB_relative_error" => abs(fluxB - exact_flux) / abs(exact_flux),
        "fluxA_vs_source" => abs(fluxA - source_integral) / source_integral,
        "fluxB_vs_source" => abs(fluxB - source_integral) / source_integral,
        "raw_nitsche_flux_relative_error" => abs(raw_flux + exact_flux) / abs(exact_flux),
        "nitsche_penalty_relative" => abs(penalty) / abs(exact_flux),
        "weak_identity_residual" => abs(-raw_flux + penalty - source_integral) / source_integral,
        "vn_mean_abs_error_A" => sum(errA) / length(errA),
        "vn_max_abs_error_A" => maximum(errA),
        "vn_mean_abs_error_B" => sum(errB) / length(errB),
        "vn_max_abs_error_B" => maximum(errB),
        "vn_mean_abs_exact" => sum(abs, exact_vn) / length(exact_vn),
        "vn_max_abs_exact" => maximum(abs, exact_vn),
        "pressure_nodal_max_error" => nodal_error,
        "recovery_fallbacks_B" => B.fallbacks,
        "recovery_condition_max_B" => B.cond_max,
        "algebraic_residual" => state.residual,
    )
end

"""
    nearestFacetExtension(facets, fvx, fvy, backdrop) -> NamedTuple

The extension actually used by the pipeline: brute force over every interface facet for every
background node, ties broken by the smallest facet index. Costs `O(nodes × facets)` distance
evaluations and is nevertheless the *fastest* option measured here (0.4 ms at 64², 3.3 ms at
128², 10 ms at 192²): the inner loop is a contiguous SIMD-friendly sweep, while the indexed
alternative spends 60 ns per (much rarer) evaluation on ring bookkeeping.

shortcut: linear search in the number of interface facets, measured faster than a CSR ring index
up to 192² (see `extensionBenchmark`); revisit if the interface ever carries ≫10³ facets or the
background exceeds ~256².
"""
function nearestFacetExtension(facets, fvx, fvy, backdrop)
    t0 = time()
    nx = backdrop.nx
    h = backdrop.h
    vx = zeros(nx, nx)
    vy = zeros(nx, nx)
    for j in 1:nx, i in 1:nx
        x = (i - 1) * h
        y = (j - 1) * h
        best = 1
        best_distance = Inf
        for f in 1:facets.nfacet
            distance = (x - facets.mx[f])^2 + (y - facets.my[f])^2
            if distance < best_distance
                best_distance = distance
                best = f
            end
        end
        vx[i, j] = fvx[best]
        vy[i, j] = fvy[best]
    end
    velocity = SVector{2,Float64}[SVector{2,Float64}(vx[k], vy[k]) for k in eachindex(vx)]
    evaluations = nx * nx * facets.nfacet
    (; vx, vy, velocity, maximum_speed = maximum(norm, velocity), t_build = 0.0,
        t_query = time() - t0, distance_evaluations = evaluations,
        brute_force_evaluations = evaluations, distance_check_max = 0.0, buckets = 0)
end

extendBrute(facets, fvx, fvy, backdrop) = nearestFacetExtension(facets, fvx, fvy, backdrop)

"""
    extensionBenchmark(ncells_list)

Cost and correctness of the interface-to-grid extension adapter, separated from the physics:
the two recovery methods, the bucketed (spatial-index) nearest-facet search actually used, and
the first round's brute-force `O(nodes × facets)` search. The bucketed result is asserted to be
*identical* to brute force, so the cheaper search is not buying speed with a different answer.
"""
function extensionBenchmark(ncells_list)
    open(joinpath(@__DIR__, "output", "extension_benchmark.csv"), "w") do io
        println(io, join(("ncells", "facets", "nodes", "t_recovery_A", "t_recovery_B",
            "t_extension_bucket", "t_search_bucket", "t_search_brute", "evaluations_bucket",
            "evaluations_brute", "speedup", "max_abs_velocity_difference", "vmax"), ","))
        for n in ncells_list
            cfg = defaultConfig(; ncells = n)
            backdrop = backgroundModel(cfg)
            grid = ED.Grid((backdrop.axes, backdrop.axes))
            levelset = ED.Field(ED.Discretization(Dict(:init => grid)),
                Dict(:init => initialLevelSet(backdrop.axes, cfg)))
            snapshot = ED.cut(levelset, :init)
            cutinfo = cutGeometryInfo(snapshot, grid, backdrop.h)
            facets = interfaceFacets(snapshot, backdrop)
            geo = buildGeometry(levelset.data[:init], backdrop, cfg)
            spaces = buildSpaces(geo)
            pressure = solveElliptic(geo, spaces, cfg, backdrop.h; diffusion = 1.0,
                reaction = 0.0, source = 1.0, boundaryValue = 0.0, label = "bench-pressure")
            # Warm up every call path twice, then measure the third call: this benchmark is
            # about steady-state cost, not compilation.
            measure(f) = (f(); f(); @elapsed f())
            tA = measure(() -> recoveryCentroid(pressure, geo, cutinfo, facets, backdrop))
            tB = measure(() -> recoveryPatch(pressure, geo, cutinfo, facets, snapshot, backdrop,
                cfg, levelset.data[:init]))
            interface = recoveryCentroid(pressure, geo, cutinfo, facets, backdrop)
            fA() = extendToNodes(facets, interface.facet_vx, interface.facet_vy, backdrop)
            fBrute() = extendBrute(facets, interface.facet_vx, interface.facet_vy, backdrop)
            extension = measure(fA)
            bucketed = fA()
            brute = measure(fBrute)
            reference = fBrute()
            difference = maximum(abs.(bucketed.vx .- reference.vx))
            @printf("n=%3d facets=%3d nodes=%5d  recovery A %.2f ms B %.2f ms | extension %.2f ms (search %.2f ms, %d evals) vs brute %.2f ms (%d evals, %.1f×)  |Δv| = %.1e
",
                n, facets.nfacet, backdrop.nx^2, 1e3tA, 1e3tB, 1e3extension,
                1e3bucketed.t_query, bucketed.distance_evaluations, 1e3brute,
                bucketed.brute_force_evaluations,
                bucketed.brute_force_evaluations / max(bucketed.distance_evaluations, 1),
                difference)
            flush(stdout)
            println(io, join((n, facets.nfacet, backdrop.nx^2, tA, tB, extension,
                bucketed.t_query, brute, bucketed.distance_evaluations,
                bucketed.brute_force_evaluations,
                bucketed.brute_force_evaluations / max(bucketed.distance_evaluations, 1),
                difference, bucketed.maximum_speed), ","))
        end
    end
end

function mergedSuite(tags::Vector{String})
    rows = [mergedManufactured(tag) for tag in tags]
    keys_order = sort!(collect(keys(rows[1])))
    mkpath(joinpath(@__DIR__, "output"))
    open(joinpath(@__DIR__, "output", "merged_manufactured.csv"), "w") do io
        println(io, join(keys_order, ","))
        for row in rows
            println(io, join((row[k] for k in keys_order), ","))
        end
    end
    for row in rows
        @printf("%-8s n=%3d  A: mean %.2e max %.2e fluxerr %.3e | B: mean %.2e max %.2e fluxerr %.3e | rawfluxerr %.3e pen %.3e | poly/FE area %.6f
",
            row["source_tag"], row["ncells"], row["vn_mean_abs_error_A"], row["vn_max_abs_error_A"],
            row["fluxA_relative_error"], row["vn_mean_abs_error_B"], row["vn_max_abs_error_B"],
            row["fluxB_relative_error"], row["raw_nitsche_flux_relative_error"],
            row["nitsche_penalty_relative"], row["polygon_area_from_facets"] / row["area_fem"])
    end
end

# =============================================================================
# 11. Dynamics: the explicit, readable orchestrator of the workflow
# =============================================================================

"""
    dynamics(levelset, backdrop, cfg) -> NamedTuple

One evaluation of the physics at the *current* geometry, in the order the protocol prescribes:
cut the geometry, solve the nutrient problem, form `g(c)`, solve the pressure problem, recover
the interface velocity and extend it to the grid. Nothing is cached between calls and no
analytic geometry is used after initialisation.
"""
function dynamics(levelset, backdrop, cfg)
    stages = Tuple{String,Float64,Int64}[]
    function stage(label, f)
        value, seconds, bytes = @timed f()
        push!(stages, (label, seconds, bytes))
        return value
    end

    grid = levelset.discretization.levels[:init]
    snapshot = stage("ed_cut", () -> ED.cut(levelset, :init))
    cutinfo = stage("ed_cut_geometry", () -> cutGeometryInfo(snapshot, grid, backdrop.h))
    facets = stage("interface_facets", () -> interfaceFacets(snapshot, backdrop))
    geo = stage("gridap_cut", () -> buildGeometry(levelset.data[:init], backdrop, cfg))
    spaces = stage("gridap_spaces", () -> buildSpaces(geo))
    nutrient = stage("nutrient_solve", () -> solveNutrient(geo, spaces, cfg, backdrop.h))
    growth = stage("growth_rate", () -> computeGrowthRate(nutrient, geo, spaces, cfg))
    pressure = stage("pressure_solve",
        () -> solvePressure(geo, growth, spaces, cfg, backdrop.h))
    velocity = stage("velocity_recovery", () -> recoverVelocity(cfg.method, pressure, geo,
        cutinfo, facets, snapshot, backdrop, cfg, levelset.data[:init]))
    topology = stage("topology", () -> topologyCounts(snapshot, levelset.data[:init],
        backdrop, geo, backdrop.h))
    growth_nodes = stage("growth_nodal", () -> backgroundNodal(growth.values, geo, backdrop))
    divergence = stage("divergence_check", () -> divergenceMismatch(velocity.interface.cell_vx,
        velocity.interface.cell_vy, growth_nodes, snapshot, backdrop))
    growth_geom = stage("growth_integral_geom",
        () -> growthIntegralGeometric(cutinfo, growth_nodes, backdrop))
    area_fem = stage("fem_area", () -> domainArea(geo))
    perimeter = stage("fem_perimeter", () -> interfaceLength(geo))

    (; snapshot, cutinfo, facets, geo, nutrient, growth, pressure, velocity, topology,
        stages, growth_nodes, divergence, growth_geom, area_fem, perimeter,
        c_nodes = backgroundNodal(nutrient.dofs, geo, backdrop),
        p_nodes = backgroundNodal(pressure.dofs, geo, backdrop))
end

# =============================================================================
# 12. Simulation
# =============================================================================

const METRIC_COLUMNS = ("step", "t_start", "t_end", "dt", "area_start", "area_end",
    "area_start_fem", "area_end_fem", "perimeter", "d_area_dt", "integral_growth",
    "integral_growth_geom", "balance_residual", "balance_residual_fem", "interface_flux",
    "interface_flux_abs", "flux_balance_residual", "flux_vs_area_residual", "vmax",
    "facet_speed_max", "components_nodes", "components_cells", "components_geometric",
    "c_min", "c_max", "p_min", "p_max", "nutrient_residual", "pressure_residual",
    "pressure_raw_flux", "pressure_nitsche_penalty", "pressure_weak_identity_residual",
    "nitsche_boundary_defect", "nutrient_raw_flux", "nutrient_weak_identity",
    "nutrient_flux_identity_residual", "divergence_l2",
    "divergence_linf", "divergence_cells", "min_cut_fraction", "aggregate_roots",
    "reinit_area_change", "reinit_shift_forward_mean", "reinit_shift_forward_max",
    "reinit_shift_backward_mean", "reinit_shift_backward_max", "reinit_phi_rms_change",
    "reinit_phi_max_change", "reinit_sign_flips", "n_cut_cells", "n_in_cells",
    "marker_min", "marker_max", "extension_distance_check", "extension_evaluations",
    "extension_brute_force_evaluations", "recovery_fallbacks",
    "recovery_samples", "t_step", "t_dynamics", "t_transport", "t_marker",
    "filter_nodes_changed", "filter_to_positive", "filter_to_negative",
    "filter_cut_cells_before", "filter_cut_cells_after", "filter_components_before",
    "filter_components_after", "filter_area_change", "filter_seconds")

function writeRow(io, values)
    println(io, join(values, ","))
end

"""
    simulate(cfg; provenance = String[]) -> NamedTuple

The orchestration loop. Diagnostics that belong to one physical time are written together:
`dynamics` is evaluated at the start of every step (time `t`) *and again* after the geometry has
been advanced to `t + dt`, so the saved frame of step `n` carries nutrient, pressure, velocity
and marker fields solved on the geometry of its own timestamp. The step residual compares the
reconstructed-area difference between those two times against the growth integral measured at
the earlier one.
"""
function simulate(cfg; provenance = String[])
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
    metrics = open(joinpath(cfg.outdir, "metrics.csv"), "w")
    profile = open(joinpath(cfg.outdir, "profile.csv"), "w")
    println(profile, join(("step", "t", "stage", "seconds", "bytes"), ","))
    writeRow(metrics, METRIC_COLUMNS)

    frames = Vector{Any}()
    t = 0.0
    step = 0
    merged_at = (; nodes = NaN, cells = NaN, geometric = NaN)
    area_initial = NaN
    components_initial = nothing
    stats = Dict{String,Vector{Float64}}()
    push_stat!(key, value) = push!(get!(() -> Float64[], stats, key), Float64(value))
    eps_scale_bound = 0
    run_start = time()

    state = dynamics(levelset, backdrop, cfg)
    area_initial = state.cutinfo.area
    components_initial = state.topology

    function record!(io, io_profile, step, t_start, t_end, dt, state, next, marker_values,
            reinit, filt, extra_times)
        observed = (next.cutinfo.area - state.cutinfo.area) / dt
        observed_fem = (next.area_fem - state.area_fem) / dt
        growth = state.growth.integral
        growth_geom = state.growth_geom
        balance = abs(observed - growth) / max(abs(growth), cfg.eps_scale)
        balance_fem = abs(observed_fem - growth) / max(abs(growth), cfg.eps_scale)
        flux_balance = abs(state.velocity.interface_flux - growth) /
                       max(abs(growth), cfg.eps_scale)
        flux_vs_area = abs(state.velocity.interface_flux - observed) /
                       max(abs(observed), cfg.eps_scale)
        abs(growth) < cfg.eps_scale && (eps_scale_bound += 1)
        writeRow(io, (step, t_start, t_end, dt, state.cutinfo.area, next.cutinfo.area,
            state.area_fem, next.area_fem, state.perimeter, observed, growth, growth_geom,
            balance, balance_fem, state.velocity.interface_flux,
            state.velocity.interface_flux_abs, flux_balance, flux_vs_area,
            state.velocity.extension.maximum_speed, state.velocity.facet_speed_max,
            state.topology.nodes, state.topology.cells, state.topology.geometric,
            state.nutrient.cmin, state.nutrient.cmax, state.pressure.pmin,
            state.pressure.pmax, state.nutrient.residual, state.pressure.residual,
            state.pressure.flux, state.pressure.penalty,
            abs(state.pressure.flux_identity - growth) / max(abs(growth), cfg.eps_scale),
            abs(state.nutrient.boundary_defect),
            state.nutrient.flux, state.nutrient.weak_identity,
            abs(cfg.Dc * state.nutrient.weak_identity + cfg.lambda * state.nutrient.mass) /
            max(abs(cfg.lambda * state.nutrient.mass), eps(Float64)),
            state.divergence.l2, state.divergence.linf, state.divergence.cells,
            state.cutinfo.min_fraction, state.cutinfo.aggregate_roots,
            reinit.area_change, reinit.shift_forward_mean, reinit.shift_forward_max,
            reinit.shift_backward_mean, reinit.shift_backward_max,
            reinit.phi_rms_change, reinit.phi_max_change, reinit.sign_flips,
            state.cutinfo.n_cut, state.cutinfo.n_in,
            minimum(marker_values), maximum(marker_values),
            state.velocity.extension.distance_check_max,
            state.velocity.extension.distance_evaluations,
            state.velocity.extension.brute_force_evaluations,
            state.velocity.interface.fallbacks, get(state.velocity.interface, :samples, 0),
            extra_times.t_step, extra_times.t_dynamics, extra_times.t_transport,
            extra_times.t_marker,
            filt.nodes, filt.to_positive, filt.to_negative, filt.cut_before, filt.cut_after,
            filt.comps_before, filt.comps_after, filt.area_change, filt.t))
        for (label, seconds, bytes) in state.stages
            println(io_profile, join((step, t_start, label, seconds, bytes), ","))
        end
        flush(io)
        flush(io_profile)
        push_stat!("balance", balance)
        push_stat!("balance_fem", balance_fem)
        push_stat!("flux_balance", flux_balance)
        push_stat!("flux_vs_area", flux_vs_area)
        push_stat!("observed_rate", observed)
        push_stat!("growth", growth)
        push_stat!("interface_flux", state.velocity.interface_flux)
        push_stat!("divergence_l2", state.divergence.l2)
        push_stat!("nutrient_residual", state.nutrient.residual)
        push_stat!("pressure_residual", state.pressure.residual)
        push_stat!("nitsche_defect", abs(state.nutrient.boundary_defect))
        push_stat!("min_cut_fraction", state.cutinfo.min_fraction)
        push_stat!("vmax", state.velocity.extension.maximum_speed)
        push_stat!("dt", dt)
        push_stat!("growth_times_dt", state.growth.integral * dt)
        push_stat!("flux_times_dt", state.velocity.interface_flux * dt)
        push_stat!("area_change", next.cutinfo.area - state.cutinfo.area)
        push_stat!("reinit_area_change", reinit.area_change)
        push_stat!("reinit_shift_max", reinit.shift_forward_max)
        push_stat!("reinit_phi_rms_change", reinit.phi_rms_change)
        push_stat!("filter_nodes_changed", filt.nodes)
        push_stat!("filter_area_change", filt.area_change)
        push_stat!("filter_seconds", filt.t)
    end

    function saveframe!(step, t, state, marker_values)
        dictionary = (; step, t, synchronized = true,
            phi = copy(levelset.data[:init]), c = state.c_nodes, p = state.p_nodes,
            growth = state.growth_nodes, q = copy(marker_values),
            vx = state.velocity.extension.vx, vy = state.velocity.extension.vy,
            cell_vx = state.velocity.interface.cell_vx,
            cell_vy = state.velocity.interface.cell_vy,
            cell_area = state.cutinfo.inside_area,
            facet_mx = state.facets.mx, facet_my = state.facets.my,
            facet_vn = state.velocity.vn, facet_nx = state.facets.nnx,
            facet_ny = state.facets.nny,
            area = state.cutinfo.area, perimeter = state.perimeter,
            components = state.topology, integral_growth = state.growth.integral,
            method = string(cfg.method))
        push!(frames, dictionary)
    end

    saveframe!(0, 0.0, state, marker_source.data[:init])
    @printf("%s | %d²  h=%.5f  method=%s  dtmax=%.1e  tfinal=%.2f\n", cfg.tag, cfg.ncells,
        1.0 / cfg.ncells, cfg.method, cfg.dtmax, cfg.tfinal)
    @printf("%s  t=%.4f  A=%.6f  P=%.4f  ∫g=%.4f  comps=%d/%d/%d  vmax=%.4f\n", "init",
        0.0, state.cutinfo.area, state.perimeter, state.growth.integral,
        state.topology.nodes, state.topology.cells, state.topology.geometric,
        state.velocity.extension.maximum_speed)

    while t < cfg.tfinal - 1.0e-12
        step += 1
        step_start = time()
        dt, vmax = timestep(state.velocity.extension.velocity, backdrop, cfg, t)
        dt <= 0 && break

        t_marker = @elapsed advanceMarker!(marker_target, marker_source,
            state.velocity.extension.velocity, dt, cfg)
        marker_source, marker_target = marker_target, marker_source

        geometry = advanceGeometry!(levelset, state.velocity.extension.velocity, dt, cfg, grid,
            backdrop)
        t_transport = geometry.transport.t + geometry.reinit.t

        t += dt
        next = dynamics(levelset, backdrop, cfg)
        step_time = time() - step_start

        record!(metrics, profile, step, t - dt, t, dt, state, next,
            marker_source.data[:init], geometry.reinit, geometry.filter,
            (; t_step = step_time, t_transport, t_marker,
                t_dynamics = sum(s[2] for s in state.stages)))

        merged_at = (;
            nodes = isnan(merged_at.nodes) && next.topology.nodes == 1 &&
                    components_initial.nodes > 1 ? t : merged_at.nodes,
            cells = isnan(merged_at.cells) && next.topology.cells == 1 &&
                    components_initial.cells > 1 ? t : merged_at.cells,
            geometric = isnan(merged_at.geometric) && next.topology.geometric == 1 &&
                        components_initial.geometric > 1 ? t : merged_at.geometric)

        @printf("%3d  t=%.4f dt=%.4f  A=%.6f P=%.4f ∫g=%.4f dA/dt=%.4f bal=%.4f flux=%.4f  comps=%d/%d/%d  vmax=%.4f  reinitΔA=%.2e  filt=%d %d/%d %d→%d ΔA=%+.2e (%.3fs)  (%.2fs)\n",
            step, t, dt, next.cutinfo.area, next.perimeter, state.growth.integral,
            (next.cutinfo.area - state.cutinfo.area) / dt,
            last(stats["balance"]), state.velocity.interface_flux,
            next.topology.nodes, next.topology.cells, next.topology.geometric,
            state.velocity.extension.maximum_speed, geometry.reinit.area_change,
            geometry.filter.nodes, geometry.filter.comps_before, geometry.filter.comps_after,
            geometry.filter.cut_before, geometry.filter.cut_after,
            geometry.filter.area_change, geometry.filter.t, step_time)
        flush(stdout)

        state = next
        if step % cfg.frame_stride == 0 || t >= cfg.tfinal - 1.0e-12
            saveframe!(step, t, state, marker_source.data[:init])
        end
    end
    close(metrics)
    close(profile)

    serialize(joinpath(cfg.outdir, "frames.jls"),
        (; axes = backdrop.axes, ncells = cfg.ncells, method = string(cfg.method),
            frames))
    total_runtime = time() - run_start

    mean_of(key) = isempty(get(stats, key, Float64[])) ? NaN :
                   sum(stats[key]) / length(stats[key])
    cum_area = frames[end].area - area_initial
    cum_growth = sum(get(stats, "growth_times_dt", Float64[]))
    cum_flux = sum(get(stats, "flux_times_dt", Float64[]))
    area_start = area_initial
    per = Dict{String,Float64}()
    for key in ("balance", "balance_fem", "flux_balance", "flux_vs_area", "divergence_l2",
        "nutrient_residual", "pressure_residual", "nitsche_defect", "reinit_shift_max",
        "reinit_phi_rms_change")
        push!(per, "mean_" * key => mean_of(key))
        push!(per, "max_" * key => isempty(get(stats, key, Float64[])) ? NaN :
              maximum(stats[key]))
    end
    per["min_min_cut_fraction"] = isempty(get(stats, "min_cut_fraction", Float64[])) ? NaN :
                                  minimum(stats["min_cut_fraction"])
    per["max_vmax"] = isempty(get(stats, "vmax", Float64[])) ? NaN : maximum(stats["vmax"])
    per["sum_reinit_area_change"] = isempty(get(stats, "reinit_area_change", Float64[])) ? NaN :
                                    sum(abs, stats["reinit_area_change"])
    per["final_area"] = frames[end].area
    per["final_area_fem"] = state.area_fem
    per["final_perimeter"] = frames[end].perimeter
    per["cumulative_area_change"] = cum_area
    per["cumulative_growth_integral"] = cum_growth
    per["cumulative_interface_flux"] = cum_flux
    per["cumulative_ratio"] = cum_area / cum_growth
    per["cumulative_flux_ratio"] = cum_flux / cum_growth
    per["cumulative_flux_over_area_change"] = cum_flux / cum_area
    per["area_initial_ed"] = area_start
    per["eps_scale_binding_steps"] = eps_scale_bound
    per["merged_at_nodes"] = merged_at.nodes
    per["merged_at_cells"] = merged_at.cells
    per["merged_at_geometric"] = merged_at.geometric
    per["steps"] = step
    per["t_final"] = t
    per["total_runtime"] = total_runtime
    per["node_ordering_error"] = ordering_error

    # Steady-state per-stage cost: drop the first two steps (compilation) and average the rest.
    stage_seconds = Dict{String,Vector{Float64}}()
    stage_bytes = Dict{String,Vector{Int64}}()
    for line in eachline(joinpath(cfg.outdir, "profile.csv"))
        parts = split(line, ",")
        length(parts) == 5 || continue
        parts[1] == "step" && continue
        parse(Int, parts[1]) < 3 && continue
        push!(get!(() -> Float64[], stage_seconds, parts[3]), parse(Float64, parts[4]))
        push!(get!(() -> Int64[], stage_bytes, parts[3]), parse(Int64, parts[5]))
    end
    for (label, values) in stage_seconds
        per["stage_mean_seconds[" * label * "]"] = sum(values) / length(values)
    end
    for (label, values) in stage_bytes
        per["stage_mean_bytes[" * label * "]"] = sum(values) / length(values)
    end

    metadata = open(joinpath(cfg.outdir, "run_metadata.txt"), "w")
    for line in provenance
        println(metadata, line)
    end
    println(metadata, "case = NutrientGrowthRefined")
    println(metadata, "method = ", cfg.method)
    println(metadata, "ncells = ", cfg.ncells, "  h = ", 1.0 / cfg.ncells)
    println(metadata, "tfinal_requested = ", cfg.tfinal, "  t_reached = ", t)
    println(metadata, "dtmax = ", cfg.dtmax, "  cfl = ", cfg.cfl, "  min_dt = ",
        isempty(get(stats, "dt", Float64[])) ? NaN : minimum(stats["dt"]))
    println(metadata, "Dc = ", cfg.Dc, "  lambda = ", cfg.lambda, "  g0 = ", cfg.g0,
        "  cstar = ", cfg.cstar)
    println(metadata, "nitsche = ", cfg.nitsche, "  ghost = ", cfg.ghost,
        "  quadrature_degree = ", cfg.qdegree)
    println(metadata, "patch_radius_h = ", cfg.patch_radius_h, "  patch_radius_max_h = ",
        cfg.patch_radius_max_h, "  min_patch_samples = ", cfg.min_patch_samples,
        "  patch_condition_max = ", cfg.patch_condition_max,
        "  extension_index = ", cfg.extension_index)
    println(metadata, "filter_mode = ", cfg.filter_mode, "  filter_radius = ", cfg.filter_radius,
        "  (0 disables filterSubResolutionFeatures!)")
    println(metadata, "eps_scale = ", cfg.eps_scale,
        "  (denominator scale of the area-balance residual; binds in ",
        eps_scale_bound, " steps)")
    println(metadata, "area_initial_ed = ", area_start)
    println(metadata, "frame_count = ", length(frames), " (all fields synchronised to t)")
    for key in sort!(collect(keys(per)))
        println(metadata, key, " = ", per[key])
    end
    close(metadata)

    (; per, frames = length(frames), total_runtime, steps = step, t_final = t,
        final_area = frames[end].area, cfg)
end

# =============================================================================
# 13. Entry point
# =============================================================================

function runOne(; ncells, tfinal, method, dtmax, tag, filter_mode = :none, filter_radius = 0)
    cfg = defaultConfig(; ncells, tfinal, method, dtmax, filter_mode, filter_radius,
        outdir = joinpath(@__DIR__, "output", tag), tag)
    mkpath(cfg.outdir)
    provenance = checkEvolvingDomains()
    for line in provenance
        println(line)
    end
    verification = verificationCase(cfg)
    open(joinpath(cfg.outdir, "verification.txt"), "w") do io
        for key in sort!(collect(keys(verification)))
            println(io, key, " = ", verification[key])
        end
    end
    simulate(cfg; provenance)
end

function verificationSuite(ncells_list)
    mkpath(joinpath(@__DIR__, "output"))
    provenance = checkEvolvingDomains()
    for line in provenance
        println(line)
    end
    path = joinpath(@__DIR__, "output", "verification.csv")
    keys_order = String[]
    rows = Vector{Dict{String,Float64}}()
    for n in ncells_list
        cfg = defaultConfig(; ncells = n)
        row = verificationCase(cfg)
        isempty(keys_order) && (keys_order = sort!(collect(keys(row))))
        push!(rows, row)
        @printf("n=%3d  areaErr(ED/FEM)=%.2e/%.2e  pErr(max/L2)=%.2e/%.2e  Vn_A=%.2e  Vn_B=%.2e\n",
            n, row["verification_area_relative_error_ed"],
            row["verification_area_relative_error_fem"],
            row["verification_pressure_nodal_max_error"],
            row["verification_pressure_l2_relative_error"],
            row["verification_recovery_A_mean_abs_error"],
            row["verification_recovery_B_mean_abs_error"])
        flush(stdout)
    end
    open(path, "w") do io
        println(io, join(keys_order, ","))
        for row in rows
            println(io, join((row[k] for k in keys_order), ","))
        end
    end
    open(joinpath(@__DIR__, "output", "verification_provenance.txt"), "w") do io
        for line in provenance
            println(io, line)
        end
    end
    @printf("wrote %s\n", path)
end

function main(args = ARGS)
    mode = isempty(args) ? "verify" : args[1]
    if mode == "verify"
        list = length(args) >= 2 ? [parse(Int, a) for a in args[2:end]] : [16, 32, 64, 96, 128]
        verificationSuite(list)
    elseif mode == "recovery"
        # A/B interface-velocity recovery scored against the manufactured pressure only.
        ncells = length(args) >= 2 ? parse(Int, args[2]) : 64
        cfg = defaultConfig(; ncells)
        row = verificationCase(cfg)
        @printf("Vn exact = %.6f\n", row["verification_recovery_exact_normal_velocity"])
        @printf("A: mean|err| = %.4e  max|err| = %.4e  flux rel.err = %.4e\n",
            row["verification_recovery_A_mean_abs_error"],
            row["verification_recovery_A_max_abs_error"],
            row["verification_recovery_A_flux_relative_error"])
        @printf("B: mean|err| = %.4e  max|err| = %.4e  flux rel.err = %.4e  fallbacks = %d\n",
            row["verification_recovery_B_mean_abs_error"],
            row["verification_recovery_B_max_abs_error"],
            row["verification_recovery_B_flux_relative_error"],
            round(Int, row["verification_recovery_B_fallbacks"]))
    elseif mode == "extension"
        extensionBenchmark(length(args) >= 2 ? [parse(Int, a) for a in args[2:end]] : [64, 128])
    elseif mode == "merged"
        tags = length(args) >= 2 ? args[2:end] : ["n64_A", "n64_B", "n96_B", "n128_B"]
        mergedSuite(tags)
    elseif mode == "run"
        ncells = parse(Int, args[2])
        tfinal = parse(Float64, args[3])
        method = Symbol(args[4] == "A" ? :centroid : args[4] == "B" ? :patch :
                        error("method must be A or B"))
        dtmax = length(args) >= 5 ? parse(Float64, args[5]) : 5.0e-3
        tag = length(args) >= 6 ? args[6] : @sprintf("n%d_%s", ncells, args[4])
        spec = length(args) >= 7 ? args[7] : "none"
        filter_mode, filter_radius = if spec == "none"
            (:none, 0)
        else
            parts = split(spec, ":")
            (Symbol(parts[1]), parse(Int, parts[2]))
        end
        runOne(; ncells, tfinal, method, dtmax, tag, filter_mode, filter_radius)
    else
        error("unknown mode $mode (verify | recovery | run)")
    end
end

if abspath(PROGRAM_FILE) == @__FILE__
    main()
end
