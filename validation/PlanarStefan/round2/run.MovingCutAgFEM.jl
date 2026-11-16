#!/usr/bin/env julia
# MovingCutAgFEM -- moving-domain elliptic solve where EvolvingDomains supplies the cut
# geometry and StaticAssembler.jl assembles the discrete operator directly.
#
#   D = [0,1]^2
#   Omega(0) = {(x,y) : (x-0.38)^2 + (y-0.46)^2 < R^2},  R = 0.23
#   c(t) = (0.38 + 0.14 t, 0.46 + 0.08 t)                 (rigid translation, v = (0.14,0.08))
#   phi_t + v.grad phi = 0                                (numerical WENO5 transport)
#   -Lap u = 4 in Omega(t),  u = 0 on dOmega(t)
#   u_exact = R^2 - [(x-cx)^2 + (y-cy)^2]
#
# Architectural claim under test: EvolvingDomains' public `Cut` snapshot is a complete
# descriptor of the moving cut geometry, so StaticAssembler can assemble the projected
# (master-only) AgFEM system B = E'A E from it -- no GridapEmbedded re-cut.
#
# Run (from this directory):
#   julia --project=. run.MovingCutAgFEM.jl                     # all stages, n = 33
#   julia --project=. run.MovingCutAgFEM.jl --n 65 --stage D
#
# Stages (see REPORT.md):
#   A  version/descriptor audit only
#   B  Cut -> CutData adapter audit
#   C  one stationary manufactured solve
#   D  full moving-domain sequence
#   E  refinement, reuse, allocation and conditioning study
#   F  Ferrite cut-element quadrature probe
#   all (default) = A B C D E F

using EvolvingDomains
using StaticAssembler
import Ferrite
using Ferrite: CellValues, PointValues, QuadratureRule, RefQuadrilateral, Quadrilateral,
    Lagrange, Vec, Node, Grid, PointEvalHandler, getnquadpoints, getdetJdV, reinit!,
    shape_gradient, shape_value
using StaticArrays: SVector
using SparseArrays: nnz, sparse, SparseMatrixCSC
using LinearAlgebra
using Printf
using Statistics

const ED = EvolvingDomains
const SA = StaticAssembler

# =============================================================================
# 1. Variational form: -Lap u = 4 with symmetric Nitsche Dirichlet data
#
#   a_h(u,v) = ∫_Ω ∇u·∇v - ∫_Γ (∂n u) v - ∫_Γ (∂n v) u + (γ/h) ∫_Γ u v
#   l_h(v)   = ∫_Ω 4 v                  (Nitsche data g = 0, so the boundary load
#                                        terms +∫(∂n v) g and +∫(γ/h) g v vanish)
#
# StaticAssembler ships no PDE: differential operators are user types implementing its two
# documented contraction hooks (`_op_pair`, `_rhs_contract!`) plus the optional boundary
# specialisations `_boundary_op_pair` / `_boundary_rhs_contract!`. These four methods are
# the entire mathematics-to-assembly interface.
# =============================================================================

struct ValueOp <: SA.AbstractDifferentialOperator end
struct GradOp <: SA.AbstractDifferentialOperator end

@inline SA._op_pair(::ValueOp, ::ValueOp, row::Int, col::Int,
                    sh_t, gx_t, gy_t, sh_c, gx_c, gy_c) = sh_t[row] * sh_c[col]

@inline SA._op_pair(::GradOp, ::GradOp, row::Int, col::Int,
                    sh_t, gx_t, gy_t, sh_c, gx_c, gy_c) =
    gx_t[row] * gx_c[col] + gy_t[row] * gy_c[col]

# (∂n v) u  and  v (∂n u): the normal arrives through the boundary point frame.
@inline function SA._boundary_op_pair(::GradOp, ::ValueOp, row::Int, col::Int,
                                      sh_t, gx_t, gy_t, sh_c, gx_c, gy_c, point)
    n1, n2 = SA.boundary_normal(point)
    return (gx_t[row] * n1 + gy_t[row] * n2) * sh_c[col]
end

@inline function SA._boundary_op_pair(::ValueOp, ::GradOp, row::Int, col::Int,
                                      sh_t, gx_t, gy_t, sh_c, gx_c, gy_c, point)
    n1, n2 = SA.boundary_normal(point)
    return sh_t[row] * (gx_c[col] * n1 + gy_c[col] * n2)
end

@inline SA._boundary_op_pair(::ValueOp, ::ValueOp, row::Int, col::Int,
                             sh_t, gx_t, gy_t, sh_c, gx_c, gy_c, point) =
    sh_t[row] * sh_c[col]

@inline function SA._rhs_contract!(fe, ::ValueOp, slot, shape, gx, gy, c, w, offset::Int)
    value = slot.scale * slot.source(c) * w
    @inbounds for row in eachindex(shape)
        fe[offset + row] += value * shape[row]
    end
    return nothing
end

"""
    ellipticForm(penalty) -> ElementVariationalForm

`-Δu = 4` on the cut phase with homogeneous Dirichlet data imposed by symmetric Nitsche,
penalty `penalty`. The volume bilinear slot carries a `ConstantCoefficient`, so
StaticAssembler precomputes the interior and cut element stiffness blocks once per geometry
refresh; only the RHS source and the boundary slots are contracted pointwise.
"""
function ellipticForm(penalty::Float64)
    space = SA.FiniteElement(SA.Q1ReferenceQuadrilateral(), 1)
    return SA.ElementVariationalForm(
        space,
        (SA.BilinearSlot(GradOp(), GradOp(), SA.ConstantCoefficient(1.0)),),
        (SA.RhsSlot(ValueOp(), _ -> 4.0, 1.0),),
        0;
        boundary_bilinear = (
            SA.BilinearSlot(GradOp(), ValueOp(), SA.ConstantCoefficient(-1.0)),   # -∫ ∂n u v
            SA.BilinearSlot(ValueOp(), GradOp(), SA.ConstantCoefficient(-1.0)),   # -∫ ∂n v u
            SA.BilinearSlot(ValueOp(), ValueOp(), SA.ConstantCoefficient(penalty)),  # +γ/h ∫ u v
        ),
    )
end

# =============================================================================
# 2. Parameters and exact solution
# =============================================================================

function defaultParams(; n = 33, nout = 10, out = "output", gamma0 = 40.0, stage = "all")
    nsteps = 200                # dt = 2.5e-3 over 0 <= t <= 0.5 (CFL ~ 0.13 at n = 65)
    return (R = 0.23, c0 = (0.38, 0.46), v = (0.14, 0.08), tf = 0.5,
        nsteps = nsteps, dt = 0.5 / nsteps,
        n = n, nout = nout, out = out, gamma0 = gamma0, stage = stage)
end

center(p, t) = (p.c0[1] + p.v[1] * t, p.c0[2] + p.v[2] * t)

uExact(p, x, y, t) = begin
    cx, cy = center(p, t)
    p.R^2 - ((x - cx)^2 + (y - cy)^2)
end

# =============================================================================
# 3. EvolvingDomains geometry: Grid, level-set Field, numerical transport
# =============================================================================

function configureCutCellsLibrary()
    override = get(ENV, "EVOLVINGDOMAINS_CUTCELLS_LIBRARY", "")
    isempty(override) || return override
    ext = Sys.iswindows() ? "dll" : Sys.isapple() ? "dylib" : "so"
    local_lib = joinpath(@__DIR__, "build", "libevolvingdomains_cut_cells." * ext)
    isfile(local_lib) || error("no CutCells library: set EVOLVINGDOMAINS_CUTCELLS_LIBRARY")
    ENV["EVOLVINGDOMAINS_CUTCELLS_LIBRARY"] = local_lib
    return local_lib
end

"""
    buildEvolvingGeometry(p; n) -> (grid, discretization, field)

Uniform `n x n` node lattice on `D` plus the level-set Field holding the *exact* signed
distance at t = 0 only. After this call the field is advanced numerically; the exact circle
is never re-imposed.
"""
function buildEvolvingGeometry(p; n = p.n)
    x = collect(range(0.0, 1.0; length = n))
    d = ED.Discretization(x, x)
    cx, cy = center(p, 0.0)
    # EvolvingDomains convention: the negative phase (CUT_IN) is phi < 0, so phi is the
    # standard signed distance of the disk, negative inside Omega.
    phi = [((xi - cx)^2 + (yj - cy)^2) - p.R^2 for xi in x, yj in x]
    field = ED.Field(d, Dict(:init => phi))
    return d.levels[:init], d, field
end

transportVelocity(p, grid) = fill(SVector(p.v[1], p.v[2]), length(grid.axes[1]) * length(grid.axes[2]))

# =============================================================================
# 4. Adapter: EvolvingDomains `Cut` -> StaticAssembler `Mesh` / `CutData`
#
# Pure data conversion except for one structural difference: EvolvingDomains emits volume
# subcells for *every* background cell (an entirely-outside cell owns its own outside
# volume), while StaticAssembler requires subcells to be attached to CUT cells only, grouped
# by ascending background-cell index. The adapter keeps the subcells of CUT cells and
# renumbers their CSR pointers. Coordinates are converted, not recomputed: the point arrays
# stay read-only downstream, so no geometry is re-derived anywhere in the primary path.
# =============================================================================

# EvolvingDomains already stores the exact descriptor type; reuse the array when it matches
# and copy only when it does not. StaticAssembler reads these arrays, never mutates them.
points2(v::Vector{NTuple{2,Float64}}) = v
points2(v) = NTuple{2,Float64}[(Float64(q[1]), Float64(q[2])) for q in v]

const ED_IN, ED_CUT, ED_OUT = ED.CUT_IN, ED.CUT_CUT, ED.CUT_OUT

# Fail loudly if the two state encodings ever drift (they are identical today).
ED_IN == SA.IN && ED_CUT == SA.CUT && ED_OUT == SA.OUT ||
    error("EvolvingDomains and StaticAssembler cut-state encodings disagree")

function buildStaticAssemblerMesh(grid)
    axes = grid.axes
    length(axes) == 2 || throw(ArgumentError("this experiment is two-dimensional"))
    x, y = axes[1], axes[2]
    hx, hy = x[2] - x[1], y[2] - y[1]
    # StaticAssembler derives gradients from a single (hx,hy); the lattice must be uniform.
    for axis in axes
        step = axis[2] - axis[1]
        maximum(abs, diff(axis) .- step) <= 1e-12 * step ||
            throw(ArgumentError("StaticAssembler requires a uniform background lattice"))
    end
    return SA.Mesh{Float64}(grid.mesh_id, length(x) - 1, length(y) - 1,
        Float64(x[1]), Float64(y[1]), Float64(hx), Float64(hy))
end

"""
    convertCutToStaticAssembler(grid, cut) -> CutData{Float64}

Convert an EvolvingDomains `Cut` snapshot into the backend-agnostic descriptor
StaticAssembler consumes. Only the CUT cells' subcells survive, read through `Cut`'s own
background-cell -> subcell index (`cutCellIds`, `subcells`) rather than by re-filtering
`cell_to_bgcell`. The cells are ascending and the producer groups each cell's subcells, so
the emitted `parents` are ascending as StaticAssembler requires.
"""
function convertCutToStaticAssembler(grid, cut::ED.Cut)
    cut_cells = ED.cutCellIds(cut)
    nsub = sum(cell -> length(ED.subcells(cut, cell)), cut_cells)
    parents = Vector{Int32}(undef, nsub)
    states = Vector{Int8}(undef, nsub)
    point_ids = Int32[]
    sizehint!(point_ids, 3nsub)
    pointers = Vector{Int32}(undef, nsub + 1)
    k = 0
    for cell in cut_cells, subcell in ED.subcells(cut, cell)
        k += 1
        parents[k] = cut.cell_to_bgcell[subcell]
        states[k] = cut.subcell_states[subcell]
        pointers[k] = Int32(length(point_ids) + 1)
        for q in Int(cut.point_pointers[subcell]):(Int(cut.point_pointers[subcell + 1]) - 1)
            push!(point_ids, cut.point_ids[q])
        end
    end
    pointers[end] = Int32(length(point_ids) + 1)

    b = cut.boundary
    boundary_data = SA.BoundaryData{Float64}(
        Vector{Int32}(b.facet_to_bgcell),
        Vector{Int32}(b.point_ids),
        Vector{Int32}(b.point_pointers),
        points2(b.normals), points2(b.physical_points), points2(b.reference_points),
    )
    return SA.CutData{Float64}(
        grid.mesh_id,
        Vector{Int8}(cut.bg_states),
        states, parents, point_ids, pointers,
        points2(cut.physical_points),
        points2(cut.reference_points),
        boundary_data,
    )
end

# =============================================================================
# 4b. Index audit: `Cut`'s own cut-cell index vs the filtered scan it replaces
# =============================================================================

"""
    checkCutSubcellIndex(cut) -> NamedTuple

Verify that traversing `cut.cutCellIds` x `subcells(cut, cell)` reproduces, in the same
order, the subcell list obtained by scanning `cell_to_bgcell` and testing
`bg_states`. Returns the counts alongside the boolean so the run can report coverage.
"""
function checkCutSubcellIndex(cut)
    scan = Int32[subcell for subcell in eachindex(cut.subcell_states)
                 if cut.bg_states[Int(cut.cell_to_bgcell[subcell])] == ED_CUT]
    cut_cells = ED.cutCellIds(cut)
    indexed = Int32[subcell for cell in cut_cells for subcell in ED.subcells(cut, cell)]
    return (identical = indexed == scan, n_cut_cells = length(cut_cells),
        n_cut_subcells = length(indexed), n_subcells = length(cut.subcell_states))
end

# =============================================================================
# 5. StaticAssembler driver: refresh, assembly, solve
# =============================================================================

"""
    refreshDiscreteOperator!(system, grid, cut) -> (cutdata, roots, structural)

Rebuilds the descriptors from the current `Cut` snapshot and refreshes the existing
`StaticSystem` in place. `structural` reports whether StaticAssembler replaced its
projection tables (E and B's sparse pattern) instead of reusing them.
"""
function refreshDiscreteOperator!(system, grid, cut)
    convert = timeit(() -> convertCutToStaticAssembler(grid, cut))
    aggregate = timeit(() -> ED.aggregateCutCells(cut, grid))
    E_before = objectid(SA.extension(system))
    refresh = timeit(() -> SA.refresh_geometry!(system, convert.value, aggregate.value))
    return (cutdata = convert.value, roots = aggregate.value,
        structural = objectid(SA.extension(system)) != E_before,
        convert = convert, aggregate = aggregate, refresh = refresh)
end

struct EllipticSolution
    u_m::Vector{Float64}
    u_bg::Vector{Float64}
    B::SparseMatrixCSC{Float64,Int}
    b::Vector{Float64}
    stage::NamedTuple
    cholesky_ok::Bool
end

timeit(f) = (s = @timed f(); (value = s.value, time = s.time, bytes = s.bytes))

"""
    solveElliptic(system) -> EllipticSolution

Assemble, factor, solve and extend: `B u_m = b_m` by sparse Cholesky (CHOLMOD) on the
projected operator, then `u_bg = E u_m`. Each stage is timed and its allocation measured.
"""
function solveElliptic(system; threaded = false)
    a = timeit(() -> SA.assemble!(system, (); threaded = threaded))
    B, b = a.value
    rhs = Vector(b)                       # the system-owned buffer is overwritten next call
    f = timeit(() -> cholesky(Symmetric(B); check = false))
    fact = f.value
    ok = issuccess(fact)
    s = timeit(() -> ok ? Vector(fact \ rhs) : Vector(B \ rhs))
    u_m = s.value
    r = timeit(() -> SA.extension(system) * u_m)
    return EllipticSolution(u_m, r.value, B, rhs,
        (assemble = a, factor = f, solve = s, recover = r), ok)
end

# --- finite-element evaluation on the background lattice ---------------------

@inline function cell_offsets(grid, cell)
    nx = length(grid.axes[1]) - 1
    return (cell - 1) % nx, (cell - 1) ÷ nx
end

"""Background lattice node index of local Q1 node `k` (1..4) of `cell`."""
@inline function q1Node(grid, cell, k)
    i, j = cell_offsets(grid, cell)
    n = length(grid.axes[1])
    base = i + j * n
    return k == 1 ? base + 1 : k == 2 ? base + 2 : k == 3 ? base + n + 1 : base + n + 2
end

"""Bilinear Q1 evaluation of a whole background node vector at `(r,s)` in `cell`."""
@inline function q1Evaluate(u, grid, cell, r, s)
    return (1 - r) * (1 - s) * u[q1Node(grid, cell, 1)] +
           r * (1 - s) * u[q1Node(grid, cell, 2)] +
           (1 - r) * s * u[q1Node(grid, cell, 3)] +
           r * s * u[q1Node(grid, cell, 4)]
end

const TRI_POINTS = ((1/6, 1/6, 2/3), (1/6, 2/3, 1/6), (2/3, 1/6, 1/6))
const TRI_WEIGHTS = (1/3, 1/3, 1/3)
const QUAD_POINTS = (0.21132486540518713, 0.78867513459481287)   # 2x2 Gauss on [0,1]^2

"""
    integrateDomain(fval, grid, cutdata) -> Float64

Integrate `fval(x, y, r, s, cell)` over the *reconstructed* numerical domain: 2x2 Gauss on
fully-interior background cells, 3-point Gauss on the inside subtriangles of cut cells.
"""
function integrateDomain(fval, grid, cutdata)
    x0, y0 = grid.axes[1][1], grid.axes[2][1]
    hx, hy = grid.axes[1][2] - grid.axes[1][1], grid.axes[2][2] - grid.axes[2][1]
    ncells = (length(grid.axes[1]) - 1) * (length(grid.axes[2]) - 1)
    total = 0.0
    for cell in 1:ncells
        cutdata.bg_states[cell] == SA.IN || continue
        i, j = cell_offsets(grid, cell)
        for r in QUAD_POINTS, s in QUAD_POINTS
            total += hx * hy * 0.25 *
                     fval(x0 + (i + r) * hx, y0 + (j + s) * hy, r, s, cell)
        end
    end
    for subcell in eachindex(cutdata.subcell_states)
        cutdata.subcell_states[subcell] == SA.IN || continue
        parent = Int(cutdata.cell_to_bgcell[subcell])
        k = Int(cutdata.point_pointers[subcell])
        i1, i2, i3 = Int(cutdata.point_ids[k]), Int(cutdata.point_ids[k + 1]),
                      Int(cutdata.point_ids[k + 2])
        r1, r2, r3 = cutdata.reference_points[i1], cutdata.reference_points[i2],
                     cutdata.reference_points[i3]
        x1, x2, x3 = cutdata.physical_points[i1], cutdata.physical_points[i2],
                     cutdata.physical_points[i3]
        area = abs((x2[1] - x1[1]) * (x3[2] - x1[2]) -
                   (x3[1] - x1[1]) * (x2[2] - x1[2])) / 2
        area > 0 || continue
        for q in 1:3
            lam1, lam2, lam3 = TRI_POINTS[q]
            r = lam1 * r1[1] + lam2 * r2[1] + lam3 * r3[1]
            s = lam1 * r1[2] + lam2 * r2[2] + lam3 * r3[2]
            x = lam1 * x1[1] + lam2 * x2[1] + lam3 * x3[1]
            y = lam1 * x1[2] + lam2 * x2[2] + lam3 * x3[2]
            total += area * TRI_WEIGHTS[q] * fval(x, y, r, s, parent)
        end
    end
    return total
end

"""Area and centroid of the reconstructed domain (interior cells + inside cut pieces)."""
function domainAreaCenter(grid, cutdata)
    area = 0.0
    mx = 0.0
    my = 0.0
    hx = grid.axes[1][2] - grid.axes[1][1]
    hy = grid.axes[2][2] - grid.axes[2][1]
    ncells = (length(grid.axes[1]) - 1) * (length(grid.axes[2]) - 1)
    for cell in 1:ncells
        cutdata.bg_states[cell] == SA.IN || continue
        i, j = cell_offsets(grid, cell)
        area += hx * hy
        mx += hx * hy * (grid.axes[1][1] + (i + 0.5) * hx)
        my += hx * hy * (grid.axes[2][1] + (j + 0.5) * hy)
    end
    for subcell in eachindex(cutdata.subcell_states)
        cutdata.subcell_states[subcell] == SA.IN || continue
        k = Int(cutdata.point_pointers[subcell])
        ids = (Int(cutdata.point_ids[k]), Int(cutdata.point_ids[k + 1]),
               Int(cutdata.point_ids[k + 2]))
        a, b, c = cutdata.physical_points[ids[1]], cutdata.physical_points[ids[2]],
                  cutdata.physical_points[ids[3]]
        s = abs((b[1] - a[1]) * (c[2] - a[2]) - (c[1] - a[1]) * (b[2] - a[2])) / 2
        s > 0 || continue
        area += s
        mx += s * (a[1] + b[1] + c[1]) / 3
        my += s * (a[2] + b[2] + c[2]) / 3
    end
    return area, (mx / area, my / area)
end

"""Smallest cut-cell area fraction among cut cells that carry quadrature."""
function minCutAreaFraction(grid, cutdata)
    hx = grid.axes[1][2] - grid.axes[1][1]
    hy = grid.axes[2][2] - grid.axes[2][1]
    areas = Dict{Int,Float64}()
    for subcell in eachindex(cutdata.subcell_states)
        cutdata.subcell_states[subcell] == SA.IN || continue
        k = Int(cutdata.point_pointers[subcell])
        ids = (Int(cutdata.point_ids[k]), Int(cutdata.point_ids[k + 1]),
               Int(cutdata.point_ids[k + 2]))
        a, b, c = cutdata.physical_points[ids[1]], cutdata.physical_points[ids[2]],
                  cutdata.physical_points[ids[3]]
        s = abs((b[1] - a[1]) * (c[2] - a[2]) - (c[1] - a[1]) * (b[2] - a[2])) / 2
        s > 0 || continue
        cell = Int(cutdata.cell_to_bgcell[subcell])
        areas[cell] = get(areas, cell, 0.0) + s
    end
    isempty(areas) && return NaN, 0
    return minimum(values(areas)) / (hx * hy), length(areas)
end

# =============================================================================
# 6. Diagnostics
# =============================================================================

"""
    checkManufacturedSolution(p, grid, cutdata, solution, t; system, roots) -> NamedTuple

Relative L2 error, interior maximum error, Nitsche boundary error, algebraic residual,
symmetry of B, DOF counts, aggregate count, cut-cell geometry quality, and the separate
interpolation error of the exact solution (the discretisation-only share of the total).
"""
function checkManufacturedSolution(p, grid, cutdata, solution::EllipticSolution, t;
                                   system, roots)
    u = solution.u_bg
    err2(x, y, r, s, cell) = (q1Evaluate(u, grid, cell, r, s) - uExact(p, x, y, t))^2
    num = integrateDomain(err2, grid, cutdata)
    den = integrateDomain((x, y, r, s, cell) -> uExact(p, x, y, t)^2, grid, cutdata)

    # discretisation-only share: the exact solution interpolated on the same lattice
    nx1 = length(grid.axes[1])
    ny1 = length(grid.axes[2])
    u_interp = Vector{Float64}(undef, nx1 * ny1)
    for j in 1:ny1, i in 1:nx1
        u_interp[i + (j - 1) * nx1] = uExact(p, grid.axes[1][i], grid.axes[2][j], t)
    end
    ierr2(x, y, r, s, cell) = (q1Evaluate(u_interp, grid, cell, r, s) - uExact(p, x, y, t))^2
    inum = integrateDomain(ierr2, grid, cutdata)

    # maximum error on safely interior lattice nodes (>= 1.5 h away from the interface)
    cx, cy = center(p, t)
    x0, y0 = grid.axes[1][1], grid.axes[2][1]
    hx = grid.axes[1][2] - grid.axes[1][1]
    hy = grid.axes[2][2] - grid.axes[2][1]
    maxinterior = 0.0
    for j in 1:ny1, i in 1:nx1
        x, y = grid.axes[1][i], grid.axes[2][j]
        sqrt((x - cx)^2 + (y - cy)^2) < p.R - 1.5max(hx, hy) || continue
        err = abs(u[i + (j - 1) * nx1] - uExact(p, x, y, t))
        err > maxinterior && (maxinterior = err)
    end

    # boundary-condition error: max |u_h| at embedded-boundary facet midpoints
    bmax = 0.0
    boundary = cutdata.boundary
    for facet in eachindex(boundary.facet_to_bgcell)
        k = Int(boundary.point_pointers[facet])
        a = boundary.physical_points[boundary.point_ids[k]]
        b = boundary.physical_points[boundary.point_ids[k + 1]]
        x, y = (a[1] + b[1]) / 2, (a[2] + b[2]) / 2
        cell = Int(boundary.facet_to_bgcell[facet])
        i, j = cell_offsets(grid, cell)
        v = abs(q1Evaluate(u, grid, cell, (x - (x0 + i * hx)) / hx, (y - (y0 + j * hy)) / hy))
        v > bmax && (bmax = v)
    end

    B, rhs, u_m = solution.B, solution.b, solution.u_m
    resid = norm(B * u_m - rhs) / max(norm(rhs), eps())
    asym = norm(B - transpose(B)) / max(norm(B), eps())
    E = SA.extension(system)
    nact = count(!iszero, vec(sum(abs, E; dims = 2)))
    area, ctr = domainAreaCenter(grid, cutdata)
    minfrac, ncut = minCutAreaFraction(grid, cutdata)
    cxe, cye = center(p, t)

    return (t = t, l2_rel = sqrt(num / den), l2_interp_rel = sqrt(inum / den),
        max_interior = maxinterior, bc_error = bmax, residual = resid, asymmetry = asym,
        n_master = size(E, 2), n_active_bg = nact,
        n_aggregates = length(unique(filter(>(0), roots))),
        min_cut_fraction = minfrac, n_cut_cells = ncut,
        area = area, area_exact = pi * p.R^2,
        area_rel_error = abs(area - pi * p.R^2) / (pi * p.R^2),
        center_x = ctr[1], center_y = ctr[2], center_error = hypot(ctr[1] - cxe, ctr[2] - cye),
        nnz_B = nnz(B), nnz_E = nnz(E), cholesky_ok = solution.cholesky_ok)
end

# =============================================================================
# 7. Ferrite.jl interoperability probe
# =============================================================================

function probeFerriteQuadrature(p; path = joinpath(@__DIR__, "output"))
    grid, _, field = buildEvolvingGeometry(p)
    cut = ED.cut(field, :init)
    hx = grid.axes[1][2] - grid.axes[1][1]
    hy = grid.axes[2][2] - grid.axes[2][1]
    x0, y0 = grid.axes[1][1], grid.axes[2][1]

    # --- pick a genuinely half-cut cell: total inside area fraction nearest 0.5 ---
    inside = Dict{Int,Float64}()
    for cell in ED.cutCellIds(cut), subcell in ED.subcells(cut, cell)
        cut.subcell_states[subcell] == ED_IN || continue
        inside[cell] = get(inside, cell, 0.0) + _subcellArea(cut, subcell)
    end
    isempty(inside) && error("no cut cell found for the Ferrite probe")
    best_cell = argmin(Dict(k => abs(v / (hx * hy) - 0.5) for (k, v) in inside))
    best_frac = inside[best_cell] / (hx * hy)

    # quadrature from EvolvingDomains geometry: inside subtriangles of that cell
    ξs = Vec{2,Float64}[]      # Ferrite reference frame, [-1,1]^2
    ws = Float64[]             # physical weights
    xs = Vec{2,Float64}[]
    for subcell in ED.subcells(cut, best_cell)
        cut.subcell_states[subcell] == ED_IN || continue
        k = Int(cut.point_pointers[subcell])
        ids = (Int(cut.point_ids[k]), Int(cut.point_ids[k + 1]), Int(cut.point_ids[k + 2]))
        r1, r2, r3 = cut.reference_points[ids[1]], cut.reference_points[ids[2]],
                     cut.reference_points[ids[3]]
        x1, x2, x3 = cut.physical_points[ids[1]], cut.physical_points[ids[2]],
                     cut.physical_points[ids[3]]
        area = abs((x2[1] - x1[1]) * (x3[2] - x1[2]) - (x3[1] - x1[1]) * (x2[2] - x1[2])) / 2
        area > 1e-14 || continue
        for q in 1:3
            l1, l2, l3 = TRI_POINTS[q]
            r = l1 * r1[1] + l2 * r2[1] + l3 * r3[1]
            s = l1 * r1[2] + l2 * r2[2] + l3 * r3[2]
            px = l1 * x1[1] + l2 * x2[1] + l3 * x3[1]
            py = l1 * x1[2] + l2 * x2[2] + l3 * x3[2]
            push!(ξs, Vec{2,Float64}((2r - 1, 2s - 1)))   # EvolvingDomains [0,1]^2 -> Ferrite [-1,1]^2
            push!(ws, area * TRI_WEIGHTS[q])
            push!(xs, Vec{2,Float64}((px, py)))
        end
    end

    # --- Ferrite: background cell in Ferrite's cyclic CCW node order --------
    i, j = cell_offsets(grid, best_cell)
    corner(r, s) = Vec{2,Float64}((x0 + (i + r) * hx, y0 + (j + s) * hy))
    # StaticAssembler/EVOLVINGDOMAINS Q1 order is (0,0),(1,0),(0,1),(1,1); Ferrite needs
    # (0,0),(1,0),(1,1),(0,1) -- the 3<->4 permutation.
    X = [corner(0.0, 0.0), corner(1.0, 0.0), corner(1.0, 1.0), corner(0.0, 1.0)]
    nb = 4
    ip = Lagrange{RefQuadrilateral,1}()
    qr = QuadratureRule{RefQuadrilateral}(ws, ξs)
    cv = CellValues(qr, ip; update_detJdV = false)   # physical weights: detJ must not be reapplied
    reinit!(cv, X)
    Kfer = zeros(nb, nb)
    for q in 1:getnquadpoints(cv), a in 1:nb, b in 1:nb
        Kfer[a, b] += ws[q] * (shape_gradient(cv, q, a) ⋅ shape_gradient(cv, q, b))
    end

    # route 2: PointValues, one reinit! per quadrature point
    pv = PointValues(Float64, ip)
    Kpv = zeros(nb, nb)
    for (q, ξq) in enumerate(ξs)
        reinit!(pv, X, ξq)
        for a in 1:nb, b in 1:nb
            Kpv[a, b] += ws[q] * (shape_gradient(pv, 1, a) ⋅ shape_gradient(pv, 1, b))
        end
    end

    # --- independent reference: hand-coded bilinear Q1 in StaticAssembler's order
    # explicit bilinear gradients, StaticAssembler Q1 order
    handgx(r, s, k) = k == 1 ? -(1 - s) / hx : k == 2 ? (1 - s) / hx :
                      k == 3 ? -s / hx : s / hx
    handgy(r, s, k) = k == 1 ? -(1 - r) / hy : k == 2 ? -r / hy :
                      k == 3 ? (1 - r) / hy : r / hy
    Khand = zeros(nb, nb)
    for (q, ξq) in enumerate(ξs)
        r, s = (ξq[1] + 1) / 2, (ξq[2] + 1) / 2
        for a in 1:nb, b in 1:nb
            Khand[a, b] += ws[q] * (handgx(r, s, a) * handgx(r, s, b) +
                                    handgy(r, s, a) * handgy(r, s, b))
        end
    end
    # permute Ferrite's (0,0),(1,0),(1,1),(0,1) matrix into (0,0),(1,0),(0,1),(1,1)
    perm = [1, 2, 4, 3]
    KferSA = Kfer[perm, perm]
    KpvSA = Kpv[perm, perm]

    # --- physical -> reference round trip through the public inverse map ----
    fgrid = Ferrite.Grid([Quadrilateral((1, 2, 3, 4))],
        [Node(X[1]), Node(X[2]), Node(X[3]), Node(X[4])])
    handler = PointEvalHandler(fgrid, xs)
    ξcheck = maximum(abs(handler.local_coords[q][d] - ξs[q][d])
                     for q in eachindex(xs), d in 1:2)
    allmissing = any(isnothing, handler.cells)

    # --- the same cell's outside pieces must complete the full-cell matrix --
    ξout = Vec{2,Float64}[]
    wout = Float64[]
    for subcell in ED.subcells(cut, best_cell)
        cut.subcell_states[subcell] == ED_OUT || continue
        k = Int(cut.point_pointers[subcell])
        ids = (Int(cut.point_ids[k]), Int(cut.point_ids[k + 1]), Int(cut.point_ids[k + 2]))
        r1, r2, r3 = cut.reference_points[ids[1]], cut.reference_points[ids[2]],
                     cut.reference_points[ids[3]]
        x1, x2, x3 = cut.physical_points[ids[1]], cut.physical_points[ids[2]],
                     cut.physical_points[ids[3]]
        area = abs((x2[1] - x1[1]) * (x3[2] - x1[2]) - (x3[1] - x1[1]) * (x2[2] - x1[2])) / 2
        area > 1e-14 || continue
        for q in 1:3
            l1, l2, l3 = TRI_POINTS[q]
            r = l1 * r1[1] + l2 * r2[1] + l3 * r3[1]
            s = l1 * r1[2] + l2 * r2[2] + l3 * r3[2]
            push!(ξout, Vec{2,Float64}((2r - 1, 2s - 1)))
            push!(wout, area * TRI_WEIGHTS[q])
        end
    end
    Kfull = zeros(nb, nb)
    if !isempty(ξout)
        cvout = CellValues(QuadratureRule{RefQuadrilateral}(wout, ξout), ip; update_detJdV = false)
        reinit!(cvout, X)
        for q in 1:getnquadpoints(cvout), a in 1:nb, b in 1:nb
            Kfull[a, b] += wout[q] * (shape_gradient(cvout, q, a) ⋅ shape_gradient(cvout, q, b))
        end
        Kfull += Kfer
    end

    # full-cell reference matrix from Ferrite's own standard rule (volume = 2x2 Gauss)
    qr2 = QuadratureRule{RefQuadrilateral}(2)   # default :legendre, 2x2 Gauss
    cv2 = CellValues(qr2, ip)
    reinit!(cv2, X)
    Kstd = zeros(nb, nb)
    for q in 1:getnquadpoints(cv2), a in 1:nb, b in 1:nb
        Kstd[a, b] += getdetJdV(cv2, q) *
                      (shape_gradient(cv2, q, a) ⋅ shape_gradient(cv2, q, b))
    end

    checks = (n = p.n, cut_cell = best_cell, area_fraction = best_frac,
        nquadpoints = length(ws),
        quad_area = sum(ws), cell_area = hx * hy,
        symmetry = maximum(abs, KferSA - KferSA'),
        min_eig = minimum(eigvals(Symmetric(KferSA))),
        max_eig = maximum(eigvals(Symmetric(KferSA))),
        constant_rowsum = maximum(abs, KferSA * ones(nb)),
        cellvalues_vs_pointvalues = maximum(abs, KferSA - KpvSA),
        ferrite_vs_handcoded = maximum(abs, KferSA - Khand),
        inverse_map_error = ξcheck,
        all_points_found = !allmissing,
        inside_plus_outside_vs_fullcell = maximum(abs, Kfull - Kstd))
    writeRows(joinpath(path, "ferrite_probe_n$(p.n).csv"), [checks])
    open(joinpath(path, "ferrite_quadrature_n$(p.n).csv"), "w") do io
        println(io, "cell,x,y,xi,eta,w")
        for q in eachindex(ws)
            @printf(io, "%d,%.17g,%.17g,%.17g,%.17g,%.17g\n", best_cell, xs[q][1], xs[q][2],
                ξs[q][1], ξs[q][2], ws[q])
        end
    end
    @printf("stage F: Ferrite %s | cut cell %d (fraction %.3f), %d quadrature points\n",
        string(pkgversion(Ferrite)), best_cell, best_frac, length(ws))
    @printf("stage F: symmetry=%.2e min_eig=%.2e rowsum=%.2e | vs(PV)=%.2e vs(hand)=%.2e | inv_map=%.2e | in+out=full %.2e\n",
        checks.symmetry, checks.min_eig, checks.constant_rowsum,
        checks.cellvalues_vs_pointvalues, checks.ferrite_vs_handcoded,
        checks.inverse_map_error, checks.inside_plus_outside_vs_fullcell)
    return checks
end

function _subcellArea(cut, subcell)
    k = Int(cut.point_pointers[subcell])
    ids = (Int(cut.point_ids[k]), Int(cut.point_ids[k + 1]), Int(cut.point_ids[k + 2]))
    a, b, c = cut.physical_points[ids[1]], cut.physical_points[ids[2]],
              cut.physical_points[ids[3]]
    return abs((b[1] - a[1]) * (c[2] - a[2]) - (c[1] - a[1]) * (b[2] - a[2])) / 2
end

# =============================================================================
# 8. CSV output
# =============================================================================

function writeRows(file, rows; digits = 17)
    isempty(rows) && return nothing
    cols = propertynames(first(rows))
    open(file, "w") do io
        println(io, join(cols, ","))
        for row in rows
            println(io, join((string(typeof(getproperty(row, c)) <: AbstractFloat ?
                                      round(getproperty(row, c); digits = digits) :
                                      getproperty(row, c)) for c in cols), ","))
        end
    end
    return nothing
end

function saveSnapshot(outdir, grid, system, sol, tag)
    E = SA.extension(system)
    open(joinpath(outdir, "E_$(tag).csv"), "w") do io
        println(io, "row,col,value")
        for col in 1:size(E, 2), k in E.colptr[col]:(E.colptr[col + 1] - 1)
            @printf(io, "%d,%d,%.17g\n", E.rowval[k], col, E.nzval[k])
        end
    end
    open(joinpath(outdir, "B_$(tag).csv"), "w") do io
        println(io, "# n=$(size(sol.B,1)) nnz=$(nnz(sol.B))")
        println(io, "row,col,value")
        for col in 1:size(sol.B, 2), k in sol.B.colptr[col]:(sol.B.colptr[col + 1] - 1)
            @printf(io, "%d,%d,%.17g\n", sol.B.rowval[k], col, sol.B.nzval[k])
        end
    end
    open(joinpath(outdir, "rhs_$(tag).csv"), "w") do io
        println(io, "i,value")
        for i in eachindex(sol.b)
            @printf(io, "%d,%.17g\n", i, sol.b[i])
        end
    end
    open(joinpath(outdir, "u_bg_$(tag).csv"), "w") do io
        println(io, "x,y,u")
        nx1 = length(grid.axes[1])
        for j in 1:length(grid.axes[2]), i in 1:nx1
            @printf(io, "%.17g,%.17g,%.17g\n", grid.axes[1][i], grid.axes[2][j],
                sol.u_bg[i + (j - 1) * nx1])
        end
    end
    return nothing
end

# =============================================================================
# 9. Stage E: reuse, allocations, conditioning
# =============================================================================

"""`λmax` by power iteration and `λmin ≈ 1/λmax(B⁻¹)` by inverse iteration."""
function spectralEstimate(B, fact; iterations = 50)
    n = size(B, 1)
    x = fill(1.0, n) ./ sqrt(n)
    λmax = 0.0
    for _ in 1:iterations
        y = B * x
        ny = norm(y)
        ny == 0 && return (NaN, NaN)
        x = y ./ ny
        λmax = dot(x, B * x)
    end
    z = fill(1.0, n) ./ sqrt(n)
    for _ in 1:iterations
        w = fact \ z
        nw = norm(w)
        nw == 0 && return (NaN, NaN)
        z = w ./ nw
    end
    λinv = dot(z, fact \ z)
    return (1 / λinv, λmax)
end

"""
    conditioningStudy(p, mesh, cutdata, roots, form) -> NamedTuple

Spectral estimate of the projected operator on the coarsest grid, compared with the same
mesh and cut geometry assembled *without* aggregation (every active cell its own root), which
is the ill-conditioned cut-FEM baseline AgFEM is meant to repair.
"""
function conditioningStudy(p, mesh, cutdata, roots, form)
    system = SA.StaticSystem(mesh, cutdata, roots; form = form)
    B, _ = SA.assemble!(system, (); threaded = false)
    pred = cholesky(Symmetric(B); check = false)
    n_agf = size(B, 1)
    eig_agf = if n_agf <= 1500
        extrema(eigvals(Symmetric(Matrix(B))))       # dense: modest grids only
    elseif issuccess(pred)
        spectralEstimate(B, pred)
    else
        (NaN, NaN)
    end

    selfroots = Int32[cutdata.bg_states[c] == SA.OUT ? 0 : c
                      for c in eachindex(cutdata.bg_states)]
    control = SA.StaticSystem(mesh, cutdata, selfroots; form = form)
    Bc, _ = SA.assemble!(control, (); threaded = false)
    predc = cholesky(Symmetric(Bc); check = false)
    n_ctl = size(Bc, 1)
    eig_ctl = if n_ctl <= 1500
        extrema(eigvals(Symmetric(Matrix(Bc))))
    elseif issuccess(predc)
        spectralEstimate(Bc, predc)
    else
        (NaN, NaN)
    end
    maxdiag = maximum(diag(Bc))

    return (n_m_agfem = n_agf, min_eig_agfem = eig_agf[1], max_eig_agfem = eig_agf[2],
        cond_agfem = eig_agf[2] / eig_agf[1], cholesky_agfem = issuccess(pred),
        min_diag_agfem = minimum(diag(B)),
        n_m_selfroot = n_ctl, min_eig_selfroot = eig_ctl[1],
        max_eig_selfroot = eig_ctl[2], cond_selfroot = eig_ctl[2] / eig_ctl[1],
        cholesky_selfroot = issuccess(predc), min_diag_selfroot = minimum(diag(Bc)),
        tiny_diag_selfroot = count(<(1e-10 * maxdiag), diag(Bc)))
end

# =============================================================================
# 10. Simulation
# =============================================================================

function simulate(p; n = p.n, stage = p.stage)
    grid, _, field = buildEvolvingGeometry(p; n = n)
    mesh = buildStaticAssemblerMesh(grid)
    velocity = transportVelocity(p, grid)
    outdir = joinpath(@__DIR__, p.out)
    isdir(outdir) || mkpath(outdir)

    if get(ENV, "MOVINGCUTAGFEM_VERSIONS", "1") == "1"
        ENV["MOVINGCUTAGFEM_VERSIONS"] = "0"
        ext = Base.get_extension(SA, :StaticAssemblerGridapExt)
        open(joinpath(outdir, "versions.txt"), "w") do io
            println(io, "julia            = ", VERSION)
            println(io, "EvolvingDomains  = ", pkgversion(ED), "  path=", pkgdir(ED))
            println(io, "StaticAssembler  = ", pkgversion(SA), "  path=", pkgdir(SA))
            println(io, "Ferrite          = ", pkgversion(Ferrite))
            println(io, "Gridap extension = ", ext === nothing ? "NOT LOADED (no GridapEmbedded in the primary path)" : string(ext))
            git(dir) = strip(read(`git -C $dir rev-parse HEAD`, String))
            println(io, "EvolvingDomains  git = ", git(pkgdir(ED)))
            println(io, "StaticAssembler  git = ", git(pkgdir(SA)))
            println(io, "EvolvingDomains  git status:")
            for line in eachline(`git -C $(pkgdir(ED)) status --porcelain`)
                println(io, "  ", line)
            end
            println(io, "StaticAssembler  git status:")
            for line in eachline(`git -C $(pkgdir(SA)) status --porcelain`)
                println(io, "  ", line)
            end
            println(io, "CutCells library = ", get(ENV, "EVOLVINGDOMAINS_CUTCELLS_LIBRARY", "artifact"))
            println(io, "n                = ", n, "   nsteps = ", p.nsteps, "   dt = ", p.dt)
            println(io, "gamma0           = ", p.gamma0, "  (penalty = gamma0 / min(hx,hy))")
            println(io, "threads          = ", Threads.nthreads())
        end
    end

    cut0 = ED.cut(field, :init)
    cutdata0 = convertCutToStaticAssembler(grid, cut0)
    roots0 = Int32.(ED.aggregateCutCells(cut0, grid))
    index0 = checkCutSubcellIndex(cut0)
    index0.identical || error("Cut subcell index disagrees with the filtered scan")

    # --- Stage A/B: descriptor audit ----------------------------------------
    if stage in ("A", "B")
        @printf("stage A/B n=%d: nodes=%d cells=%d cut=%d in=%d out=%d | kept subcells=%d/%d facets=%d\n",
            n, length(grid.axes[1])^2, length(cut0.bg_states),
            count(==(ED_CUT), cut0.bg_states), count(==(ED_IN), cut0.bg_states),
            count(==(ED_OUT), cut0.bg_states), length(cutdata0.cell_to_bgcell),
            length(cut0.subcell_states), length(cutdata0.boundary.facet_to_bgcell))
        @printf("stage A/B n=%d: aggregate roots nonzero=%d\n",
            n, length(unique(filter(>(0), roots0))))
        @printf("stage A/B n=%d: subcell index == scan (%d cut cells, %d/%d subcells kept)\n",
            n, index0.n_cut_cells, index0.n_cut_subcells, index0.n_subcells)
        return nothing
    end

    h = min(grid.axes[1][2] - grid.axes[1][1], grid.axes[2][2] - grid.axes[2][1])
    form = ellipticForm(p.gamma0 / h)

    # --- Stage C: cold construction and one stationary solve -----------------
    cold_jit = timeit(() -> SA.StaticSystem(mesh, cutdata0, roots0; form = form))
    cold = timeit(() -> SA.StaticSystem(mesh, cutdata0, roots0; form = form))
    system = cold.value
    sol = solveElliptic(system)
    diag0 = checkManufacturedSolution(p, grid, cutdata0, sol, 0.0; system = system, roots = roots0)
    @printf("stage C n=%d: cold StaticSystem %.1f ms / %.2f MB (first call incl. JIT: %.1f ms / %.2f MB)\n",
        n, 1e3cold.time, cold.bytes / 2^20, 1e3cold_jit.time, cold_jit.bytes / 2^20)
    @printf("stage C n=%d: t=0 L2rel=%.4e (interp-only %.4e) interior=%.4e bc=%.3e asym=%.2e res=%.2e\n",
        n, diag0.l2_rel, diag0.l2_interp_rel, diag0.max_interior, diag0.bc_error,
        diag0.asymmetry, diag0.residual)
    @printf("stage C n=%d: n_m=%d active_bg=%d aggregates=%d min_cut_fraction=%.3e area_relerr=%.3e\n",
        n, diag0.n_master, diag0.n_active_bg, diag0.n_aggregates, diag0.min_cut_fraction,
        diag0.area_rel_error)
    stage == "C" && (saveSnapshot(outdir, grid, system, sol, "n$(n)_t000"); return [diag0])

    # --- reuse baseline: refresh with an *unchanged* geometry ----------------
    E_id = objectid(SA.extension(system))
    warm = timeit(() -> SA.refresh_geometry!(system, cutdata0, roots0))
    warm_structural = objectid(SA.extension(system)) != E_id
    warm_asm = timeit(() -> SA.assemble!(system, (); threaded = false))
    SA.assemble!(system, (); threaded = true)          # warm up: first threaded call compiles
    thr = timeit(() -> SA.assemble!(system, (); threaded = true))
    @printf("stage E n=%d: warm refresh (same geometry) %.2f ms / %.1f kB structural=%s | assemble %.2f ms / %.1f kB (threaded %.2f ms)\n",
        n, 1e3warm.time, warm.bytes / 1024, string(warm_structural),
        1e3warm_asm.time, warm_asm.bytes / 1024, 1e3thr.time)

    # --- Stage D/E: moving-domain sequence ---------------------------------
    rows = Any[diag0]
    steps = NamedTuple[]
    diag_every = max(1, round(Int, p.nsteps / max(p.nout, 1)))
    for step in 0:p.nsteps
        if step > 0
            ED.eulerianTransport!(field, velocity, p.dt; level_id = :init)
        end
        t = step * p.dt
        c = timeit(() -> ED.cut(field, :init))
        cut = c.value
        r = refreshDiscreteOperator!(system, grid, cut)
        cutdata, roots_now, structural = r.cutdata, r.roots, r.structural
        sol = solveElliptic(system)
        push!(steps, (step = step, t = t,
            t_cut = c.time, b_cut = c.bytes,
            t_convert = r.convert.time, b_convert = r.convert.bytes,
            t_aggregate = r.aggregate.time, b_aggregate = r.aggregate.bytes,
            t_refresh = r.refresh.time, b_refresh = r.refresh.bytes,
            t_assemble = sol.stage.assemble.time, b_assemble = sol.stage.assemble.bytes,
            t_factor = sol.stage.factor.time, b_factor = sol.stage.factor.bytes,
            t_solve = sol.stage.solve.time, b_solve = sol.stage.solve.bytes,
            t_recover = sol.stage.recover.time, b_recover = sol.stage.recover.bytes,
            structural = structural, n_master = length(sol.u_m), nnz_B = nnz(sol.B),
            cholesky_ok = sol.cholesky_ok))
        if step % diag_every == 0 || step == p.nsteps
            diag = checkManufacturedSolution(p, grid, cutdata, sol, t;
                system = system, roots = roots_now)
            push!(rows, diag)
            @printf("  n=%d t=%.3f L2rel=%.4e interp=%.4e interior=%.4e bc=%.2e area=%.2e ctr=%.2e n_m=%d agg=%d\n",
                n, t, diag.l2_rel, diag.l2_interp_rel, diag.max_interior, diag.bc_error,
                diag.area_rel_error, diag.center_error, diag.n_master, diag.n_aggregates)
        end
        step == p.nsteps && saveSnapshot(outdir, grid, system, sol, "n$(n)_t500")
    end
    writeRows(joinpath(outdir, "diagnostics_n$(n).csv"), rows)
    writeRows(joinpath(outdir, "steps_n$(n).csv"), steps)

    body0 = steps[2:end]
    nstruct = count(s -> s.structural, body0)
    @printf("stage E n=%d: %d/%d refreshes structural\n", n, nstruct, length(body0))
    for flag in (true, false)
        sel = [s for s in body0 if s.structural == flag]
        isempty(sel) && continue
        @printf("stage E n=%d: refresh_geometry! structural=%-5s n=%3d mean=%6.2f ms max=%6.2f ms alloc mean=%8.1f kB\n",
            n, string(flag), length(sel), 1e3mean(s.t_refresh for s in sel),
            1e3maximum(s.t_refresh for s in sel), mean(s.b_refresh for s in sel) / 1024)
    end
    if stage in ("D", "all")
        # --- aggregate timing over the moving run ---------------------------
        body = steps[2:end]
        function stat(f)
            v = [f(s) for s in body]
            return (sum(v), mean(v), minimum(v), maximum(v))
        end
        for (label, t, b) in (("cut (EvolvingDomains)", stat(s -> s.t_cut), stat(s -> s.b_cut)),
                              ("adapter convert", stat(s -> s.t_convert), stat(s -> s.b_convert)),
                              ("aggregateCutCells", stat(s -> s.t_aggregate), stat(s -> s.b_aggregate)),
                              ("refresh_geometry!", stat(s -> s.t_refresh), stat(s -> s.b_refresh)),
                              ("assemble!", stat(s -> s.t_assemble), stat(s -> s.b_assemble)),
                              ("cholesky", stat(s -> s.t_factor), stat(s -> s.b_factor)),
                              ("solve", stat(s -> s.t_solve), stat(s -> s.b_solve)),
                              ("recover E*u", stat(s -> s.t_recover), stat(s -> s.b_recover)))
            @printf("stage E n=%d: %-26s total=%8.1f ms  mean=%6.2f ms  min=%6.2f  max=%7.2f  alloc mean=%9.1f kB\n",
                n, label, 1e3t[1], 1e3t[2], 1e3t[3], 1e3t[4], b[2] / 1024)
        end
    end

    # --- Stage E: conditioning comparison ----------------------------------
    if stage in ("E", "all")
        cond = conditioningStudy(p, mesh, cutdata0, roots0, form)
        writeRows(joinpath(outdir, "conditioning_n$(n).csv"), [cond])
        @printf("stage E n=%d: cond(B) agfem=%.3e (min eig %.3e, chol ok=%s) vs self-rooted=%.3e (min eig %.3e, chol ok=%s)\n",
            n, cond.cond_agfem, cond.min_eig_agfem, string(cond.cholesky_agfem),
            cond.cond_selfroot, cond.min_eig_selfroot, string(cond.cholesky_selfroot))
        @printf("stage E n=%d: n_m agfem=%d vs self-rooted=%d\n",
            n, cond.n_m_agfem, cond.n_m_selfroot)
    end

    # --- Stage F: Ferrite probe --------------------------------------------
    if stage in ("F", "all")
        probeFerriteQuadrature(p; path = outdir)
    end
    return rows
end

# =============================================================================
# 11. Entry point
# =============================================================================

function main(argv)
    p0 = defaultParams()
    n, nout, stage, out, gamma0 = p0.n, p0.nout, p0.stage, p0.out, p0.gamma0
    i = 1
    while i <= length(argv)
        a = argv[i]
        a == "--n" && (n = parse(Int, argv[i + 1]); i += 2; continue)
        a == "--nout" && (nout = parse(Int, argv[i + 1]); i += 2; continue)
        a == "--stage" && (stage = argv[i + 1]; i += 2; continue)
        a == "--out" && (out = argv[i + 1]; i += 2; continue)
        a == "--gamma0" && (gamma0 = parse(Float64, argv[i + 1]); i += 2; continue)
        error("unknown argument $a")
    end
    return defaultParams(; n = n, nout = nout, out = out, gamma0 = gamma0, stage = stage)
end

if abspath(PROGRAM_FILE) == @__FILE__
    p = main(ARGS)
    configureCutCellsLibrary()
    simulate(p)
end
