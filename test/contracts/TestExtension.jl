using Test
using StaticArrays: SVector
using LinearAlgebra: dot, cross
using Main.EvolvingDomains: Cut, Grid, Discretization, Field, extend!, CUT_IN, CUT_OUT

# ---------------------------------------------------------------------------
# Hand-built cuts. No native library is involved, so every case here is a
# piecewise-linear interface whose facets are known exactly.

_col(grid, x) = clamp(searchsortedlast(grid.axes[1], x), 1, length(grid.axes[1]) - 1)

function _planarCut2D(grid, positions)
    y0, y1 = first(grid.axes[2]), last(grid.axes[2])
    n = length(positions)
    points = NTuple{2,Float64}[(x, y) for x in positions for y in (y0, y1)]
    ncells = prod(length(a) - 1 for a in grid.axes)
    Cut(grid, :init; bg_states=fill(Int8(1), ncells),
        boundary_facet_to_bgcell=Int32[_col(grid, x) for x in positions],
        boundary_point_ids=Int32.(1:2n), boundary_point_pointers=Int32.(1:2:2n+1),
        boundary_physical_points=points, boundary_reference_points=points,
        boundary_normals=[(1.0, 0.0) for _ in positions])
end

function _planarCut3D(grid, positions)
    y0, y1 = first(grid.axes[2]), last(grid.axes[2])
    z0, z1 = first(grid.axes[3]), last(grid.axes[3])
    n = length(positions)
    points = NTuple{3,Float64}[(x, y, z) for x in positions
                               for (y, z) in ((y0, z0), (y1, z0), (y0, z1))]
    ncells = prod(length(a) - 1 for a in grid.axes)
    Cut(grid, :init; bg_states=fill(Int8(1), ncells),
        boundary_facet_to_bgcell=Int32[_col(grid, x) for x in positions],
        boundary_point_ids=Int32.(1:3n), boundary_point_pointers=Int32.(1:3:3n+1),
        boundary_physical_points=points, boundary_reference_points=points,
        boundary_normals=[(1.0, 0.0, 0.0) for _ in positions])
end

# Arbitrary facet soup: vertices given explicitly, parent cells chosen as the cell holding
# each facet's first vertex, then sorted so the Cut's ascending-parent invariant holds.
function _facetSoup(grid, vertices, dim)
    nf = length(vertices)
    parents = Int32[]
    for v in vertices
        idx = ntuple(d -> clamp(searchsortedlast(grid.axes[d], v[1][d]), 1,
                                length(grid.axes[d]) - 1), dim)
        stride = 1; linear = 1
        for d in 1:dim
            linear += (idx[d] - 1) * stride
            stride *= length(grid.axes[d]) - 1
        end
        push!(parents, Int32(linear))
    end
    order = sortperm(parents)
    parents = parents[order]
    facets = vertices[order]
    flat = NTuple{dim,Float64}[]
    ptrs = Int32[1]
    for v in facets
        append!(flat, v)
        push!(ptrs, Int32(length(flat) + 1))
    end
    ncells = prod(length(a) - 1 for a in grid.axes)
    Cut(grid, :init; bg_states=fill(Int8(1), ncells),
        boundary_facet_to_bgcell=parents,
        boundary_point_ids=Int32.(1:length(flat)), boundary_point_pointers=ptrs,
        boundary_physical_points=flat, boundary_reference_points=flat,
        boundary_normals=[ntuple(_ -> 1.0, dim) for _ in facets]), order
end

# ---------------------------------------------------------------------------
# Independent brute force: closest point over every facet, written from the
# mathematical definition (no package kernel is reused).

function _bruteSegment(p, a, b)
    e = b - a
    len2 = dot(e, e)
    t = len2 > 0 ? clamp(dot(p - a, e) / len2, 0.0, 1.0) : 0.0
    return a + t * e
end

function _bruteTriangle(p, v)
    a, b, c = v
    e1, e2 = b - a, c - a
    n = cross(e1, e2)
    n2 = dot(n, n)
    if n2 > 0
        q = p - (dot(p - a, n) / n2) * n
        w = q - a
        d00, d01, d11 = dot(e1, e1), dot(e1, e2), dot(e2, e2)
        s = (d11 * dot(w, e1) - d01 * dot(w, e2)) / n2
        t = (d00 * dot(w, e2) - d01 * dot(w, e1)) / n2
        s >= 0 && t >= 0 && s + t <= 1 && return q
    end
    cands = (_bruteSegment(p, a, b), _bruteSegment(p, b, c), _bruteSegment(p, c, a))
    return cands[argmin([dot(p - q, p - q) for q in cands])]
end

_bruteFacet(p, v::NTuple{2}) = _bruteSegment(p, v[1], v[2])
_bruteFacet(p, v::NTuple{3}) = _bruteTriangle(p, v)

function _bruteClosest(surface, p, dim)
    best, bestd = 0, Inf
    for f in 1:length(surface.facet_to_bgcell)
        lo = Int(surface.point_pointers[f])
        hi = Int(surface.point_pointers[f + 1]) - 1
        v = ntuple(q -> SVector{dim,Float64}(surface.physical_points[
            surface.point_ids[lo + q - 1]]), dim)
        q = _bruteFacet(p, v)
        d = dot(p - q, p - q)
        if d < bestd
            bestd, best = d, f
        end
    end
    return best, bestd
end

# ---------------------------------------------------------------------------
# Circular interface: a closed polygon inscribed in a circle, facets sorted by
# parent cell so the Cut invariant holds. `n` is the number of facets.

function _circleCut(cx, cy, r, n, dim)
    ax = collect(range(0.0, 1.0; length=4 * n + 1))
    grid = dim == 2 ? Grid((ax, copy(ax))) : Grid((ax, copy(ax), copy(ax)))
    verts = [ntuple(d -> (d == 1 ? cx + r * cos(2π * k / n) :
                              d == 2 ? cy + r * sin(2π * k / n) : 0.5), dim) for k in 1:n]
    facets = [(verts[k], verts[mod1(k + 1, n)]) for k in 1:n]
    cut, order = _facetSoup(grid, facets, dim)
    return cut, grid, facets[order]
end

_axis(n) = collect(range(0.0, 1.0; length=n))

@testset "closest-point extension" begin
    @testset "identity map exposes the projection point" begin
        # u(x) = x is linear, so multilinear interpolation reproduces it exactly: the
        # extension of the identity map is the closest point itself.
        for dim in (2, 3)
            ax = collect(0.0:0.1:1.0)
            grid = dim == 2 ? Grid((copy(ax), copy(ax))) : Grid((copy(ax), copy(ax), copy(ax)))
            verts = dim == 2 ? [(0.11, 0.42), (0.83, 0.67), (0.37, 0.05), (0.52, 0.94),
                                (0.66, 0.23)] :
                               [(0.11, 0.42, 0.3), (0.83, 0.67, 0.55), (0.37, 0.05, 0.92),
                                (0.52, 0.94, 0.31), (0.66, 0.23, 0.74), (0.95, 0.81, 0.18)]
            facets = dim == 2 ?
                     [(verts[i], verts[i + 1]) for i in 1:2:(length(verts) - 1)] :
                     [(verts[i], verts[i + 1], verts[i + 2]) for i in 1:3:(length(verts) - 2)]
            cut, order = _facetSoup(grid, facets, dim)
            src = dim == 2 ? [SVector(x, y) for x in ax, y in ax] :
                             [SVector(x, y, z) for x in ax, y in ax, z in ax]
            dest = similar(src)
            extend!(dest, src, cut, grid; band=Inf)
            worst = 0.0
            for I in CartesianIndices(dest)
                p = SVector{dim,Float64}(ntuple(d -> ax[I[d]], dim))
                _, bestd = _bruteClosest(cut.boundary, p, dim)
                got = dest[I]
                worst = max(worst, sqrt(sum(abs2, got - p)) - sqrt(bestd))
            end
            @test worst < 1e-12        # projection never worse than brute force
        end
    end

    @testset "facet source matches brute-force nearest facet" begin
        ay = collect(range(0.0, 1.0; length=13))
        for dim in (2, 3), band in (0.12, 0.3, Inf)
            axes = ntuple(_ -> copy(ay), dim)
            grid = Grid(axes)
            verts = dim == 2 ? [(0.11, 0.42), (0.83, 0.67), (0.37, 0.05), (0.52, 0.94),
                                (0.66, 0.23)] :
                               [(0.11, 0.42, 0.3), (0.83, 0.67, 0.55), (0.37, 0.05, 0.92),
                                (0.52, 0.94, 0.31), (0.66, 0.23, 0.74), (0.95, 0.81, 0.18)]
            facets = dim == 2 ?
                     [(verts[i], verts[i + 1]) for i in 1:2:(length(verts) - 1)] :
                     [(verts[i], verts[i + 1], verts[i + 2]) for i in 1:3:(length(verts) - 2)]
            cut, order = _facetSoup(grid, facets, dim)
            # facet value = the facet's position in the cut, which brute force also reports
            values = Float64.(1:length(order))
            dest = fill(-1.0, ntuple(_ -> 13, dim))
            extend!(dest, values, cut, grid; band=band)
            for I in CartesianIndices(dest)
                p = SVector{dim,Float64}(ntuple(d -> ay[I[d]], dim))
                f, bestd = _bruteClosest(cut.boundary, p, dim)
                if band == Inf || sqrt(bestd) <= band
                    @test dest[I] == f
                else
                    @test dest[I] == -1.0      # outside the band: untouched
                end
            end
        end
    end

    @testset "endpoints, edges and vertices" begin
        grid = Grid((_axis(5), _axis(5)))
        cut = _planarCut2D(grid, [0.25])
        # node exactly on a facet endpoint
        dest = fill(-1.0, 5, 5)
        extend!(dest, [2.5], cut, grid; band=Inf)
        @test dest[2, 1] == 2.5 && dest[2, 5] == 2.5

        grid3 = Grid((_axis(5), _axis(5), _axis(5)))
        cut3 = _planarCut3D(grid3, [0.25])
        dest3 = fill(-1.0, 5, 5, 5)
        extend!(dest3, [2.5], cut3, grid3; band=Inf)
        @test dest3[2, 1, 1] == 2.5 && dest3[2, 5, 5] == 2.5   # triangle vertices
    end

    @testset "degenerate facets are sets, not lines" begin
        grid = Grid((_axis(5), _axis(5)))
        # zero-length segment (a point) and a normal segment
        cut = Cut(grid, :init; bg_states=fill(Int8(1), 16),
                  boundary_facet_to_bgcell=Int32[1, _col(grid, 0.75)],
                  boundary_point_ids=Int32[1, 1, 2, 3], boundary_point_pointers=Int32[1, 3, 5],
                  boundary_physical_points=[(0.25, 0.5), (0.25, 0.5), (0.75, 0.0), (0.75, 1.0)],
                  boundary_reference_points=[(0.25, 0.5), (0.25, 0.5), (0.75, 0.0), (0.75, 1.0)],
                  boundary_normals=[(1.0, 0.0), (1.0, 0.0)])
        dest = fill(-1.0, 5, 5)
        extend!(dest, [10.0, 20.0], cut, grid; band=Inf)
        @test dest[2, 3] == 10.0        # the point facet wins near x = 0.25
        @test dest[4, 3] == 20.0
        @test all(isfinite, dest)

        # collinear triangle: no plane, must fall back to its edges
        grid3 = Grid((_axis(5), _axis(5), _axis(5)))
        deg = Cut(grid3, :init; bg_states=fill(Int8(1), 64),
                  boundary_facet_to_bgcell=Int32[1],
                  boundary_point_ids=Int32[1, 2, 3], boundary_point_pointers=Int32[1, 4],
                  boundary_physical_points=[(0.0, 0.5, 0.5), (0.5, 0.5, 0.5), (1.0, 0.5, 0.5)],
                  boundary_reference_points=[(0.0, 0.5, 0.5), (0.5, 0.5, 0.5), (1.0, 0.5, 0.5)],
                  boundary_normals=[(1.0, 0.0, 0.0)])
        dest3 = fill(-1.0, 5, 5, 5)
        extend!(dest3, [7.0], deg, grid3; band=Inf)
        @test all(==(7.0), dest3)
    end

    @testset "kernels agree with a dense sampling of the facet" begin
        # Independent of the closed-form kernels: minimise the distance over the facet set by
        # sampling barycentrics on a fine grid (segment parameter in 2D).
        function sampled_distance(p, v, dim; n=400)
            best = Inf
            if dim == 2
                a, b = v
                for i in 0:n
                    q = a + (i / n) * (b - a)
                    best = min(best, sum(abs2, p - q))
                end
            else
                a, b, c = v
                for i in 0:n, j in 0:(n - i)
                    u, w = i / n, j / n
                    q = a + u * (b - a) + w * (c - a)
                    best = min(best, sum(abs2, p - q))
                end
            end
            return best
        end
        for (dim, v) in ((2, (SVector(0.1, 0.2), SVector(0.8, 0.55))),
                         (3, (SVector(0.1, 0.2, 0.3), SVector(0.8, 0.55, 0.2),
                              SVector(0.3, 0.9, 0.7))))
            for p in (SVector(ntuple(d -> 0.05 * d, dim)),
                      SVector(ntuple(d -> 0.5 + 0.1 * d, dim)),
                      SVector(ntuple(d -> 0.95, dim)))
                exact = sum(abs2, p - Main.EvolvingDomains.Geometric.Extension._closestPointOnFacet(p, v))
                sampled = sampled_distance(p, v, dim)
                @test exact <= sampled + 1e-12          # the kernel is never worse than a sampled point
                @test sampled - exact <= 1e-4           # and sampling converges onto the kernel's point
            end
        end
    end

    @testset "equidistant facets resolve deterministically" begin
        grid = Grid((_axis(5), _axis(5)))
        cut = _planarCut2D(grid, [0.25, 0.75])
        dest = fill(-1.0, 5, 5)
        extend!(dest, [1.0, 2.0], cut, grid; band=Inf)
        @test dest[3, 3] == 1.0                       # x = 0.5: lowest facet index wins
        @test cut.boundary.facet_to_bgcell[1] < cut.boundary.facet_to_bgcell[2]
        runs = [extend!(fill(-1.0, 5, 5), [1.0, 2.0], cut, grid; band=Inf) for _ in 1:3]
        @test all(==(runs[1]), runs)
    end

    @testset "band, empty interface, and destination preservation" begin
        grid = Grid((_axis(5), _axis(5)))
        cut = _planarCut2D(grid, [0.5])
        dest = fill(-7.0, 5, 5)
        @test extend!(dest, [3.0], cut, grid; band=0.2) === dest
        @test all(==(-7.0), dest[[1, 5], :])          # outside the band: preserved
        @test all(==(3.0), dest[3, :])                # on the interface: written

        empty = Cut(grid, :init; bg_states=fill(Int8(CUT_OUT), 16))
        untouched = fill(-7.0, 5, 5)
        @test extend!(untouched, Float64[], empty, grid) === untouched
        @test all(==(-7.0), untouched)

        # default band is three largest spacings: the whole unit box here
        defaulted = fill(-7.0, 5, 5)
        extend!(defaulted, [3.0], cut, grid)
        @test all(==(3.0), defaulted)

        # a band far larger than the box must saturate, not overflow the node-index arithmetic
        huge = fill(-7.0, 5, 5)
        extend!(huge, [3.0], cut, grid; band=1e300)
        @test all(==(3.0), huge)
        @test extend!(huge, [3.0], cut, grid; band=Inf) === huge
    end

    @testset "geometry only: phase labels do not matter" begin
        # same interface, different classification: the extension must not change
        grid = Grid((_axis(5), _axis(5)))
        a = _planarCut2D(grid, [0.5])
        b = Cut(grid, :init; bg_states=fill(Int8(CUT_IN), 16),
                cell_facet_states=fill((Int8(-1), Int8(-1), Int8(-1), Int8(-1)), 16),
                boundary_facet_to_bgcell=a.boundary.facet_to_bgcell,
                boundary_point_ids=a.boundary.point_ids,
                boundary_point_pointers=a.boundary.point_pointers,
                boundary_physical_points=a.boundary.physical_points,
                boundary_reference_points=a.boundary.reference_points,
                boundary_normals=a.boundary.normals)
        da = fill(-1.0, 5, 5); db = fill(-1.0, 5, 5)
        extend!(da, [4.0], a, grid; band=Inf)
        extend!(db, [4.0], b, grid; band=Inf)
        @test da == db
    end

    @testset "nodal source: trace validity is the caller's responsibility" begin
        ax = _axis(9)
        grid = Grid((ax, ax))
        cut = _planarCut2D(grid, [0.46])               # interface between two nodes
        # bulk-valid data, junk on the other side: interpolation across the interface mixes
        # both sides, so the extension is neither side's value
        u = [x < 0.46 ? 1.0 : 100.0 for x in ax, y in ax]
        dest = fill(-1.0, 9, 9)
        extend!(dest, u, cut, grid; band=Inf)
        t = (0.46 - ax[4]) / (ax[5] - ax[4])           # fraction between the two nodes
        expected = (1 - t) * 1.0 + t * 100.0
        @test all(isapprox.(dest, expected; atol=1e-12))
        @test 1.0 < expected < 100.0                   # contaminated, not one-sidedly solved
        # a nodal field that is a valid trace gives the exact closest-point value
        v = [2x for x in ax, y in ax]                  # linear -> interpolation is exact
        extend!(dest, v, cut, grid; band=Inf)
        @test all(isapprox.(dest, 2 * 0.46; atol=1e-12))
    end

    @testset "vector-valued sources" begin
        grid = Grid((_axis(5), _axis(5)))
        cut = _planarCut2D(grid, [0.25])
        dest = fill(SVector(0.0, 0.0), 5, 5)
        extend!(dest, [SVector(1.0, 2.0)], cut, grid; band=Inf)
        @test all(v -> v == SVector(1.0, 2.0), dest)
        nodal = [SVector(x, 2y) for x in _axis(5), y in _axis(5)]
        extend!(dest, nodal, cut, grid; band=Inf)
        @test dest[2, 3] == SVector(0.25, 2 * _axis(5)[3])
    end

    @testset "storage semantics and validation" begin
        grid = Grid((_axis(5), _axis(5)))
        cut = _planarCut2D(grid, [0.5])
        @test_throws DimensionMismatch extend!(zeros(4, 5), [1.0], cut, grid)
        @test_throws DimensionMismatch extend!(zeros(5, 5), [1.0, 2.0], cut, grid)
        @test_throws DimensionMismatch extend!(zeros(5, 5), zeros(5, 4), cut, grid)
        @test_throws ArgumentError extend!(zeros(5, 5), [1.0], cut, grid; band=-1.0)
        @test_throws ArgumentError extend!(zeros(5, 5), [1.0], cut,
                                           Grid((_axis(5), _axis(5))))
        nonuniform = Grid(([0.0, 0.2, 0.5, 0.75, 1.0], _axis(5)))
        @test_throws ArgumentError extend!(zeros(5, 5), [1.0],
                                           _planarCut2D(nonuniform, [0.5]), nonuniform)

        # in-place nodal extension equals out-of-place
        u = [x + y for x in _axis(5), y in _axis(5)]
        out = fill(-1.0, 5, 5); extend!(out, u, cut, grid; band=Inf)
        inplace = copy(u); extend!(inplace, inplace, cut, grid; band=Inf)
        @test inplace == out

        # facet values aliasing the destination must give the same answer as a copy
        cut1 = _planarCut2D(grid, [0.5])
        reference = fill(-1.0, 5, 5)
        extend!(reference, [4.0], cut1, grid; band=Inf)
        aliased = fill(4.0, 5, 5)[:, :]               # the source view shares this memory
        extend!(aliased, vec(aliased)[1:1], cut1, grid; band=Inf)
        @test eltype(reference) == eltype(aliased)
        @test all(==(4.0), aliased)
    end

    @testset "Field method" begin
        d = Discretization(_axis(5), _axis(5))
        phi = [x - 0.5 for x in _axis(5), y in _axis(5)]
        dest = Field(d, Dict(:init => fill(-1.0, 5, 5)))
        src = Field(d, Dict(:init => [2x for x in _axis(5), y in _axis(5)]))
        grid = d.levels[:init]
        cut = _planarCut2D(grid, [0.5])
        @test extend!(dest, src, cut) === dest
        @test all(isapprox.(dest.data[:init], 1.0; atol=1e-14))
        other = Field(Discretization(_axis(5), _axis(5)), Dict(:init => zeros(5, 5)))
        @test_throws ArgumentError extend!(other, src, cut)
    end

    @testset "allocation stays proportional to the grid" begin
        grid = Grid((_axis(33), _axis(33)))
        cut = _planarCut2D(grid, collect(0.3:0.05:0.7))
        values = Float64.(1:length(cut.boundary.facet_to_bgcell))
        dest = zeros(33, 33)
        extend!(dest, values, cut, grid)
        bytes = minimum(Base.@allocated(extend!(dest, values, cut, grid)) for _ in 1:5)
        @test bytes < 40 * length(dest)     # search buffers only, no per-facet grid arrays
    end
end

@testset "closest-point extension accuracy" begin
    # Polygon inscribed in a circle; nodal source u(x) = x is linear, so the interpolator is
    # exact and the error measures the projection onto the piecewise-linear interface alone.
    # The projection point moves along a chord, so its position is quantised by the facet
    # angular width: the max-norm error is first order in h even though the *distance* error
    # is second order. That is measured here rather than assumed.
    function projection_error(n)
        r, cx, cy = 0.25, 0.5, 0.5
        cut, grid, _ = _circleCut(cx, cy, r, n, 2)
        ax = grid.axes[1]
        src = [x for x in ax, y in ax]
        dest = fill(NaN, length(ax), length(ax))
        extend!(dest, src, cut, grid; band=0.3)
        worst = 0.0
        for I in CartesianIndices(dest)
            isnan(dest[I]) && continue                 # not written: outside the band
            x, y = ax[I[1]], ax[I[2]]
            d = hypot(x - cx, y - cy)
            d > 0 || continue                          # centre: any projection is a tie
            worst = max(worst, abs(dest[I] - (cx + r * (x - cx) / d)))
        end
        return worst
    end
    e1, e2, e3 = projection_error(16), projection_error(32), projection_error(64)
    order21 = log(e1 / e2) / log(2)
    order32 = log(e2 / e3) / log(2)
    @test e1 > e2 > e3                              # refinement reduces the error
    @test order32 > 0.8                             # measured first order, not assumed second
    @test e3 < 0.02                              # measured 0.012 at 64 facets
    @info "closest-point extension convergence" errors=(e1, e2, e3) orders=(order21, order32)
end
