using Test
using StaticArrays: SVector
using LinearAlgebra: dot, cross
using Main.EvolvingDomains: Cut, Grid, extendCutValues, extendCutValues!

# Hand-built boundary-only cuts: every facet spans the whole off-axis extent of a plane
# perpendicular to axis 1, at the given positions, in ascending order. Facets therefore sit
# in the cut's own order (`cut.boundary.facet_to_bgcell` ascending), which is the order the
# per-facet values are given in — `extendCutValues` needs no extra mapping.
function _column(grid, x)
    axis = grid.axes[1]
    clamp(searchsortedlast(axis, x), 1, length(axis) - 1)
end

function _planarCut2D(grid, positions)
    y0, y1 = first(grid.axes[2]), last(grid.axes[2])
    n = length(positions)
    points = NTuple{2,Float64}[(x, y) for x in positions for y in (y0, y1)]
    ncells = (length(grid.axes[1]) - 1) * (length(grid.axes[2]) - 1)
    Cut(grid, :init; bg_states=fill(Int8(1), ncells),
        boundary_facet_to_bgcell=Int32[_column(grid, x) for x in positions],
        boundary_point_ids=Int32.(1:2n), boundary_point_pointers=Int32.(1:2:2n+1),
        boundary_physical_points=points, boundary_reference_points=points,
        boundary_normals=[(1.0, 0.0) for _ in positions])
end

function _planarCut3D(grid, positions)    y0, y1 = first(grid.axes[2]), last(grid.axes[2])
    z0, z1 = first(grid.axes[3]), last(grid.axes[3])
    n = length(positions)
    points = NTuple{3,Float64}[(x, y, z) for x in positions
                               for (y, z) in ((y0, z0), (y1, z0), (y0, z1))]
    ncells = prod(length(axis) - 1 for axis in grid.axes)
    Cut(grid, :init; bg_states=fill(Int8(1), ncells),
        boundary_facet_to_bgcell=Int32[_column(grid, x) for x in positions],
        boundary_point_ids=Int32.(1:3n), boundary_point_pointers=Int32.(1:3:3n+1),
        boundary_physical_points=points, boundary_reference_points=points,
        boundary_normals=[(1.0, 0.0, 0.0) for _ in positions])
end

_whole(row, value) = all(isequal(value), row)

# Independent brute-force nearest-facet search: every node against every facet, same
# definition of distance (segment in 2D, triangle in 3D). Irrelevant where a facet is
# axis-aligned, which is why it exercises slanted facets.
function _bruteForce(surface, grid, band)
    d = length(grid.axes)
    origin = ntuple(k -> first(grid.axes[k]), d)
    h = ntuple(k -> (last(grid.axes[k]) - first(grid.axes[k])) / (length(grid.axes[k]) - 1), d)
    dims = ntuple(k -> length(grid.axes[k]), d)
    out = Array{Int}(undef, dims)
    for I in CartesianIndices(dims)
        p = SVector{d,Float64}(ntuple(k -> origin[k] + (I[k] - 1) * h[k], d))
        best, bestd = 0, Inf
        for f in 1:length(surface.facet_to_bgcell)
            v = [SVector{d,Float64}(surface.physical_points[surface.point_ids[i]])
                 for i in Int(surface.point_pointers[f]):Int(surface.point_pointers[f + 1]) - 1]
            dd = d == 2 ? _bruteSegment(p, v[1], v[2]) : _bruteTriangle(p, v)
            if dd < bestd
                bestd, best = dd, f
            end
        end
        out[I] = best > 0 && bestd <= band^2 ? best : 0
    end
    return out
end

function _bruteSegment(p, a, b)
    e = b - a
    len2 = dot(e, e)
    t = len2 > 0 ? clamp(dot(p - a, e) / len2, 0.0, 1.0) : 0.0
    return sum(abs2, p - (a + t * e))
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
        s >= 0 && t >= 0 && s + t <= 1 && return sum(abs2, p - q)
    end
    return min(_bruteSegment(p, a, b), _bruteSegment(p, b, c), _bruteSegment(p, c, a))
end

# Interior points of the unit cube, used as slanted facets of both dimensions.
const _SLANT_2D = [(0.11, 0.42), (0.83, 0.67), (0.37, 0.05), (0.52, 0.94),
                   (0.66, 0.23), (0.95, 0.81)]
const _SLANT_3D = [(0.11, 0.42, 0.03), (0.83, 0.67, 0.55), (0.37, 0.05, 0.92),
                   (0.52, 0.94, 0.31), (0.66, 0.23, 0.74), (0.95, 0.81, 0.18)]

@testset "nearest-facet cut-value extension" begin
    axis = [0.0, 0.25, 0.5, 0.75, 1.0]
    grid = Grid((copy(axis), copy(axis)))
    grid3 = Grid((copy(axis), copy(axis), copy(axis)))
    xs = grid.axes[1]

    @testset "constant facets extend the constant" begin
        cut = _planarCut2D(grid, [0.25])
        dest = extendCutValues([3.0], cut, grid; band=Inf)
        @test size(dest) == (5, 5) && eltype(dest) === Float64
        @test all(==(3.0), dest)
        @test all(==(3.0), vec(dest))            # transport operators take the linear view
    end

    @testset "band limits the extension and fillvalue fills the rest" begin
        cut = _planarCut2D(grid, [0.25])
        dest = extendCutValues([2.0], cut, grid; band=0.3, fillvalue=-1.0)
        for i in eachindex(xs)
            @test _whole(view(dest, i, :), abs(xs[i] - 0.25) <= 0.3 ? 2.0 : -1.0)
        end
        # The default band is three of the largest node spacing: every node is covered here.
        @test all(==(2.0), extendCutValues([2.0], cut, grid))
        # Vector-valued per-facet data keeps its type and its components.
        vectors = extendCutValues([SVector(2.0, -1.0)], cut, grid; band=Inf)
        @test eltype(vectors) === SVector{2,Float64}
        @test all(v -> v == SVector(2.0, -1.0), vectors)
    end

    @testset "each node carries its closest facet's value" begin
        cut = _planarCut2D(grid, [0.25, 0.75])
        dest = extendCutValues([1.0, 2.0], cut, grid; band=Inf)
        for i in eachindex(xs)
            expected = xs[i] > 0.5 ? 2.0 : 1.0        # x = 0.5 is a tie: lowest index wins
            @test _whole(view(dest, i, :), expected)
        end
        vectors = extendCutValues([SVector(1.0, 2.0), SVector(-1.0, 0.5)], cut, grid; band=Inf)
        for i in eachindex(xs)
            expected = xs[i] > 0.5 ? SVector(-1.0, 0.5) : SVector(1.0, 2.0)
            @test all(v -> v == expected, view(vectors, i, :))
        end
        # Equidistant facets: the lowest facet index wins, deterministically and repeatably.
        first_run = extendCutValues([1.0, 2.0], cut, grid; band=Inf)
        @test all(extendCutValues([1.0, 2.0], cut, grid; band=Inf) == first_run for _ in 1:3)
        @test xs[3] == 0.5 && _whole(view(first_run, 3, :), 1.0)
    end

    @testset "3D facets and 3D values" begin
        cut3 = _planarCut3D(grid3, [0.25, 0.75])
        scalars = extendCutValues([4.0, -4.0], cut3, grid3; band=Inf)
        @test size(scalars) == (5, 5, 5)
        for i in eachindex(axis)
            @test _whole(view(scalars, i, :, :), axis[i] > 0.5 ? -4.0 : 4.0)
        end
        vectors = extendCutValues([SVector(1.0, 2.0, 3.0), SVector(-1.0, 0.0, 1.0)], cut3,
                                  grid3; band=Inf)
        @test eltype(vectors) === SVector{3,Float64}
        @test all(v -> v == SVector(1.0, 2.0, 3.0), view(vectors, 1, :, :))
        @test all(v -> v == SVector(-1.0, 0.0, 1.0), view(vectors, 5, :, :))
        # A facet lying in a node plane: nodes on the facet get its value at zero distance.
        dest = Array{SVector{3,Float64}}(undef, 5, 5, 5)
        @test extendCutValues!(dest, [SVector(1.0, 2.0, 3.0)], _planarCut3D(grid3, [0.5]),
                               grid3; band=Inf) === dest
        @test all(v -> v == SVector(1.0, 2.0, 3.0), dest)
    end

    @testset "slanted facets match a brute-force search" begin
        ay = collect(range(0.0, 1.0; length=13))
        for (dim, verts) in ((2, _SLANT_2D), (3, _SLANT_3D))
            axes = ntuple(_ -> copy(ay), dim)
            g = Grid(axes)
            nf = length(verts) ÷ dim
            npts = dim * nf
            pts = NTuple{dim,Float64}[(verts[dim * (f - 1) + q]) for f in 1:nf for q in 1:dim]
            ncells = prod(length(axis) - 1 for axis in axes)
            cut = Cut(g, :init; bg_states=fill(Int8(1), ncells),
                      boundary_facet_to_bgcell=Int32.(1:nf),
                      boundary_point_ids=Int32.(1:npts),
                      boundary_point_pointers=Int32.(1:dim:(npts + 1)),
                      boundary_physical_points=pts, boundary_reference_points=pts,
                      boundary_normals=[ntuple(_ -> 1.0, dim) for _ in 1:nf])
            for band in (0.12, 0.3, Inf)
                values = Float64.(1:nf)
                got = extendCutValues(values, cut, g; band=band)
                expected = _bruteForce(cut.boundary, g, band)
                @test all(I -> got[I] == (expected[I] == 0 ? 0.0 : expected[I]),
                          CartesianIndices(size(got)))
            end
        end
    end

    @testset "no interface, validation, and in-place allocation" begin
        lone = Grid((copy(axis), copy(axis)))
        @test all(==(0.0), extendCutValues(Float64[], _planarCut2D(lone, Float64[]), lone))
        cut = _planarCut2D(grid, [0.25])
        @test_throws DimensionMismatch extendCutValues([1.0, 2.0], cut, grid)  # facet count
        @test_throws DimensionMismatch extendCutValues!(zeros(4, 5), [1.0], cut, grid)
        @test_throws ArgumentError extendCutValues([1.0], cut, grid; band=-1.0)
        @test_throws ArgumentError extendCutValues([1.0], cut, Grid((copy(axis), copy(axis))))

        # Repeated in-place extension must not rebuild a grid-sized destination of values:
        # what it allocates is the nearest-facet record, one Float64 per node.
        values = [SVector(1.0, 2.0, 3.0), SVector(-1.0, 0.0, 1.0)]
        cut3 = _planarCut3D(grid3, [0.25, 0.75])
        dest = Array{SVector{3,Float64}}(undef, 5, 5, 5)
        extendCutValues!(dest, values, cut3, grid3)
        bytes = minimum(Base.@allocated(extendCutValues!(dest, values, cut3, grid3)) for _ in 1:6)
        @test bytes < prod(size(dest)) * sizeof(eltype(dest))
    end
end
