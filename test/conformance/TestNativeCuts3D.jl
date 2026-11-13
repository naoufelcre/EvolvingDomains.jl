using Test
using Main.EvolvingDomains: Grid, Discretization, Field, cut, subcells, cut_cell_ids

function _tetra_volume3(cut, cell)
    ids = Int(cut.point_pointers[cell]):Int(cut.point_pointers[cell + 1])-1
    p, q, r, s = cut.physical_points[Int.(cut.point_ids[ids])]
    a, b, c = ntuple(d -> q[d]-p[d], 3), ntuple(d -> r[d]-p[d], 3), ntuple(d -> s[d]-p[d], 3)
    return abs(a[1]*(b[2]*c[3]-b[3]*c[2])-a[2]*(b[1]*c[3]-b[3]*c[1])+a[3]*(b[1]*c[2]-b[2]*c[1]))/6
end

_triangle_points3(cut, face) = cut.boundary.physical_points[
    Int.(cut.boundary.point_ids[Int(cut.boundary.point_pointers[face]):Int(cut.boundary.point_pointers[face + 1])-1])]

# cross product of the first two edges: the normal implied by the winding
function _triangle_cross3(points)
    p, q, r = points
    a, b = ntuple(d -> q[d]-p[d], 3), ntuple(d -> r[d]-p[d], 3)
    return (a[2]*b[3]-a[3]*b[2], a[3]*b[1]-a[1]*b[3], a[1]*b[2]-a[2]*b[1])
end

_triangle_area3(points) = hypot(_triangle_cross3(points)...) / 2
_triangle_cross3(cut, face) = _triangle_cross3(_triangle_points3(cut, face))
_triangle_area3(cut, face) = _triangle_area3(_triangle_points3(cut, face))


@inline _dot3(a, b) = a[1]*b[1] + a[2]*b[2] + a[3]*b[3]

# every emitted triangle must wind so that its cross product follows the stored normal
_winding_matches_normal(cut) = all(
    f -> _dot3(_triangle_cross3(cut, f), cut.boundary.normals[f]) > 0,
    eachindex(cut.boundary.facet_to_bgcell))

# Non-negligible facet normals are the unit plane normals, not unrelated gradients.
function _normal_matches_cross(cut)
    for f in eachindex(cut.boundary.facet_to_bgcell)
        area = _triangle_area3(cut, f)
        area > 1e-8 || continue
        cx, cy, cz = _triangle_cross3(cut, f)
        n = cut.boundary.normals[f]
        isapprox(n[1], cx / (2area); atol=1e-10) &&
            isapprox(n[2], cy / (2area); atol=1e-10) &&
            isapprox(n[3], cz / (2area); atol=1e-10) || return false
    end
    return true
end

@testset "3D native CutCells backend" begin
    axes = ([2.0,3.0], [-1.0,2.0], [5.0,7.0])
    grid = Grid(axes)
    phi = [x - 2.5 for x in axes[1], y in axes[2], z in axes[3]]
    field = Field(Discretization(Dict(:plane => grid)), Dict(:plane => phi))
    snapshot = cut(field, :plane)
    @test !isempty(snapshot.subcell_states)
    @test all(diff(snapshot.point_pointers) .== 4)
    @test Set(snapshot.subcell_states) == Set(Int8[-1, 1])
    @test isapprox(sum(_tetra_volume3(snapshot, c) for c in eachindex(snapshot.subcell_states)), 6.0; atol=1e-12)
    @test isapprox(sum(_tetra_volume3(snapshot, c) for c in eachindex(snapshot.subcell_states)
                       if snapshot.subcell_states[c] == -1), 3.0; atol=1e-12)
    @test isapprox(sum(_tetra_volume3(snapshot, c) for c in eachindex(snapshot.subcell_states)
                       if snapshot.subcell_states[c] == 1), 3.0; atol=1e-12)
    @test !isempty(snapshot.boundary.facet_to_bgcell)
    @test all(diff(snapshot.boundary.point_pointers) .== 3)
    @test _winding_matches_normal(snapshot)
    @test isapprox(sum(_triangle_area3(snapshot, f) for f in eachindex(snapshot.boundary.facet_to_bgcell)), 6.0; atol=1e-12)
    @test all(isapprox.(first.(snapshot.boundary.normals), 1.0; atol=1e-12))
    @test all(isapprox.(getindex.(snapshot.boundary.normals, 2), 0.0; atol=1e-12))
    @test all(isapprox.(getindex.(snapshot.boundary.normals, 3), 0.0; atol=1e-12))
    @test all(p -> all(0.0 .<= p .<= 1.0), snapshot.reference_points)
    @test all(p -> all(isfinite, p), snapshot.physical_points)
    @test all(isapprox(snapshot.physical_points[q][1], 2.0 + p[1]) &&
              isapprox(snapshot.physical_points[q][2], -1.0 + 3p[2]) &&
              isapprox(snapshot.physical_points[q][3], 5.0 + 2p[3])
              for q in eachindex(snapshot.reference_points) for p in (snapshot.reference_points[q],))
    @test all(isapprox(snapshot.boundary.physical_points[q][1], 2.5) &&
              isapprox(snapshot.boundary.physical_points[q][2], -1.0 + 3snapshot.boundary.reference_points[q][2]) &&
              isapprox(snapshot.boundary.physical_points[q][3], 5.0 + 2snapshot.boundary.reference_points[q][3])
              for q in eachindex(snapshot.boundary.reference_points))
    @test all(snapshot.cell_facet_states[c][2] == snapshot.cell_facet_states[c+1][1]
              for c in 1:(length(snapshot.cell_facet_states)-1))
    # precomputed background-cell -> subcell index against a fresh scan
    @test cut_cell_ids(snapshot) == Int32.(findall(==(Int8(0)), snapshot.bg_states))
    @test sum(c -> length(subcells(snapshot, c)), 1:length(snapshot.bg_states)) ==
          length(snapshot.subcell_states)
    @test all(c -> collect(subcells(snapshot, c)) ==
                   Int32.(findall(==(Int32(c)), snapshot.cell_to_bgcell)),
              1:length(snapshot.bg_states))

    aligned_axes = ([0.0,0.5,1.0], [0.0,1.0], [0.0,1.0])
    aligned_grid = Grid(aligned_axes)
    aligned_phi = [x - 0.5 for x in aligned_axes[1], y in aligned_axes[2], z in aligned_axes[3]]
    aligned_field = Field(Discretization(Dict(:aligned => aligned_grid)), Dict(:aligned => aligned_phi))
    aligned = cut(aligned_field, :aligned)
    @test aligned.bg_states == Int8[0,1]
    @test isapprox(sum(_tetra_volume3(aligned, c) for c in eachindex(aligned.subcell_states)
                       if aligned.subcell_states[c] == -1), 0.5; atol=1e-12)
    @test isapprox(sum(_tetra_volume3(aligned, c) for c in eachindex(aligned.subcell_states)
                       if aligned.subcell_states[c] == 1), 0.5; atol=1e-12)
    @test !isempty(aligned.boundary.facet_to_bgcell)
    @test all(==(Int32(1)), aligned.boundary.facet_to_bgcell)
    @test isapprox(sum(_triangle_area3(aligned, f) for f in eachindex(aligned.boundary.facet_to_bgcell)), 1.0; atol=1e-12)
    @test all(isapprox(n[1], 1.0; atol=1e-12) for n in aligned.boundary.normals)
    @test _winding_matches_normal(aligned)
    facets = aligned.cell_facet_states
    @test facets[1][2] == facets[2][1] == Int8(0)

    # z-oriented interface: winding and stored normal follow the level-set gradient
    zaxes = ([0.0,1.0], [0.0,1.0], [0.0,0.5,1.0])
    zphi = [z - 0.25 for x in zaxes[1], y in zaxes[2], z in zaxes[3]]
    zfield = Field(Discretization(Dict(:z => Grid(zaxes))), Dict(:z => zphi))
    zcut = cut(zfield, :z)
    @test !isempty(zcut.boundary.facet_to_bgcell)
    @test all(all(isapprox(n[d], d == 3 ? 1.0 : 0.0; atol=1e-12) for d in 1:3) for n in zcut.boundary.normals)
    @test _winding_matches_normal(zcut)
    @test isapprox(sum(_triangle_area3(zcut, f) for f in eachindex(zcut.boundary.facet_to_bgcell)), 1.0; atol=1e-12)

    # oblique plane: normal is the normalized diagonal
    daxes = ([0.0,0.5,1.0], [0.0,0.5,1.0], [0.0,0.5,1.0])
    dphi = [x + y + z - 1.0 for x in daxes[1], y in daxes[2], z in daxes[3]]
    dfield = Field(Discretization(Dict(:diag => Grid(daxes))), Dict(:diag => dphi))
    dcut = cut(dfield, :diag)
    @test !isempty(dcut.boundary.facet_to_bgcell)
    @test _winding_matches_normal(dcut)
    @test _normal_matches_cross(dcut)
    # Sliver facets have ill-conditioned cross products; check every non-negligible facet.
    for f in eachindex(dcut.boundary.facet_to_bgcell)
        _triangle_area3(dcut, f) > 1e-8 || continue
        n = dcut.boundary.normals[f]
        @test all(isapprox(n[d], inv(sqrt(3.0)); atol=1e-10) for d in 1:3)
    end

    # non-planar interface: every triangle winds outwards from the negative region
    saxes = (collect(0.0:0.25:1.0), collect(0.0:0.25:1.0), collect(0.0:0.25:1.0))
    sphi = [(x-0.5)^2 + (y-0.5)^2 + (z-0.5)^2 - 0.16 for x in saxes[1], y in saxes[2], z in saxes[3]]
    sfield = Field(Discretization(Dict(:sphere => Grid(saxes))), Dict(:sphere => sphi))
    scut = cut(sfield, :sphere)
    @test !isempty(scut.boundary.facet_to_bgcell)
    @test _winding_matches_normal(scut)
    @test _normal_matches_cross(scut)
    @test all(f -> begin
        ids = Int(scut.boundary.point_pointers[f]):Int(scut.boundary.point_pointers[f + 1])-1
        p, q, r = scut.boundary.physical_points[Int.(scut.boundary.point_ids[ids])]
        radial = ntuple(d -> (p[d] + q[d] + r[d]) / 3 - 0.5, 3)
        radius = hypot(radial...)
        radius > 0 ? _dot3(ntuple(d -> radial[d] / radius, 3), scut.boundary.normals[f]) > 0.9 : false
    end, eachindex(scut.boundary.facet_to_bgcell))
    @test isapprox(sum(_tetra_volume3(scut, c) for c in eachindex(scut.subcell_states)), 1.0; atol=1e-12)

    allzero = Field(Discretization(Dict(:zero => Grid(([0.0,1.0],[0.0,1.0],[0.0,1.0])))),
        Dict(:zero => zeros(2,2,2)))
    @test_throws ArgumentError cut(allzero, :zero)
end

@testset "3D interface winding follows the level-set gradient" begin
    binding = Main.EvolvingDomains.CutCellsBinding
    axes = ([0.0,1.0], [0.0,1.0], [0.0,1.0])
    values = (-0.5,-0.5,0.5,0.5,-0.5,-0.5,0.5,0.5)   # phi = y - 0.5, gradient +y
    a, b, c = (0.0,0.5,0.0), (1.0,0.5,0.0), (0.0,0.5,1.0)
    emit = function (triangle)
        physical, reference = NTuple{3,Float64}[], NTuple{3,Float64}[]
        parents, ids, pointers, normals = Int32[], Int32[], Int32[1], NTuple{3,Float64}[]
        binding._append3DBoundary!(physical, reference, parents, ids, pointers, normals,
            triangle, 1, values, axes, (1,1,1); orient_by_levelset=true)
        return physical, ids, pointers, normals
    end
    ccw = (collect(Float64, (a..., b..., c...)), Int32[0,3], Int32[0,1,2], Int32[2])
    cw = (collect(Float64, (a..., c..., b...)), Int32[0,3], Int32[0,1,2], Int32[2])
    degenerate = (collect(Float64, (a..., b..., b...)), Int32[0,3], Int32[0,1,2], Int32[2])
    for triangle in (ccw, cw)   # the winding swap must repair either input order
        physical, ids, pointers, normals = emit(triangle)
        @test length(normals) == 1 && pointers == Int32[1,4] && length(ids) == 3
        @test normals[1] == (0.0, 1.0, 0.0)
        points = physical[Int.(ids)]
        @test _dot3(_triangle_cross3(points), normals[1]) > 0
        @test isapprox(_triangle_area3(points), 0.5; atol=1e-12)
    end
    physical, ids, pointers, normals = emit(degenerate)
    @test isempty(ids) && pointers == Int32[1] && isempty(normals)
end

@testset "3D facets are physically nondegenerate and physically normalised" begin
    # A sphere whose poles land on node planes makes the cutting plane tangent to a Kuhn
    # tet: the triangle collapses during the reference -> physical map, so a degeneracy
    # check done in reference coordinates misses it entirely.
    ax = collect(range(0.0, 1.0; length=17))
    g = Grid((ax, copy(ax), copy(ax)))
    phi = [hypot(x - 0.5, y - 0.5, z - 0.5) - 0.25 for x in ax, y in ax, z in ax]
    s = cut(Field(Discretization(Dict(:init => g)), Dict(:init => phi)), :init)
    @test !isempty(s.boundary.facet_to_bgcell)
    @test all(f -> _triangle_area3(s, f) > 0, eachindex(s.boundary.facet_to_bgcell))

    # Unequal spacings: the stored normal must be the physical one. Computed in reference
    # coordinates it is the normal of a differently-shaped triangle (measured 33.7 deg off).
    axx = collect(range(0.0, 1.0; length=11))    # hx = 0.1
    axy = collect(range(0.0, 0.4; length=11))    # hy = 0.04
    g2 = Grid((axx, axy, copy(axy)))
    plane = [x + y - 0.5 for x in axx, y in axy, z in axy]
    s2 = cut(Field(Discretization(Dict(:init => g2)), Dict(:init => plane)), :init)
    expected = (1 / sqrt(2), 1 / sqrt(2), 0.0)
    @test !isempty(s2.boundary.normals)
    # Angular criterion: a reference-space normal is 33.7 deg off here, so 1e-6 rad still
    # separates the two by ten orders of magnitude. Thin facets carry ~1e-12 direction noise
    # from rounding of their own vertices, which is why this is not tested bitwise.
    @test all(n -> sum(n[d] * expected[d] for d in 1:3) > 1 - 1e-12, s2.boundary.normals)
end
