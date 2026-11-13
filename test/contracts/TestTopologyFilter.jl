using Test
using Main.EvolvingDomains
using Main.EvolvingDomains.Geometric.TopologyFilter: filterSmallPhaseIslands!

# Uniform grid over the unit box: 9 nodes per axis give spacing exactly 1/8.
function _island_field(dims, phi)
    axes = ntuple(d -> collect(range(0.0, 1.0; length=dims[d])), length(dims))
    grid = EvolvingDomains.Grid(axes)
    return EvolvingDomains.Field(EvolvingDomains.Discretization(Dict(:g => grid)),
                                 Dict(:g => phi))
end

@testset "small-phase-island filter" begin
    # An isolated 4-node blob is flipped; a 9-node region of the same phase survives.
    phi = ones(9, 9)
    phi[4:5, 4:5] .= -1.0
    phi[7:9, 1:3] .= -1.0
    stats = filterSmallPhaseIslands!(phi; min_nodes=5)
    @test stats == (n_components=2, n_removed_components=1, n_flipped_nodes=4)
    @test all(==(1.0), phi[4:5, 4:5])              # magnitude defaults to one(T)
    @test all(==(-1.0), phi[7:9, 1:3])
    @test count(==(1.0), phi) == 81 - 9

    # The threshold is strict: exactly `min_nodes` nodes is a kept component.
    phi = ones(5, 5)
    phi[2:4, 3] .= -1.0
    phi[2, 2] = phi[2, 4] = -1.0
    before = copy(phi)
    stats = filterSmallPhaseIslands!(phi; min_nodes=5)
    @test stats.n_removed_components == 0 && stats.n_flipped_nodes == 0
    @test phi == before

    # Zero belongs to the negative phase only, so a positive filter never sees it.
    phi = fill(-1.0, 9, 9)
    phi[5, 5] = 1.0
    phi[3, 3] = 0.0
    stats = filterSmallPhaseIslands!(phi; phase=:positive, min_nodes=2, magnitude=0.5)
    @test stats == (n_components=1, n_removed_components=1, n_flipped_nodes=1)
    @test phi[5, 5] == -0.5 && phi[3, 3] == 0.0

    # A single-signed field is one component of its own phase: nothing is small.
    zeros_phi = zeros(9, 9)
    @test filterSmallPhaseIslands!(zeros_phi; min_nodes=5) ==
          (n_components=1, n_removed_components=0, n_flipped_nodes=0)

    # 3D uses six face neighbours: a diagonal pair is two components, not one.
    phi3 = ones(5, 5, 5)
    phi3[2, 2, 2] = -1.0
    phi3[4, 4, 4] = -1.0
    stats = filterSmallPhaseIslands!(phi3; min_nodes=2, magnitude=0.25)
    @test stats == (n_components=2, n_removed_components=2, n_flipped_nodes=2)
    @test count(==(0.25), phi3) == 2 && phi3[2, 2, 2] == phi3[4, 4, 4] == 0.25
    phi3 = ones(5, 5, 5)
    phi3[2, 2, 2] = phi3[3, 2, 2] = -1.0
    @test filterSmallPhaseIslands!(phi3; min_nodes=2).n_removed_components == 0  # not < 2

    # The Field method defaults the flip magnitude to the node spacing.
    phi9 = ones(9, 9)
    phi9[4, 4] = -1.0
    field = _island_field((9, 9), phi9)
    stats = filterSmallPhaseIslands!(field, :g)
    @test stats.n_flipped_nodes == 1 && phi9[4, 4] == 1 / 8

    @test_throws ArgumentError filterSmallPhaseIslands!(ones(3, 3); phase=:both)
    @test_throws ArgumentError filterSmallPhaseIslands!(ones(3, 3); min_nodes=0)
    @test_throws ArgumentError filterSmallPhaseIslands!(ones(3, 3); magnitude=0.0)
    @test_throws ArgumentError filterSmallPhaseIslands!(ones(3, 3); magnitude=NaN)
    @test_throws DimensionMismatch filterSmallPhaseIslands!(_island_field((9, 9), ones(8, 8)), :g)
end

# A uniform grid over the unit box; 17 nodes per axis give spacing exactly 1/16.
function _subres_field(dims, phi)
    axes = ntuple(d -> collect(range(0.0, 1.0; length=dims[d])), length(dims))
    grid = EvolvingDomains.Grid(axes)
    return EvolvingDomains.Field(EvolvingDomains.Discretization(Dict(:g => grid)),
                                 Dict(:g => phi))
end

@testset "sub-resolution feature filter" begin
    # The exact nodal law: an axis-aligned band of t nodes is removed iff t <= 2k.
    for k in (1, 2), t in (2k - 1, 2k, 2k + 1)
        phi = ones(15, 15)
        phi[8:8 + t - 1, :] .= -1.0
        filterSubResolutionFeatures!(phi; radius=k, magnitude=1.0)
        @test (sum(phi .< 0) == 0) == (t <= 2k)
    end

    # Oblique bands follow the lattice threshold (measured anisotropy of the digital ball).
    function _band(half)
        phi = ones(15, 15)
        for I in CartesianIndices(phi)
            abs(I[1] - I[2]) <= half && (phi[I] = -1.0)
        end
        return phi
    end
    phi = _band(0)                                    # ~0.7h wide
    filterSubResolutionFeatures!(phi; radius=1, magnitude=1.0)
    @test all(phi .> 0)
    phi = _band(1)                                    # ~2.1h wide
    filterSubResolutionFeatures!(phi; radius=1, magnitude=1.0)
    @test sum(phi .< 0) == 43
    phi = _band(1)                                    # same band, larger element
    filterSubResolutionFeatures!(phi; radius=2, magnitude=1.0)
    @test all(phi .> 0)

    # An isolated island goes; an axis-aligned block survives with its corners rounded (a
    # corner never contains the ball), which is exactly four nodes at radius 1.
    phi = ones(17, 17)
    phi[2:11, 2:11] .= -1.0
    phi[14, 14] = -1.0
    stats = filterSubResolutionFeatures!(phi; radius=1, magnitude=0.5)
    @test phi[14, 14] == 0.5 && stats.n_flipped_to_negative == 0
    @test stats.n_flipped_to_positive == 5            # island + 4 block corners
    @test sum(phi .< 0) == 10 * 10 - 4
    @test all(phi[2:11, 2:11] .< 0) == false
    @test all(phi[3:10, 3:10] .< 0)                   # interior untouched
    @test phi[2, 2] == 0.5 && phi[11, 11] == 0.5      # corners rounded

    # A thin bridge splits its region only if it is longer than the element: erosion deletes
    # every bridge node, dilation restores the `k` nodes nearest each root, so a neck of
    # nodal length <= 2k survives and one longer than 2k is severed.
    function _neck(length)
        phi = ones(21, 21)
        phi[3:7, 3:7] .= -1.0
        phi[3:7, 8 + length:12 + length] .= -1.0
        phi[5, 8:7 + length] .= -1.0
        return phi
    end
    for (length, severed) in ((1, false), (2, false), (3, true), (5, true))
        phi = _neck(length)
        stats = filterSubResolutionFeatures!(phi; radius=1, magnitude=1.0)
        @test stats.n_components_before == 1
        @test stats.n_components_after == (severed ? 2 : 1)
    end
    phi = _neck(5)
    filterSubResolutionFeatures!(phi; radius=1, magnitude=1.0)
    @test phi[5, 11] > 0                               # the middle of the neck is gone

    phi = ones(17, 17)
    phi[3:8, 3:7] .= -1.0
    phi[3:8, 9:13] .= -1.0
    stats = filterSubResolutionFeatures!(phi; radius=1, magnitude=1.0, mode=:closing)
    @test stats.n_components_before == 2 && stats.n_components_after == 1
    @test stats.n_flipped_to_positive == 0

    # A one-node cavity is filled by closing, and left alone by opening.
    phi = fill(-1.0, 11, 11)
    phi[6, 6] = 1.0
    @test filterSubResolutionFeatures!(copy(phi); radius=1, mode=:opening).n_changed_nodes == 0
    stats = filterSubResolutionFeatures!(phi; radius=1, mode=:closing, magnitude=0.25)
    @test phi[6, 6] == -0.25 && stats.n_flipped_to_positive == 0

    # The positive phase is filtered through the same duality.
    phi = fill(-1.0, 11, 11)
    phi[6, 6] = 1.0
    stats = filterSubResolutionFeatures!(phi; phase=:positive, mode=:opening, radius=1,
                                         magnitude=0.25)
    @test phi[6, 6] == -0.25 && stats.n_flipped_to_negative == 1

    # Invariants: opening is anti-extensive and idempotent, closing extensive and idempotent.
    phi = ones(17, 17)
    phi[4:9, 4:9] .= -1.0
    phi[6, 11:13] .= -1.0
    phibefore = copy(phi)
    filterSubResolutionFeatures!(phi; radius=1, magnitude=1.0)
    @test all(phi .>= phibefore)                      # no node entered the phase
    @test filterSubResolutionFeatures!(phi; radius=1, magnitude=1.0).n_changed_nodes == 0

    # 3D: the same element is a six-neighbour ball, so a lone node is removed and a plane
    # is kept.
    phi3 = ones(5, 5, 5)
    phi3[3, 3, 3] = -1.0
    filterSubResolutionFeatures!(phi3; radius=1, magnitude=1.0)
    @test all(phi3 .> 0)
    phi3 = ones(5, 5, 5)
    phi3[3, :, :] .= -1.0                             # a one-node-thick slab is sub-resolution
    filterSubResolutionFeatures!(phi3; radius=1, magnitude=1.0)
    @test all(phi3 .> 0)
    # At radius 1 the 3D ball is the six-neighbour cross, so the retained part of a 5^3 block
    # is the cross-dilation of its 3^3 core: the centre and the face centres of the cube
    # survive, its corners and face diagonals do not (44 of 125 nodes changed).
    phi3 = ones(7, 7, 7)
    phi3[2:6, 2:6, 2:6] .= -1.0
    stats = filterSubResolutionFeatures!(phi3; radius=1, magnitude=1.0)
    @test stats.n_changed_nodes == 44
    @test phi3[4, 4, 4] < 0 && phi3[2, 4, 4] < 0
    @test phi3[2, 2, 2] > 0 && phi3[2, 2, 4] > 0

    # The Field method writes |φ| = h and leaves untouched nodes bit-for-bit.
    phi = ones(17, 17)
    phi[10, 10] = -1.0
    before = copy(phi)
    field = _subres_field((17, 17), phi)
    stats = filterSubResolutionFeatures!(field, :g; radius=1)
    @test phi[10, 10] == 1 / 16 && stats.n_changed_nodes == 1
    @test all(phi[I] == before[I] for I in CartesianIndices(phi) if I != CartesianIndex(10, 10))

    @test_throws ArgumentError filterSubResolutionFeatures!(ones(5, 5); phase=:both)
    @test_throws ArgumentError filterSubResolutionFeatures!(ones(5, 5); mode=:both)
    @test_throws ArgumentError filterSubResolutionFeatures!(ones(5, 5); radius=0)
    @test_throws ArgumentError filterSubResolutionFeatures!(ones(5, 5); magnitude=0.0)
    @test_throws ArgumentError filterSubResolutionFeatures!(ones(5, 5); magnitude=Inf)
    @test_throws DimensionMismatch filterSubResolutionFeatures!(_subres_field((9, 9), ones(8, 8)), :g)
end
