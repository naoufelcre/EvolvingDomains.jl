using Test
using Main.EvolvingDomains
using Main.EvolvingDomains: refine!, coarsen!, attach!, detach!, immediateCoarser,
    immediateFiner, coarsest, finest, transfer!, linearInterpolator

_axis(n) = collect(range(0.0, 1.0; length=n))

@testset "Multilevel discretization" begin
    @testset "shorthand construction" begin
        d = Discretization(range(0.0, 1.0; length=9), range(0.0, 1.0; length=9))
        @test collect(keys(d.levels)) == [:init]
        @test length(d.levels[:init].axes[1]) == 9
    end

    @testset "shared-domain invariant" begin
        g9, g17 = Grid((_axis(9), _axis(9))), Grid((_axis(17), _axis(17)))
        d = Discretization(Dict(:coarse => g9, :fine => g17))
        @test length(d.levels) == 2          # same box, different resolution
        @test_throws ArgumentError Discretization(Dict(:a => g9, :b => Grid((1.0 .+ _axis(9), _axis(9)))))
        @test_throws ArgumentError Discretization(Dict(:a => g9, :b => Grid((_axis(9), 2.0 .* _axis(9)))))
        @test_throws ArgumentError Discretization(Dict(:a => g9, :b => Grid((_axis(9), _axis(9), _axis(9)))))
    end

    @testset "refine!/coarsen! preserve the domain and nest" begin
        d = Discretization(Dict(:init => Grid((_axis(9), _axis(9)))))
        mid = refine!(d, :init; factor=2)
        @test mid == :init_r2
        @test length(d.levels[mid].axes[1]) == 17
        @test first(d.levels[mid].axes[1]) == 0.0 && last(d.levels[mid].axes[1]) == 1.0
        fine = refine!(d, mid, :fine)                       # explicit name
        @test fine == :fine && length(d.levels[:fine].axes[2]) == 33
        clash = refine!(d, :init)                           # auto name on collision
        @test clash == :init_r2_2
        @test_throws ArgumentError refine!(d, :init, :fine)  # explicit name already taken
        @test_throws ArgumentError refine!(d, :init; factor=1)
        @test_throws ArgumentError refine!(d, :missing)

        back = coarsen!(d, :fine; factor=2)
        @test back == :fine_c2
        @test d.levels[back].axes[1] == d.levels[mid].axes[1]                     # exact subsample
        bad = Discretization(Dict(:only => Grid((_axis(10), _axis(9)))))
        @test_throws ArgumentError coarsen!(bad, :only; factor=2)                 # 9 intervals
        @test_throws ArgumentError coarsen!(bad, :only; factor=4)                 # 8 not divisible by 4 for axis 2
    end

    @testset "attach!/detach! validate and retain" begin
        d = Discretization(Dict(:init => Grid((_axis(9), _axis(9)))))
        f = Field(d, Dict(:init => zeros(9, 9)))
        refine!(d, :init, :fine)
        fine = zeros(17, 17)
        @test attach!(f, :fine, fine) === f
        @test f.data[:fine] === fine                        # retained, not copied
        @test_throws ArgumentError attach!(f, :nope, zeros(5, 5))
        @test_throws ArgumentError attach!(f, :fine, zeros(Float32, 17, 17))
        @test_throws ArgumentError attach!(f, :fine, zeros(17))
        @test_throws DimensionMismatch attach!(f, :fine, zeros(9, 17))
        @test detach!(f, :fine) === f
        @test !haskey(f.data, :fine) && haskey(d.levels, :fine)
        @test_throws ArgumentError detach!(f, :fine)
    end

    @testset "navigation by resolution" begin
        d = Discretization(Dict(:init => Grid((_axis(9), _axis(9)))))
        refine!(d, :init, :mid)
        refine!(d, :mid, :fine)
        @test coarsest(d) == :init
        @test finest(d) == :fine
        @test immediateFiner(d, :init) == :mid
        @test immediateFiner(d, :mid) == :fine
        @test immediateFiner(d, :fine) === nothing
        @test immediateCoarser(d, :fine) == :mid
        @test immediateCoarser(d, :init) === nothing
        @test_throws ArgumentError immediateFiner(d, :missing)
        # equal resolution is ordered by name, so navigation stays total: coarsening back
        # onto an existing resolution is a normal workflow, not an error
        tied = Discretization(Dict(:a => Grid((_axis(9), _axis(9))), :b => Grid((_axis(9), _axis(9)))))
        @test coarsest(tied) == :a && finest(tied) == :b
        @test immediateFiner(tied, :a) == :b
        @test immediateFiner(tied, :b) === nothing
    end

    @testset "transfer! reproduces a linear field exactly" begin
        d = Discretization(Dict(:init => Grid((_axis(9), _axis(9)))))
        refine!(d, :init, :fine)
        exact(x, y) = 2.0 + 3.0x - 4.0y
        f = Field(d, Dict(:init => [exact(x, y) for x in _axis(9), y in _axis(9)]))
        attach!(f, :fine, zeros(17, 17))

        transfer!(f, :init, :fine)                          # prolongation
        @test maximum(abs.([f.data[:fine][i, j] - exact(_axis(17)[i], _axis(17)[j])
                            for i in 1:17, j in 1:17])) < 1e-14

        fill!(f.data[:init], 0.0)
        transfer!(f, :fine, :init)                          # restriction by sampling
        @test maximum(abs.([f.data[:init][i, j] - exact(_axis(9)[i], _axis(9)[j])
                            for i in 1:9, j in 1:9])) < 1e-14

        @test_throws ArgumentError transfer!(f, :init, :init)
        detach!(f, :fine)
        @test_throws ArgumentError transfer!(f, :init, :fine)   # target not attached

        # non-nested levels over the same box: interpolation would work, transfer refuses
        d2 = Discretization(Dict(:a => Grid((_axis(9), _axis(9))), :b => Grid((_axis(12), _axis(12)))))
        f2 = Field(d2, Dict(:a => zeros(9, 9)))
        attach!(f2, :b, zeros(12, 12))
        @test_throws ArgumentError transfer!(f2, :a, :b)
    end

    @testset "the interpolator is pluggable" begin
        d = Discretization(Dict(:coarse => Grid((_axis(9), _axis(9)))))
        refine!(d, :coarse, :fine)
        f = Field(d, Dict(:coarse => fill(7.0, 9, 9)))
        attach!(f, :fine, zeros(17, 17))

        transfer!(f, :coarse, :fine; itp=(x, values, grid) -> 42.0)     # a plain function
        @test all(==(42.0), f.data[:fine])

        transfer!(f, :coarse, :fine; itp=(x, values, grid) -> x[1] + 2x[2])
        @test all(isapprox(f.data[:fine][i, j], _axis(17)[i] + 2 * _axis(17)[j]; atol=1e-15)
                  for i in 1:17, j in 1:17)

        # a closure over per-grid precomputation has the same shape
        displacement = 0.25
        transfer!(f, :coarse, :fine;
                  itp=(x, values, grid) -> linearInterpolator((x[1] - displacement, x[2]), values, grid))
        @test f.data[:fine][9, 9] ≈ 7.0

        # the default reproduces the source on the nodes it shares with the target
        transfer!(f, :coarse, :fine)
        @test all(isapprox(f.data[:fine][i, j], 7.0; atol=1e-15)
                  for i in 1:2:17, j in 1:2:17)
    end

    @testset "permissive container, strict use" begin
        # The split this package commits to: mutating `d.levels` is not intercepted;
        # *using* a level is validated, once, through `_levelgrid`/`_levelvalues`.
        d = Discretization(range(0.0, 1.0; length=9), range(0.0, 1.0; length=9))
        f = Field(d, Dict(:init => zeros(9, 9)))

        # 1. a foreign grid can be grafted in by hand, silently
        foreign = Grid((1.0 .+ _axis(9), _axis(9)))       # box shifted off the domain
        d.levels[:foreign] = foreign
        @test haskey(d.levels, :foreign)

        # 2. on its own that level is a consistent single-level problem, so it is usable:
        #    the checks apply to the level being used, not to its siblings.
        @test attach!(f, :foreign, zeros(9, 9)) === f
        f.data[:foreign] .= [x - 1.5 for x in (1.0 .+ _axis(9)), y in _axis(9)]
        @test length(cut(f, :foreign).bg_states) == 64

        # 3. as soon as it must agree with another level, it is refused
        @test_throws ArgumentError transfer!(f, :init, :foreign)
        @test_throws ArgumentError transfer!(f, :foreign, :init)

        # 4. navigation only orders by resolution, so a foreign box does not panic it
        @test immediateFiner(d, :foreign) == :init         # equal resolution, ordered by name
        @test coarsest(d) == :foreign && finest(d) == :init

        # 5. shape is validated against the level's own grid, at use not at mutation
        d.levels[:bad] = Grid((_axis(17), _axis(17)))
        @test_throws DimensionMismatch attach!(f, :bad, zeros(9, 9))   # shape -> DimensionMismatch (D4 convention)
        @test_throws ArgumentError detach!(f, :bad)        # nothing attached there yet
    end
end

@testset "nestedness tolerates representation" begin
    # The same mesh written two ways can differ in the last bit; it is still one mesh.
    d = Discretization(range(0.0, 1.0; length=17), range(0.0, 1.0; length=17))
    refined = refine!(d, :init; factor=3)                     # 49 nodes, by refine!
    d.levels[:written] = Grid((collect(range(0.0, 1.0; length=49)),
                               collect(range(0.0, 1.0; length=49))))
    @test maximum(abs.(d.levels[refined].axes[1] .- d.levels[:written].axes[1])) > 0

    exact(x, y) = 2x - 3y
    ax = d.levels[:init].axes
    f = Field(d, Dict(:init => [exact(x, y) for x in ax[1], y in ax[2]]))
    attach!(f, refined, zeros(49, 49))
    attach!(f, :written, zeros(49, 49))
    @test transfer!(f, :init, :written) === f
    axw = d.levels[:written].axes
    @test maximum(abs.(f.data[:written] .- [exact(x, y) for x in axw[1], y in axw[2]])) < 1e-14
end

@testset "whole-number meshes refine refuses, coarsen does not" begin
    # A grid of integers is a legal container (geometry on it is exact), but a midpoint is
    # not an integer: refinement refuses where the fraction is needed, rather than letting
    # an internal conversion error surface.
    d = Discretization(Dict(:init => Grid(([0, 1, 2], [0, 1, 2]))))
    @test eltype(d.levels[:init].axes[1]) === Int
    @test_throws ArgumentError refine!(d, :init)
    coarse = coarsen!(d, :init)                       # subsampling needs no fractions
    @test d.levels[coarse].axes[1] == [0, 2]
end
