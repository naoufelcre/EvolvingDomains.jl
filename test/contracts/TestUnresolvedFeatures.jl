using Test
using Main.EvolvingDomains
using Main.EvolvingDomains: Grid, Discretization, Field, filterUnresolvedFeatures!,
    CUT_IN, CUT_OUT, CUT_CUT

# Contract of `filterUnresolvedFeatures!` (Luo, Hu & Adams 2016, §3; project adaptations in
# src/Geometric/UnresolvedFeatures.jl). No native CutCells library is used: the `Cut` argument
# is hand-built. Behaviour is observed through the nodal array only.
#
# eps = 0.75h on every grid below. Values are exact analytic signed distances, so a resolved
# node is never written and exact equality (not a tolerance) is the correct expectation.

_axis(n) = collect(range(0.0, 1.0; length=n))
_hstep(n) = 1.0 / (n - 1)
_eps(n) = 0.75 * _hstep(n)

function _field2(n, f)
    ax = _axis(n)
    Field(Discretization(ax, ax), Dict(:init => [f(x, y) for x in ax, y in ax]))
end

function _field3(n, f)
    ax = _axis(n)
    Field(Discretization(Dict(:init => Grid((copy(ax), copy(ax), copy(ax))))),
          Dict(:init => [f(x, y, z) for x in ax, y in ax, z in ax]))
end

# Run the filter and return the arrays before and after.
function _run!(field, level=:init; kwargs...)
    before = copy(field.data[level])
    returned = filterUnresolvedFeatures!(field, level; kwargs...)
    @test returned === field                 # documented return value: the field itself
    return before, field.data[level]
end

# Grid edges whose endpoints have strictly opposite signs (a zero crossing of phi = 0).
function _crossings(phi::AbstractArray{T,N}) where {T,N}
    n = 0
    for I in CartesianIndices(phi), d in 1:N
        J = CartesianIndex(ntuple(k -> k == d ? I[k] + 1 : I[k], N))
        checkbounds(Bool, phi, J) || continue
        a, b = phi[I], phi[J]
        ((a > 0 && b < 0) || (a < 0 && b > 0)) && (n += 1)
    end
    return n
end

# Background-cell states computed independently of the implementation: a cell is CUT when its
# corner nodes are not all on one side of zero. Linear index: first axis fastest, as documented.

@testset "unresolved-feature filter contract" begin

    @testset "resolved geometry is not written" begin
        # A resolved cut cell has both auxiliary crossings (phi = ±eps) inside its
        # 3-neighbourhood, so it is classified S2 and never enters either correction path.
        # No node is written, so the arrays must agree bit for bit (exact equality, no
        # tolerance): every field below is an exact distance function with |∇phi| = 1 and the
        # smallest resolved structure is several cells wide.
        cases = ["circle r=0.3, 33²" => _field2(33, (x, y) -> hypot(x - 0.5, y - 0.5) - 0.3),
                 "plane x=0.5, 33²" => _field2(33, (x, y) -> x - 0.5),
                 "square corner, 33²" => _field2(33, (x, y) -> max(abs(x - 0.5), abs(y - 0.5)) - 0.25),
                 "sphere r=0.3, 17³" => _field3(17, (x, y, z) -> hypot(x - 0.5, y - 0.5, z - 0.5) - 0.3),
                 "plane x=0.5, 17³" => _field3(17, (x, y, z) -> x - 0.5)]
        for (name, field) in cases
            before, after = _run!(field)
            @info "resolved geometry" name maxchange=maximum(abs.(after .- before)) changed=count(!=(0.0), after .- before)
            @test after == before
            @test maximum(abs.(after .- before)) == 0.0     # the reported number, not "small"
        end
    end

    @testset "isolated islands of either sign are removed" begin
        # 0 < max|phi| = a < eps = 0.0234: the required auxiliary level set does not exist in
        # the search region, so the project island rule (not Eq. 15) replaces the island nodes.
        h33, a = _hstep(33), 0.02
        @test a < _eps(33)
        pos = _field2(33, (x, y) -> a - hypot(x - 0.5, y - 0.5))
        neg = _field2(33, (x, y) -> hypot(x - 0.5, y - 0.5) - a)
        bpos, apos = _run!(pos)
        bneg, aneg = _run!(neg)
        @info "2D island" eps=_eps(33) maxpositive=maximum(bpos) maxnegative=-minimum(bneg) changed=(count(!=(0.0), apos .- bpos), count(!=(0.0), aneg .- bneg))
        @test count(>(0), bpos) == 1 && count(<(0), bneg) == 1        # one-node islands
        @test count(>(0), apos) == 0                                   # completely removed
        @test count(<(0), aneg) == 0
        @test _crossings(apos) == 0 && _crossings(aneg) == 0           # no residual crossings
        @test all(isfinite, apos) && all(isfinite, aneg)
        # Two phases behave symmetrically: negating the input negates the output.
        @test aneg == -apos

        # 3D mirror of the same check (eps = 0.0469 on 17³).
        a3 = 0.01
        @test a3 < _eps(17)
        pos3 = _field3(17, (x, y, z) -> a3 - hypot(x - 0.5, y - 0.5, z - 0.5))
        neg3 = _field3(17, (x, y, z) -> hypot(x - 0.5, y - 0.5, z - 0.5) - a3)
        bpos3, apos3 = _run!(pos3)
        bneg3, aneg3 = _run!(neg3)
        @info "3D island" eps=_eps(17) maxpositive=maximum(bpos3) changed=(count(!=(0.0), apos3 .- bpos3), count(!=(0.0), aneg3 .- bneg3))
        @test count(>(0), bpos3) == 1 && count(<(0), bneg3) == 1
        @test count(>(0), apos3) == 0 && count(<(0), aneg3) == 0
        @test _crossings(apos3) == 0 && _crossings(aneg3) == 0
        @test all(isfinite, apos3) && all(isfinite, aneg3)
        @test aneg3 == -apos3

        # An island spanning two nodes is one connected component: 65², centre placed between
        # the nodes at x = 0.5 and x = 0.5 + h, a = 0.011 < eps = 0.0117.
        h65, a65 = _hstep(65), 0.011
        @test a65 < _eps(65)
        pair = _field2(65, (x, y) -> a65 - hypot(x - (0.5 + h65 / 2), y - 0.5))
        bpair, apair = _run!(pair)
        @info "two-node island" eps=_eps(65) positives=count(>(0), bpair) changed=count(!=(0.0), apair .- bpair)
        @test count(>(0), bpair) == 2
        @test count(>(0), apair) == 0
        @test count(!=(0.0), apair .- bpair) == 2                      # both island nodes rewritten
        @test _crossings(apair) == 0
        @test all(isfinite, apair)

        # Unrelated geometry is untouched: a drop inside a resolved disk (33²). phi = max(...)
        # keeps the disk's zero set and adds one positive node inside it.
        drop = _field2(33, (x, y) -> max(hypot(x - 0.5, y - 0.5) - 0.3, a - hypot(x - 0.5, y - 0.5)))
        bdrop, adrop = _run!(drop)
        changed = [(i, j) for i in 1:33, j in 1:33 if adrop[i, j] != bdrop[i, j]]
        @info "island inside a resolved disk" changed=length(changed) maxphi=maximum(bdrop)
        @test length(changed) == 1                                     # only the island node
        @test maximum(adrop) == maximum(bdrop)                         # disk geometry intact
        mask = trues(33, 33)
        for (i, j) in changed
            mask[i, j] = false
        end
        @test adrop[mask] == bdrop[mask]                               # byte-identical elsewhere
        @test count(>(0), adrop) == count(>(0), bdrop) - 1             # the drop is gone
    end

    @testset "thin connected features" begin
        h, w = _hstep(33), 0.5 * _hstep(33)
        capsule(cx, cy, L, halfwidth) =
            (x, y) -> halfwidth - hypot(max(abs(x - cx) - L, 0.0), y - cy)

        # A 1-node-wide filament made of a positive phase: |phi| <= w = 0.5h < eps everywhere,
        # so no +eps crossing exists and the cut cells fail the Eq. (9) test. The published
        # reconstruction has no source, and the project island rule deletes the filament.
        filament = _field2(33, capsule(0.5, 0.5, 0.2, w))
        bfil, afil = _run!(filament)
        @info "1-node filament" eps=_eps(33) halfwidth=w changed=count(!=(0.0), afil .- bfil) crossings=(_crossings(bfil), _crossings(afil))
        @test count(>(0), bfil) > 10
        @test count(>(0), afil) == 0
        @test _crossings(afil) == 0
        @test all(isfinite, afil)

        # The published locality criterion tolerates a structure that is two nodes wide: with
        # half-width h the filament reaches |phi| = h > eps, so its cut cells still have a
        # +eps crossing in their 3-neighbourhood and are classified resolved (S2). Nothing is
        # written. This is Eq. (10), not a bug: the threshold is eps = 0.75h, not h.
        wide = _field2(33, capsule(0.5, 0.5, 0.2, h))
        bwide, awide = _run!(wide)
        @info "2-node filament kept" threshold=_eps(33) maxpositive=maximum(bwide)
        @test awide == bwide
        @test maximum(bwide) > _eps(33)

        # A 1-node-wide protrusion of the negative phase attached to a resolved disk. Its cut
        # cells miss the -eps crossing, so the published Eq. (15) re-distances them from the
        # disk's +eps surface: phi_s = eps - d, positive. The finger is erased and the values
        # away from the disk become eps - d > 0, not a bolt-on structure.
        finger = _field2(33, (x, y) -> min(hypot(x - 0.5, y - 0.5) - 0.2,
                                          -(w - hypot(max(abs(x - 0.7) - 0.1, 0.0), y - 0.5))))
        bfing, afing = _run!(finger)
        ax = _axis(33)
        touched = [(i, j) for i in 1:33, j in 1:33 if afing[i, j] != bfing[i, j]]
        after_vals = [afing[i, j] for (i, j) in touched]
        before_vals = [bfing[i, j] for (i, j) in touched]
        @info "attached protrusion" changed=length(touched) maxchange=maximum(abs.(afing .- bfing)) xrange=(minimum(ax[i] for (i, j) in touched), maximum(ax[i] for (i, j) in touched))
        @test !isempty(touched)
        @test all(ax[i] >= 0.7 for (i, j) in touched)                  # disk side untouched
        @test all(abs(v) < _eps(33) for v in before_vals)              # sub-eps structure before
        @test count(<(0), before_vals) > 0                             # the finger is negative
        @test all(>(0), after_vals)                                    # erased: the phase fills in
        @test maximum(afing) == maximum(bfing)                         # no new extremum

        # Thin bridge between two resolved disks: a 1-node-wide negative neck (|phi| = 0.5h
        # < eps) is broken for a 3-cell and a 4-cell gap. The correction is local, so the
        # discrete negative run through the neck is interrupted but the two junction nodes
        # (shared with resolved cells) stay at -0.5h: sub-eps stubs survive. Both facts are
        # asserted below as measured; the disks themselves stay byte-identical.
        let r = 0.15
            neck(gap, halfwidth) = begin
                c1, c2 = 0.5 - (r + gap / 2), 0.5 + (r + gap / 2)
                (x, y) -> min(hypot(x - c1, y - 0.5) - r, hypot(x - c2, y - 0.5) - r,
                              -(halfwidth - hypot(max(abs(x - 0.5) - (r + gap / 2), 0.0), y - 0.5)))
            end
            for gapcells in (3, 4)
                bridge = _field2(33, neck(gapcells * h, w))
                bbr, abr = _run!(bridge)
                touched = [(i, j) for i in 1:33, j in 1:33 if abr[i, j] != bbr[i, j]]
                # midline of the neck: row index of y = 0.5 is 17; the neck spans i = 16..18
                midline_before = [bbr[i, 17] for i in 16:18]
                midline_after = [abr[i, 17] for i in 16:18]
                @info "thin bridge" gapcells=gapcells changed=length(touched) midline_before=midline_before midline_after=midline_after
                @test all(<(0), midline_before)                        # the neck was continuous
                # Measured: the middle node of the neck is re-distanced into the positive phase,
                # so the discrete negative run — and with it the bridge — is broken. The two
                # nodes shared with the resolved disks are left at -0.5h < -eps: the correction
                # is local, so a sub-resolution stub survives at each junction. Asserted as
                # observed, not assumed.
                @test all(abs(v) < _eps(33) for v in midline_before)
                @test any(>(0), midline_after)                         # rewritten to the other phase
                @test !all(<(0), midline_after)                        # no continuous bridge left
                @test count(<(0), midline_after) <= 2                  # ... only junction stubs remain
                @test 0 < length(touched) <= 9                         # middle column (measured 3)
                @test all(hypot(ax[i] - 0.5, ax[j] - 0.5) < 0.12 for (i, j) in touched)
                @test maximum(abr) == maximum(bbr)                     # disk geometry intact
                @test all(isfinite, abr)
            end
            # A neck of one or two cells between the same disks is NOT flagged: its cut cells
            # still see both auxiliary crossings inside the 3-neighbourhood. This is the
            # published locality criterion (Eq. 10) with eps = 0.75h, not a bug and not hidden
            # here: the non-detection is asserted as the observed behaviour.
            for gapcells in (1, 2)
                bridge = _field2(33, neck(gapcells * h, w))
                bbr, abr = _run!(bridge)
                @info "short neck not flagged" gapcells=gapcells changed=count(!=(0.0), abr .- bbr) crossings=_crossings(bbr)
                @test _crossings(bbr) > 0            # there is a thin bridge to be missed
                @test abr == bbr
            end
        end

        # A small cavity: an isolated negative region of depth 0.02 < eps inside a resolved
        # positive phase. Same island rule as above; the surrounding phase fills it in.
        cavity = _field2(33, (x, y) -> hypot(x - 0.5, y - 0.5) - 0.02)
        bcav, acav = _run!(cavity)
        @info "small cavity" depth=-minimum(bcav) eps=_eps(33) changed=count(!=(0.0), acav .- bcav)
        @test minimum(bcav) > -_eps(33)
        @test count(<(0), acav) == 0
        @test _crossings(acav) == 0
        @test all(isfinite, acav)
    end

    @testset "idempotence, finiteness, locality" begin
        # Filtering a resolved field is a no-op, and filtering a filtered field is too: the
        # second pass sees only resolved cells.
        circle = _field2(33, (x, y) -> hypot(x - 0.5, y - 0.5) - 0.3)
        _, once = _run!(circle)
        _, twice = _run!(circle)
        @test twice == once
        island = _field2(33, (x, y) -> 0.02 - hypot(x - 0.5, y - 0.5))
        _, once = _run!(island)
        _, twice = _run!(island)
        @test twice == once
        @test all(isfinite, twice)

        # Local change only: a resolved field with one island changes exactly the island node
        # and no other node's sign structure.
        mixed = _field2(33, (x, y) -> max(hypot(x - 0.5, y - 0.5) - 0.3, 0.02 - hypot(x - 0.5, y - 0.5)))
        bmix, amix = _run!(mixed)
        changed = [(i, j) for i in 1:33, j in 1:33 if amix[i, j] != bmix[i, j]]
        @test length(changed) == 1
        keep = [CartesianIndex(i, j) for i in 1:33, j in 1:33 if !((i, j) in changed)]
        @test sign.(amix[keep]) == sign.(bmix[keep])                    # sign structure preserved
        @test amix[keep] == bmix[keep]                                 # and the values themselves
    end

    @testset "error paths" begin
        ax9 = _axis(9)
        uniform = Grid((copy(ax9), copy(ax9)))
        d9 = Discretization(Dict(:init => uniform))

        # Non-uniform spacing: the operator converts |phi| to distance, which needs h.
        nonuniform = Grid(([0.0, 0.1, 0.3, 0.6, 1.0], copy(ax9)))
        dnu = Discretization(Dict(:init => nonuniform))
        @test_throws ArgumentError filterUnresolvedFeatures!(Field(dnu, Dict(:init => zeros(5, 9))), :init)

        # Float32 values.
        @test_throws ArgumentError filterUnresolvedFeatures!(Field(d9, Dict(:init => zeros(Float32, 9, 9))), :init)

        # Non-finite values: NaN and Inf.
        @test_throws ArgumentError filterUnresolvedFeatures!(Field(d9, Dict(:init => fill(NaN, 9, 9))), :init)
        infinite = fill(1.0, 9, 9)
        infinite[3, 3] = Inf
        @test_throws ArgumentError filterUnresolvedFeatures!(Field(d9, Dict(:init => infinite)), :init)

        # Unknown level, and a level present in the discretization with no attached values.
        attached = Field(d9, Dict(:init => zeros(9, 9)))
        @test_throws ArgumentError filterUnresolvedFeatures!(attached, :missing)
        d9.levels[:extra] = Grid((copy(ax9), copy(ax9)))
        @test_throws ArgumentError filterUnresolvedFeatures!(attached, :extra)
        # Neither level was attached, so nothing was mutated anywhere.
        @test all(iszero, attached.data[:init])
    end

    @testset "no global reinitialisation" begin
        # Far away from the interface the field is deliberately NOT a unit-slope distance
        # function (values scaled by 1.5 beyond radius 0.4). A reinitialisation would rewrite
        # those nodes; the filter must leave every node it does not resolve untouched.
        ax = _axis(33)
        far(x, y) = begin
            d = hypot(x - 0.5, y - 0.5) - 0.3
            d > 0.1 ? 1.5d : d
        end
        field = Field(Discretization(ax, ax),
                      Dict(:init => [max(far(x, y), 0.02 - hypot(x - 0.5, y - 0.5)) for x in ax, y in ax]))
        before, after = _run!(field)
        h = _hstep(33)
        r = [hypot(x - 0.5, y - 0.5) for x in ax, y in ax]
        @test all(after[r .> 0.4] .== before[r .> 0.4])
        @test maximum(abs.(after[r .> 0.4] .- before[r .> 0.4])) == 0.0
        # The island inside the disk is gone (its four crossings), the rest is untouched.
        @test count(>(0), after[r .< 0.25]) == 0
        @test _crossings(after) == _crossings(before) - 4
        ring = (r .> 0.3 - 3h) .&& (r .< 0.3 + 3h)
        @test after[ring] == before[ring]
        @test all(isfinite, after)
    end
end
