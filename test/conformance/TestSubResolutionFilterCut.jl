using Test
using Main.EvolvingDomains
using Main.EvolvingDomains: cutCellIds, CUT_IN, CUT_OUT, CUT_CUT
using Main.EvolvingDomains.Geometric.TopologyFilter: _squaredDistanceTransform

# A signed level set with |φ| = O(h) at the interface, which is what the filter's write-back
# assumes: the node values are node distances to the other phase, offset by half a spacing so
# that every zero crossing sits at a cell midpoint and no node is exactly zero.
function _signedphi(mask, h)
    return [mask[I] ? -h * (sqrt(_squaredDistanceTransform(map(!, mask))[I]) + 0.5) :
                       h * (sqrt(_squaredDistanceTransform(mask)[I]) + 0.5)
            for I in CartesianIndices(mask)]
end

function _subarea(snapshot, subcell)
    ids = Int(snapshot.point_pointers[subcell]):(Int(snapshot.point_pointers[subcell + 1]) - 1)
    a, b, c = Int.(snapshot.point_ids[ids])
    p, q, r = snapshot.physical_points[[a, b, c]]
    return abs((q[1] - p[1]) * (r[2] - p[2]) - (q[2] - p[2]) * (r[1] - p[1])) / 2
end

# Area of the negative phase of the reconstructed geometry, from the subcell measures.
function _negative_area(snapshot)
    return sum(_subarea(snapshot, i) for i in eachindex(snapshot.subcell_states)
               if snapshot.subcell_states[i] == CUT_IN; init=0.0)
end

# Face-adjacency components of the cells the negative phase occupies at all (CUT_IN or CUT_CUT).
function _cell_components(snapshot, ncells)
    occupied = BitVector(getproperty(snapshot, :bg_states) .!= CUT_OUT)
    visited = falses(ncells)
    stack = Int[]
    components = 0
    for start in 1:ncells
        (visited[start] || !occupied[start]) && continue
        components += 1
        push!(stack, start)
        visited[start] = true
        while !isempty(stack)
            cell = pop!(stack)
            for side in 1:4
                axis = (side - 1) ÷ 2 + 1
                step = side % 2 == 1 ? -1 : 1
                dims = (16, 16)
                index = (cell - 1) ÷ prod(dims[1:axis - 1]) % dims[axis] + 1
                neighbour = index + step
                1 <= neighbour <= dims[axis] || continue
                next = cell + step * prod(dims[1:axis - 1])
                (visited[next] || !occupied[next]) && continue
                visited[next] = true
                push!(stack, next)
            end
        end
    end
    return components
end

@testset "sub-resolution filter, reconstructed cut" begin
    n = 17
    h = 1 / (n - 1)
    grid = Grid((collect(range(0.0, 1.0; length=n)), collect(range(0.0, 1.0; length=n))))

    function build(mask)
        phi = reshape(_signedphi(mask, h), size(mask))
        field = Field(Discretization(Dict(:g => grid)), Dict(:g => phi))
        return field, phi
    end

    # A thin bridge is severed in the *reconstructed* geometry: before filtering the negative
    # cells form one component, after it two, and the cut loses cells.
    mask = falses(n, n)
    mask[3:7, 3:7] .= true
    mask[3:7, 12:16] .= true
    mask[5, 8:11] .= true
    field, phi = build(mask)
    before = cut(field, :g)
    stats = filterSubResolutionFeatures!(field, :g; mode=:opening, radius=1)
    after = cut(field, :g)
    @test _cell_components(before, 16 * 16) == 1
    @test _cell_components(after, 16 * 16) == 2
    @test stats.n_components_before == 1 && stats.n_components_after == 2
    @test length(cutCellIds(after)) <= length(cutCellIds(before))
    # The neck's own cells are entirely outside the phase now: no sliver of it survives.
    neck_cells = [i + 16 * 3 for i in 8:11]        # cells (row 4, cols 8:11)
    @test all(after.bg_states[neck_cells] .== CUT_OUT)
    # Every cell the cut does keep has positive measure; nothing degenerates.
    @test minimum(_subarea(after, i) for i in eachindex(after.subcell_states)) > 0.0
    # Area balance: the removed neck is 4 nodes long and about one cell wide.
    @test _negative_area(before) - _negative_area(after) > 0.9 * h^2

    # An isolated node leaves no trace in the cut geometry.
    mask = falses(n, n)
    mask[3:7, 3:7] .= true
    mask[11, 11] = true
    field, phi = build(mask)
    before = cut(field, :g)
    filterSubResolutionFeatures!(field, :g; mode=:opening, radius=1)
    after = cut(field, :g)
    @test phi[11, 11] > 0
    hole_cells = [i + 16 * (j - 1) for i in 10:11 for j in 9:10]
    @test all(after.bg_states[hole_cells] .== CUT_OUT)
    @test length(cutCellIds(after)) < length(cutCellIds(before))
    @test _negative_area(before) > _negative_area(after)

    # A one-node cavity in a bulk that reaches the domain boundary is filled: the field is
    # single-signed afterwards, so the snapshot carries no interface at all.
    mask = trues(n, n)
    mask[9, 9] = false
    field, phi = build(mask)
    before = cut(field, :g)
    @test length(cutCellIds(before)) == 4          # the four cells around the hole
    filterSubResolutionFeatures!(field, :g; mode=:closing, radius=1)
    after = cut(field, :g)
    @test isempty(after.boundary.facet_to_bgcell)
    @test all(after.bg_states .== CUT_IN)
    @test _negative_area(after) ≈ 16 * 16 * h^2
end
