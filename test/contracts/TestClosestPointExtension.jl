using Test
using Main.EvolvingDomains
using Main.EvolvingDomains.Geometric.ClosestPointExtension: extendClosestPoint!

# A planar front has |∇φ| = 1 exactly, so the closest-point map is exact there: the
# extension of a bulk constant must be that constant everywhere in the void.
function _planar_case(dims)
    axes = ntuple(d -> collect(range(0.0, 1.0; length=dims[d])), length(dims))
    grid = EvolvingDomains.Grid(axes)
    phi = Array{Float64}(undef, dims)
    D = length(dims)
    for index in CartesianIndices(phi)
        phi[index] = grid.axes[D][index[D]] - 0.5
    end
    field = EvolvingDomains.Field(EvolvingDomains.Discretization(Dict(:g => grid)),
                                 Dict(:g => phi))
    return grid, phi, field
end

@testset "closest-point extension" begin
    grid, phi, field = _planar_case((9, 9))
    values = [phi[index] <= 0 ? 7.0 : 0.0 for index in CartesianIndices(phi)]
    before = copy(values)
    extendClosestPoint!(values, field, :g)

    @test all(isapprox.(values, 7.0; atol=1e-12))
    @test all(values[index] == before[index] for index in CartesianIndices(phi) if phi[index] <= 0)

    # Bulk values are per node and must come back untouched, not merely rescattered.
    values .= [phi[index] <= 0 ? 1.0 + index[1] + 10index[2] : -1.0 for index in CartesianIndices(phi)]
    before = copy(values)
    extendClosestPoint!(values, field, :g)
    @test all(values[index] == before[index]
              for index in CartesianIndices(phi) if phi[index] <= 0)

    # Void values must be finite: interpolation weighs every corner of the target cell.
    bad = copy(before)
    bad[1, 1] = NaN
    @test_throws ArgumentError extendClosestPoint!(bad, field, :g)
    @test_throws DimensionMismatch extendClosestPoint!(zeros(8, 8), field, :g)   # wrong size

    _, phi3, field3 = _planar_case((5, 5, 5))
    values3 = [phi3[index] <= 0 ? 3.0 : 0.0 for index in CartesianIndices(phi3)]
    extendClosestPoint!(values3, field3, :g)
    @test all(isapprox.(values3, 3.0; atol=1e-12))
end
