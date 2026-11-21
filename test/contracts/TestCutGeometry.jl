using Test
using LinearAlgebra: norm
using Main.EvolvingDomains: Cut, Grid, Discretization, Field, reinitialize!

function _cut_levelset(grid; level_id=:init, data=zeros(length(grid.axes[1]), length(grid.axes[2])))
    Field(Discretization(Dict(level_id => grid)), Dict(level_id => data))
end

function _cut_for_cells(grid, level_id, states)
    Cut(grid, level_id; bg_states=states)
end

@testset "Cut-based geometry operations" begin
    @testset "identity, level, shape, and finite-value validation" begin
        g = Grid((collect(range(0.0, 1.0; length=5)), collect(range(0.0, 1.0; length=4))))
        f = _cut_levelset(g)
        cut = _cut_for_cells(g, :init, zeros(Int8, 12))
        @test_throws ArgumentError reinitialize!(f, _cut_for_cells(
            Grid((copy(g.axes[1]), copy(g.axes[2]))), :init, zeros(Int8, 12)))

        missing_values_field = Field(Discretization(Dict(:init => g, :fine => g)),
            Dict(:init => zeros(5, 4)))
        @test_throws ArgumentError reinitialize!(missing_values_field,
            _cut_for_cells(g, :fine, zeros(Int8, 12)))

        wrong_shape = _cut_levelset(g; data=zeros(4, 5))
        @test_throws DimensionMismatch reinitialize!(wrong_shape, cut)
        uneven = Grid(([0.0, 0.2, 0.55, 0.8, 1.0], copy(g.axes[2])))
        @test_throws ArgumentError reinitialize!(_cut_levelset(uneven),
            _cut_for_cells(uneven, :init, zeros(Int8, 12)))
        nonfinite = zeros(5, 4)
        nonfinite[2, 2] = NaN
        @test_throws ArgumentError reinitialize!(_cut_levelset(g; data=nonfinite), cut)
    end

    @testset "reinitialization preserves grid-edge zero crossings" begin
        xs, ys = collect(range(-1.0, 1.0; length=9)), collect(range(-1.0, 1.0; length=7))
        g = Grid((xs, ys))
        phi = [((x - 0.13) * (1 + 0.12sin(3y))) for x in xs, y in ys]
        f = _cut_levelset(g; data=phi)
        cell_states = fill(Int8(1), (length(xs) - 1) * (length(ys) - 1))
        cut_i = findfirst(i -> xs[i] <= 0.13 <= xs[i + 1], 1:length(xs)-1)
        for j in 1:length(ys)-1
            cell_states[cut_i + (j - 1) * (length(xs) - 1)] = 0
        end
        cut = _cut_for_cells(g, :init, cell_states)
        before = copy(phi)
        reinitialize!(f, cut)
        after = f.data[:init]
        for j in axes(phi, 2)
            left = findfirst(i -> before[i, j] <= 0 <= before[i + 1, j], 1:length(xs)-1)
            @test isapprox((0 - before[left, j]) / (before[left + 1, j] - before[left, j]),
                           (0 - after[left, j]) / (after[left + 1, j] - after[left, j]); atol=1e-12)
        end
    end

end
