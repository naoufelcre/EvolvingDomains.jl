using Test
using Main.EvolvingDomains

@testset "Grid, Discretization, and Field data model" begin
    Grid, Discretization, Field = EvolvingDomains.Grid, EvolvingDomains.Discretization, EvolvingDomains.Field

    @testset "Cartesian axes and identity" begin
        g = Grid(([0.0, 0.5, 1.0], [-1.0, 2.0]))
        same_coordinates = Grid(([0.0, 0.5, 1.0], [-1.0, 2.0]))
        @test length(g.axes) == 2
        @test g.axes == same_coordinates.axes
        @test g.mesh_id === g.mesh_id
        @test g.mesh_id !== same_coordinates.mesh_id
        @test_throws ArgumentError Grid(([0.0, 0.0, 1.0],))
        @test_throws ArgumentError Grid(([0.0, -1.0],))
        @test_throws ArgumentError Grid((Float64[],))
        @test_throws ArgumentError Grid(([1.0],))
        @test_throws ArgumentError Grid(([0.0, Inf],))
    end

    @testset "Stable level lookup" begin
        g0 = Grid(([0.0, 1.0], [0.0, 1.0]))
        g1 = Grid(([0.0, 0.5, 1.0], [0.0, 1.0]))
        d = Discretization(Dict(:init => g0, :fine => g1))
        @test d.levels[:init] === g0
        @test d.levels[:fine] === g1
        @test_throws KeyError d.levels[:missing]
    end

    @testset "Mutable per-level arrays and validation" begin
        g0 = Grid(([0.0, 1.0], [0.0, 1.0]))
        g1 = Grid(([0.0, 0.5, 1.0], [0.0, 1.0]))
        d = Discretization(Dict(:init => g0, :fine => g1))
        # Layout is method-specific: neither shape is constrained by grid axes.
        u0, u1 = zeros(2, 3), ones(4, 1)
        f = Field(d, Dict(:init => u0, :fine => u1))

        @test f.discretization === d
        @test f.data[:init] === u0
        @test f.data[:fine] === u1
        f.data[:init][2, 1] = 4.0
        f.data[:fine][4, 1] = -2.0
        @test u0[2, 1] == 4.0
        @test u1[4, 1] == -2.0
        @test_throws KeyError f.data[:missing]
        @test_throws ArgumentError Field(d, Dict(:unknown => zeros(1, 1)))
        @test_throws ArgumentError Field(d, Dict(:init => zeros(2), :fine => zeros(1, 2)))
        @test_throws ArgumentError Field(d, Dict(:init => zeros(2, 2), :fine => ones(Int, 2, 2)))

        # Rank is also layout-specific and may exceed the Grid dimension.
        a0, a1 = zeros(2, 1, 4), ones(1, 3, 2)
        high_rank = Field(d, Dict(:init => a0, :fine => a1))
        @test ndims(high_rank.data[:init]) == 3
        @test ndims(high_rank.data[:fine]) == 3
        high_rank.data[:fine][1, 2, 1] = 7.0
        @test a1[1, 2, 1] == 7.0
    end
end
