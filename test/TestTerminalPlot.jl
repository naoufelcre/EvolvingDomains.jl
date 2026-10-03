using Test
using EvolvingDomains
using Gridap
import TPlot

@testset "TPlot geometry integration" begin
    grid = CartesianDiscreteModel((0, 1, 0, 1), (1, 1))
    geom = EvolvingDiscreteGeometry([-1.0, 1.0, 1.0, -1.0], grid)
    phi = reshape(geom.levelset, 2, 2)
    @test EvolvingDomains.plot === TPlot.plot
    @test EvolvingDomains.plot_curves === TPlot.plot_curves

    io = IOBuffer()
    plot(geom; io=io, size=(4, 6))
    @test String(take!(io)) == "+----+\n|  ██|\n|██  |\n+----+\n"
    @test_throws DimensionMismatch plot_geometry(geom; io=io, field=[1.0])
    tiny = CartesianDiscreteModel((0, 1, 0, 1), (0, 1))
    @test_throws ArgumentError plot(EvolvingDiscreteGeometry(tiny); io=io)

    for tty in (false, true), y in ([0.0, 1.0], [0.0 1.0; 1.0 0.0])
        a = IOContext(IOBuffer(), :terminal => tty)
        b = IOContext(IOBuffer(), :terminal => tty)
        plot(geom, [0.0, 1.0], y; io=a, size=(12, 70), field=collect(1.0:4.0))
        TPlot.plot(phi, [0.0, 1.0], y; io=b, size=(12, 70), field=collect(1.0:4.0))
        @test take!(a.io) == take!(b.io)
    end

    panel = TPlot.Geometry(geom)
    scene = TPlot.Row(panel, TPlot.Column(
        TPlot.Curves([0.0, 1.0], [0.0, 1.0]; title="density"),
        TPlot.Curves([0.0, 1.0], [2.0, 3.0]; title="stress")))
    TPlot.render(scene; io=io, size=(20, 70))
    before = String(take!(io))
    @test occursin("density", before) && occursin("stress", before)
    set_levelset!(geom, ones(4))
    TPlot.render(scene; io=io, size=(20, 70))
    @test String(take!(io)) != before
end
