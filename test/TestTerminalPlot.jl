using Test
using EvolvingDomains
using Gridap

grid = CartesianDiscreteModel((0, 1, 0, 1), (1, 1))
geom = EvolvingDiscreteGeometry([-1.0, 1.0, 1.0, -1.0], grid)
io = IOBuffer()
plot(geom; io=io, size=(4, 6))

@test String(take!(io)) == "+----+\n|  ██|\n|██  |\n+----+"

plot(geom; io=io, size=(5, 6), label="Iteration 1")
@test startswith(String(take!(io)), "Iterat\n+----+\n")

plot(geom; io=io, size=(5, 6), label="long label\nignored")
@test startswith(String(take!(io)), "long l\n")

tiny = CartesianDiscreteModel((0, 1, 0, 1), (0, 1))
@test_throws ArgumentError plot(EvolvingDiscreteGeometry(tiny); io=io)

plot(geom; io=io, size=(4, 6), field=collect(1.0:4.0))
@test String(take!(io)) == "+----+\n|  ██|\n|██  |\n+----+"
@test_throws DimensionMismatch plot(geom; io=io, field=[1.0])
