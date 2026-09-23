using Test
using EvolvingDomains
using Gridap

grid = CartesianDiscreteModel((0, 1, 0, 1), (1, 1))
geom = EvolvingDiscreteGeometry([-1.0, 1.0, 1.0, -1.0], grid)
io = IOBuffer()
plot(geom; io=io, size=(4, 6))

@test String(take!(io)) == "+----+\n|  ██|\n|██  |\n+----+\n"

plot(geom; io=io, size=(5, 6), label="Iteration 1")
@test startswith(String(take!(io)), "Iterat\n+----+\n")

plot(geom; io=io, size=(5, 6), label="long label\nignored")
@test startswith(String(take!(io)), "long l\n")

tiny = CartesianDiscreteModel((0, 1, 0, 1), (0, 1))
@test_throws ArgumentError plot(EvolvingDiscreteGeometry(tiny); io=io)

plot(geom; io=io, size=(4, 6), field=collect(1.0:4.0))
@test String(take!(io)) == "+----+\n|  ██|\n|██  |\n+----+\n"
@test_throws DimensionMismatch plot(geom; io=io, field=[1.0])

# plot_geometry keeps the legacy behavior through the plot dispatch above.
plot_geometry(geom; io=io, size=(4, 6))
@test String(take!(io)) == "+----+\n|  ██|\n|██  |\n+----+\n"

# Single curve renders a title, a box with data markers, and a legend.
plot_curves([0.0, 1.0, 2.0], [0.0, 1.0, 2.0]; io=io, size=(7, 12))
out = String(take!(io))
@test occursin("+", out) && occursin("*", out)
@test !occursin("\e[", out)
@test startswith(out, "trace")
@test endswith(out, "\n")

# Three series use distinct glyphs, a T(t) title, and a ranged legend.
plot_curves([0.0, 1.0], [[0.0, 1.0], [1.0, 0.0], [0.5, 0.5]]; io=io, size=(9, 30),
    labels=["min", "avg", "max"], ylabel="T")
lines = split(String(take!(io)), '\n')
@test occursin("T(t)", lines[1])
@test any(s -> occursin("*", s) && occursin("+", s) && occursin("x", s), lines) ||
    (any(s -> occursin("*", s), lines) && any(s -> occursin("+", s), lines) && any(s -> occursin("x", s), lines))
@test occursin(r"T, t ∈ \[\s*0\.0\d?, \s*1\.0\d?\]:", lines[end-1])
@test occursin("min", lines[end-1]) && occursin("avg", lines[end-1]) && occursin("max", lines[end-1])

# Ticks use two decimals.
@test any(s -> occursin(r"-?\d+\.\d{2}", s), lines)

# plot(t, y) dispatches to plot_curves.
plot([0.0, 1.0], [0.0, 1.0]; io=io, size=(7, 12))
@test occursin("*", String(take!(io)))

# Dual layout: titles, two boxes side by side, legend indented under curves.
plot(geom, [0.0, 1.0], [[0.0, 1.0], [1.0, 0.0], [0.5, 0.5]]; io=io, size=(11, 60),
    labels=["min", "avg", "max"], ylabel="T")
out = String(take!(io))
@test endswith(out, "\n")
lines = split(out, '\n')
@test occursin("φ ≤ 0", lines[1]) && occursin("T(t)", lines[1])
@test occursin("+  +", lines[2])
@test startswith(lines[end-1], " ") && occursin("t ∈", lines[end-1])
@test lines[end-1] |> s -> occursin("min", s) && occursin("avg", s) && occursin("max", s)

# Empty series draws empty axes instead of throwing.
plot_curves(Float64[], Float64[]; io=io, size=(6, 12))
@test occursin("+", String(take!(io)))

@test_throws DimensionMismatch plot_curves([0.0, 1.0], [0.0]; io=io)
@test_throws DimensionMismatch plot(geom, [0.0, 1.0], [0.0]; io=io)

# Mid ticks of both panels share one row (regression: even heights were off by one).
# Needs a TTY-marked buffer since the colorbar only renders on terminals.
tio = IOContext(IOBuffer(), :terminal => true)
plot(geom, [0.0, 1.0], [[1.0, 4.0]]; io=tio, size=(12, 70), field=collect(1.0:4.0),
    colorrange=(1.0, 4.0), yrange=(1.0, 4.0), labels=["s"])
midlines = filter(l -> occursin("2.50", l), split(String(take!(tio.io)), '\n'))
@test length(midlines) == 1 && length(collect(eachmatch(r"2\.50", midlines[1]))) == 2

# Curve panel width is capped on very wide terminals.
plot([0.0, 1.0], [0.0, 1.0]; io=io, size=(10, 200))
@test length(split(String(take!(io)), '\n')[2]) <= 62
plot(geom, [0.0, 1.0], [0.0, 1.0]; io=io, size=(12, 200), labels=["s"])
top = split(String(take!(io)), '\n')[2]
@test length(split(top, "  ")[2]) <= 62
