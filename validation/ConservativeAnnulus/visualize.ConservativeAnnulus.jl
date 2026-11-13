# ConservativeAnnulus figures. Reads output/results.jls written by
# run.ConservativeAnnulus.jl and writes five PNGs next to it.
#
#   julia validation/ConservativeAnnulus/visualize.ConservativeAnnulus.jl
#
# Reconstructed geometry (from the Cut) and exact geometry are drawn in different styles
# everywhere: solid blue = what the cut representation says, dashed orange = analytic.

using Pkg
Pkg.activate(joinpath(@__DIR__, "..", "env"); io = devnull)

using CairoMakie
using Serialization
using Printf

const OUTDIR = joinpath(@__DIR__, "output")
const RECONSTRUCTED = :dodgerblue
const EXACT = :orangered

results = deserialize(joinpath(OUTDIR, "results.jls"))
problem = results.problem
quick, fine = results.quick, results.fine

exact_radius(t; inner) = inner ? problem.RIN : problem.ROUT
exact_centre(t) = (problem.C0[1] + t * problem.VEL[1], problem.C0[2] + t * problem.VEL[2])

# ---------------------------------------------------------------- geometry helpers

segment_points(segments) =
    reduce(vcat, [[Point2f(s[1][1], s[1][2]), Point2f(s[2][1], s[2][2])] for s in segments];
           init = Point2f[])

function triangle_points(triangles)
    points = Point2f[]
    for (a, b, c) in triangles
        for (p, q) in ((a, b), (b, c), (c, a))
            push!(points, Point2f(p[1], p[2]))
            push!(points, Point2f(q[1], q[2]))
        end
    end
    return points
end

function exact_circle(t, radius; samples = 400)
    cx, cy = exact_centre(t)
    return [Point2f(cx + radius * cos(θ), cy + radius * sin(θ))
            for θ in range(0, 2π; length = samples)]
end

function draw_grid!(ax, nodes, window)
    lo, hi = window
    inner = nodes[(nodes .> lo[1]) .& (nodes .< hi[1])]
    inner_y = nodes[(nodes .> lo[2]) .& (nodes .< hi[2])]
    for x in inner[1:2:end]
        lines!(ax, [x, x], [lo[2], hi[2]]; color = (:gray, 0.35), linewidth = 0.6)
    end
    for y in inner_y[1:2:end]
        lines!(ax, [lo[1], hi[1]], [y, y]; color = (:gray, 0.35), linewidth = 0.6)
    end
end

# ---------------------------------------------------------------- figure 1: cut geometry

function figure_geometry()
    n = quick.n
    nodes = collect(range(0.0, 1.0; length = n))
    snaps = quick.snapshots
    fig = Figure(size = (1250, 1550))
    for (row, snap) in enumerate(snaps)
        t = snap.t
        cx, cy = exact_centre(t)
        # centre the zoom on the outer interface at 45°, where the cut cells are
        half = 0.05
        zx, zy = cx + 0.7071 * problem.ROUT, cy + 0.7071 * problem.ROUT
        window = ((zx - half, zy - half), (zx + half, zy + half))

        ax = Axis(fig[row, 1]; title = @sprintf("t = %.2f — cut cells on the outer interface", t),
                  xlabel = "x", ylabel = "y", aspect = DataAspect(),
                  limits = (window[1][1], window[2][1], window[1][2], window[2][2]))
        draw_grid!(ax, nodes, window)
        linesegments!(ax, triangle_points(snap.triangles);
                      color = (:dodgerblue, 0.45), linewidth = 0.7,
                      label = "in-phase cut subcells (triangles)")
        linesegments!(ax, segment_points(snap.segments);
                      color = RECONSTRUCTED, linewidth = 2.5,
                      label = "reconstructed interface (Cut)")
        lines!(ax, exact_circle(t, problem.RIN); color = EXACT, linewidth = 2, linestyle = :dash,
               label = "exact interface")
        lines!(ax, exact_circle(t, problem.ROUT); color = EXACT, linewidth = 2, linestyle = :dash)

        ax2 = Axis(fig[row, 2]; title = @sprintf("t = %.2f — full box", t),
                   xlabel = "x", ylabel = "y", aspect = DataAspect())
        linesegments!(ax2, segment_points(snap.segments); color = RECONSTRUCTED, linewidth = 2)
        lines!(ax2, exact_circle(t, problem.RIN); color = EXACT, linewidth = 2, linestyle = :dash)
        lines!(ax2, exact_circle(t, problem.ROUT); color = EXACT, linewidth = 2, linestyle = :dash)
        hidedecorations!(ax2; label = false, ticks = false)
        if row == 1
            axislegend(ax; position = :lt, framevisible = true, labelsize = 11)
        end
    end
    Label(fig[0, :], "Cut geometry: reconstructed (solid blue) vs exact (dashed orange)"; fontsize = 18)
    save(joinpath(OUTDIR, "1_cut_geometry.png"), fig)
    return joinpath(OUTDIR, "1_cut_geometry.png")
end

# ---------------------------------------------------------------- figure 2: density

function figure_density()
    snaps = quick.snapshots
    nodes = collect(range(0.0, 1.0; length = quick.n))
    figure_ = Figure(size = (1500, 520))
    low, high = 0.90, 1.25
    plot = nothing
    for (col, snap) in enumerate(snaps)
        field = copy(snap.rho)
        field[snap.phi .> 0] .= NaN          # outside the material domain: not plotted
        ax = Axis(figure_[1, col]; aspect = DataAspect(),
                  title = @sprintf("remapped ρ, t = %.2f", snap.t), xlabel = "x", ylabel = "y")
        plot = heatmap!(ax, nodes, nodes, field; colormap = :viridis, colorrange = (low, high),
                        nan_color = :transparent)
        linesegments!(ax, segment_points(snap.segments); color = RECONSTRUCTED, linewidth = 1.5)
        lines!(ax, exact_circle(snap.t, problem.RIN); color = EXACT, linewidth = 1.5, linestyle = :dash)
        lines!(ax, exact_circle(snap.t, problem.ROUT); color = EXACT, linewidth = 1.5, linestyle = :dash)
    end
    Colorbar(figure_[1, 4], plot)
    Label(figure_[0, :], "Density on the material domain, masked by the level set (φ > 0 hidden)";
          fontsize = 18)
    save(joinpath(OUTDIR, "2_density.png"), figure_)
    return joinpath(OUTDIR, "2_density.png")
end

# ---------------------------------------------------------------- figure 3: area and centre

function figure_geometry_errors()
    fig = Figure(size = (1250, 480))
    ax = Axis(fig[1, 1]; xlabel = "t", ylabel = "relative area error",
              title = "Annulus area vs exact π(r_out² - r_in²)")
    for (run, color) in ((quick, :black), (fine, :seagreen))
        ts = [r.t for r in run.records]
        lines!(ax, ts, [r.area_error for r in run.records]; color, label = run.label)
        scatter!(ax, ts, [r.area_error for r in run.records]; color, markersize = 5)
    end
    hlines!(ax, [0.0]; color = EXACT, linestyle = :dash)
    axislegend(ax; position = :rb)

    ax2 = Axis(fig[1, 2]; xlabel = "t", ylabel = "|centre - (c₀ + t v)|", yscale = log10,
               title = "Geometric centre error (area-weighted subcell centroids)")
    for (run, color) in ((quick, :black), (fine, :seagreen))
        ts = [r.t for r in run.records]
        lines!(ax2, ts, [r.centre_error for r in run.records]; color, label = run.label)
        scatter!(ax2, ts, [r.centre_error for r in run.records]; color, markersize = 5)
    end
    axislegend(ax2; position = :rb)
    Label(fig[0, :], "Geometry reconstructed from the cut, compared with the analytic annulus";
          fontsize = 18)
    save(joinpath(OUTDIR, "3_geometry_errors.png"), fig)
    return joinpath(OUTDIR, "3_geometry_errors.png")
end

# ---------------------------------------------------------------- figure 4: conservation

function figure_conservation()
    fig = Figure(size = (1300, 500))
    ax = Axis(fig[1, 1]; xlabel = "t", ylabel = "drift in ulps of Σq",
              title = "Nodal sum (the operator's invariant): |Σq(t) - Σq(0)| / (eps·Σq)")
    for (run, color) in ((quick, :black), (fine, :seagreen))
        ts = [r.t for r in run.records]
        scatter!(ax, ts, [abs(r.nodal_drift) / eps(Float64) for r in run.records]; color,
                 markersize = 6, label = run.label)
    end
    hlines!(ax, [0.0]; color = EXACT, linestyle = :dash)
    ylims!(ax, -0.2, 4.5)
    axislegend(ax; position = :rt)

    ax2 = Axis(fig[1, 2]; xlabel = "t", ylabel = "M_h(t)/M(0) - 1",
               title = "Physical mass on the cut geometry")
    for (run, color) in ((quick, :black), (fine, :seagreen))
        ts = [r.t for r in run.records]
        lines!(ax2, ts, [r.mass_error for r in run.records]; color,
               label = "M_h, " * run.label)
        lines!(ax2, ts, [r.mass_exact_field_error for r in run.records];
               color, linestyle = :dot, label = "exact ρ on the same cut, " * run.label)
    end
    hlines!(ax2, [0.0]; color = EXACT, linestyle = :dash)
    axislegend(ax2; position = :lb, labelsize = 10)
    Label(fig[0, :], "Conserved nodal sum ≠ conserved physical mass"; fontsize = 18)
    save(joinpath(OUTDIR, "4_conservation.png"), fig)
    return joinpath(OUTDIR, "4_conservation.png")
end

# ---------------------------------------------------------------- figure 5: transfer

function figure_transfer()
    nf = fine.n
    nodes = collect(range(0.0, 1.0; length = nf))
    fig = Figure(size = (1300, 520))

    ax = Axis(fig[1, 1]; aspect = DataAspect(), xlabel = "x", ylabel = "y",
              title = @sprintf("transfer! :coarse → :fine, error (L∞ = %.2e)", results.transfer.steps[1].error_inf))
    hm = heatmap!(ax, nodes, nodes, results.transfer.error_fine; colormap = :balance,
                  colorrange = (-5e-3, 5e-3))
    Colorbar(fig[1, 2], hm)

    diff = results.prolongation.values .- fine.final_rho
    diff[fine.final_phi .>= 0] .= NaN       # only nodes inside the material are meaningful
    ax2 = Axis(fig[1, 3]; aspect = DataAspect(), xlabel = "x", ylabel = "y",
               title = @sprintf("prolongated coarse ρ - fine ρ at t = 1 (L∞ = %.2e)",
                                results.prolongation.error_inf))
    hm2 = heatmap!(ax2, nodes, nodes, diff; colormap = :balance, colorrange = (-5e-2, 5e-2))
    Colorbar(fig[1, 4], hm2)
    Label(fig[0, :], "Coarse/fine interpolation, measured against the analytic field and between runs";
          fontsize = 18)
    save(joinpath(OUTDIR, "5_transfer.png"), fig)
    return joinpath(OUTDIR, "5_transfer.png")
end

# ---------------------------------------------------------------- run

written = String[]
for f in (figure_geometry, figure_density, figure_geometry_errors, figure_conservation,
          figure_transfer)
    path = f()
    push!(written, path)
    @printf("wrote %s\n", path)
end
