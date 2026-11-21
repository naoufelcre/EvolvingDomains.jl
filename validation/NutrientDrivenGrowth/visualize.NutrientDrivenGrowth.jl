#!/usr/bin/env julia
#
# visualize.NutrientDrivenGrowth.jl
# =================================
# Figures for the NutrientDrivenGrowth benchmark, from the files written by
# `run.NutrientDrivenGrowth.jl`: `fields.jls` (per-step fields) and `metrics.csv`.
#
#   julia --project=../env --startup-file=no visualize.NutrientDrivenGrowth.jl [tag] [outdir]
#
# Produces, in ./output/<tag>/:
#   fig_domain_evolution.png   montage: two components growing, approaching, merging
#   fig_nutrient.png           nutrient concentration on the physical domain
#   fig_pressure_velocity.png  pressure, velocity magnitude and velocity arrows
#   fig_marker.png             passive material marker
#   fig_metrics.png            area, area-growth rates, component counts vs time
#
# Conventions: the material boundary (level-set zero contour) is always a thick dark line;
# the computational grid boundary (the background box) is a thin dashed grey rectangle; field
# values are masked (NaN) outside the physical domain wherever the field only lives there.

using CairoMakie
using Serialization: deserialize
using Printf: @sprintf, @printf

const OUTDIR = length(ARGS) >= 2 ? ARGS[2] : joinpath(@__DIR__, "output",
    length(ARGS) >= 1 ? ARGS[1] : "n64")

"""Mask `values` outside the physical domain for plotting (NaN hides the entry)."""
function masked(values, phi)
    out = copy(values)
    for k in eachindex(out)
        (isfinite(out[k]) && phi[k] < 0) || (out[k] = NaN)
    end
    return out
end

"""Draw the material boundary (solid) and the computational grid boundary (dashed)."""
function domains!(ax, axes, phi)
    contour!(ax, axes, axes, phi; levels = [0.0], color = (:black, 1.0), linewidth = 2.5,
        label = "material boundary Γ")
    box = [0.0, 1.0, 1.0, 0.0, 0.0]
    lines!(ax, box, box[end:-1:1]; color = (:grey45, 1.0), linewidth = 1.2,
        linestyle = :dash, label = "computational box ∂D")
    return ax
end

function panel(ax, title; aspect = DataAspect())
    ax.title = title
    ax.aspect = aspect
    ax.xgridvisible = false
    ax.ygridvisible = false
    hidespines!(ax)
    return ax
end

"""Evenly spaced frame indices, always including the first and the last frame."""
function select_frames(frames, count)
    count = min(count, length(frames))
    return unique(round.(Int, range(1, length(frames); length = count)))
end

function figure_domain(frames, axes, outdir)
    chosen = select_frames(frames, 8)
    fig = Figure(size = (1500, 700))
    for (k, index) in enumerate(chosen)
        frame = frames[index]
        ax = Axis(fig[1 + (k - 1) ÷ 4, 1 + (k - 1) % 4])
        # Depth of the level set inside the phase (negative), masked outside it.
        depth = masked(min.(frame.phi, 0.0), frame.phi)
        heatmap!(ax, axes, axes, depth; colormap = cgrad(:Blues, 6), colorrange = (-0.35, 0.0),
            nan_color = :transparent)
        domains!(ax, axes, frame.phi)
        panel(ax, @sprintf("t = %.3f   |Ω| = %.4f   %d comp.", frame.t, frame.area,
            frame.components_nodes))
        xlims!(ax, -0.02, 1.02)
        ylims!(ax, -0.02, 1.02)
    end
    Label(fig[0, :], "Nutrient-driven growth: physical domain Ω(t), level-set φ = min(d₁,d₂) − R₀",
        fontsize = 19)
    save(joinpath(outdir, "fig_domain_evolution.png"), fig)
    return fig
end

function figure_nutrient(frames, axes, outdir)
    chosen = select_frames(frames, 6)
    c_range = extrema(reduce(vcat, [frame.c[isfinite.(frame.c)] for frame in frames]))
    fig = Figure(size = (1500, 620))
    for (k, index) in enumerate(chosen)
        frame = frames[index]
        ax = Axis(fig[1 + (k - 1) ÷ 3, 1 + (k - 1) % 3])
        values = masked(frame.c, frame.phi)
        heatmap!(ax, axes, axes, values; colormap = :viridis, colorrange = c_range,
            nan_color = :transparent)
        domains!(ax, axes, frame.phi)
        panel(ax, @sprintf("t = %.3f", frame.t))
        xlims!(ax, 0, 1)
        ylims!(ax, 0, 1)
    end
    Colorbar(fig[1:2, 4], colorrange = c_range, label = "nutrient c", width = 16)
    Label(fig[0, :], "Nutrient concentration on Ω(t): −D_c Δc + λc = 0, c = 1 on Γ",
        fontsize = 19)
    save(joinpath(outdir, "fig_nutrient.png"), fig)
    return fig
end

function figure_pressure_velocity(frames, axes, outdir)
    chosen = select_frames(frames, 4)
    # |v| of the smooth per-cell recovery (the Darcy field seen by the diagnostics), and the
    # nodal transport velocity drawn only in the interface band that actually drives Γ.
    speed_of(frame) = [hypot(frame.cell_vx[min(i, size(frame.cell_vx, 1)),
                                          min(j, size(frame.cell_vx, 2))],
                             frame.cell_vy[min(i, size(frame.cell_vy, 1)),
                                          min(j, size(frame.cell_vy, 2))])
                       for i in eachindex(axes), j in eachindex(axes)]
    Ωmax = maximum(maximum(speed_of(frames[index])[masked(frames[index].p,
        frames[index].phi) .|> isfinite]) for index in chosen)
    ncx = length(axes) - 1
    centres = [(axes[i] + axes[i + 1]) / 2 for i in 1:ncx]
    stride = max(1, ncx ÷ 14)
    fig = Figure(size = (1300, 1150))
    for (k, index) in enumerate(chosen)
        frame = frames[index]
        ax = Axis(fig[1 + (k - 1) ÷ 2, 1 + (k - 1) % 2])
        speed = masked(speed_of(frame), frame.phi)
        heatmap!(ax, axes, axes, speed; colormap = :magma, colorrange = (0.0, Ωmax),
            nan_color = :transparent)
        contour!(ax, axes, axes, masked(frame.p, frame.phi); levels = 5,
            color = (:white, 0.6), linewidth = 1)
        band = [(hypot(frame.vx[i, j], frame.vy[i, j]), i, j) for i in 1:stride:length(axes),
            j in 1:stride:length(axes) if abs(frame.phi[i, j]) < 1.5 / ncx]
        if !isempty(band)
            arrows2d!(ax, [axes[b[2]] for b in band], [axes[b[3]] for b in band],
                [frame.vx[b[2], b[3]] for b in band], [frame.vy[b[2], b[3]] for b in band];
                color = (:black, 0.85), tiplength = 7, tipwidth = 6, shaftwidth = 1.2,
                normalize = false)
        end
        domains!(ax, axes, frame.phi)
        panel(ax, @sprintf("t = %.3f   max|v| on Ω = %.3f", frame.t,
            maximum(speed[isfinite.(speed)])))
        xlims!(ax, 0, 1)
        ylims!(ax, 0, 1)
    end
    Colorbar(fig[1:2, 3], colorrange = (0.0, Ωmax), label = "|v| = |∇p| (per-cell recovery)",
        width = 16)
    Label(fig[0, :], "Darcy velocity and pressure: −Δp = g(c) in Ω, p = 0 on Γ   " *
                     "(arrows: nodal velocity in the interface band)", fontsize = 17)
    save(joinpath(outdir, "fig_pressure_velocity.png"), fig)
    return fig
end

function figure_marker(frames, axes, outdir)
    chosen = select_frames(frames, 6)
    fig = Figure(size = (1500, 620))
    limits = (0.0, maximum(maximum(frame.q) for frame in frames))
    for (k, index) in enumerate(chosen)
        frame = frames[index]
        ax = Axis(fig[1 + (k - 1) ÷ 3, 1 + (k - 1) % 3])
        heatmap!(ax, axes, axes, frame.q; colormap = :balance, colorrange = (0, limits[2]))
        domains!(ax, axes, frame.phi)
        panel(ax, @sprintf("t = %.3f", frame.t))
        xlims!(ax, 0, 1)
        ylims!(ax, 0, 1)
    end
    Colorbar(fig[1:2, 4], colorrange = limits, label = "material marker q", width = 16)
    Label(fig[0, :], "Passive material marker advected by the same velocity as Γ (CIP)",
        fontsize = 19)
    save(joinpath(outdir, "fig_marker.png"), fig)
    return fig
end

function figure_metrics(metrics, outdir)
    t = [row["t"] for row in metrics]
    area = [row["area_cut"] for row in metrics]
    perimeter = [row["perimeter"] for row in metrics]
    predicted = [row["integral_growth"] for row in metrics]
    observed = [row["d_area_dt"] for row in metrics]
    interface_flux = [row["interface_flux"] for row in metrics]
    nodes = [row["components_nodes"] for row in metrics]
    cells = [row["components_cells"] for row in metrics]

    fig = Figure(size = (1400, 950))
    ax1 = Axis(fig[1, 1]; ylabel = "area  |Ω(t)|", xlabel = "time t",
        title = "physical area |Ω(t)| (cut reconstruction)")
    lines!(ax1, t, area; color = :steelblue, linewidth = 2.5)

    ax2 = Axis(fig[1, 2]; ylabel = "perimeter  |Γ(t)|", xlabel = "time t",
        title = "material boundary length |Γ(t)|")
    lines!(ax2, t, perimeter; color = :darkorange, linewidth = 2.5)
    merged = [row["t"] for row in metrics if row["components_nodes"] < 2]
    if !isempty(merged)
        vlines!(ax2, [first(merged)]; color = :grey40, linestyle = :dash,
            label = @sprintf("first single-component step (t = %.3f)", first(merged)))
        axislegend(ax2; position = :lb)
    end

    ax3 = Axis(fig[2, 1]; ylabel = "area growth rate", xlabel = "time t",
        title = "area-growth balance: d|Ω|/dt = ∫_Ω g(c)")
    lines!(ax3, t, predicted; color = :black, linewidth = 2.5, label = "∫_Ω g(c)  (source)")
    lines!(ax3, t, observed; color = :steelblue, linewidth = 2, linestyle = :dash,
        label = "Δ|Ω|/Δt  (reconstructed geometry)")
    lines!(ax3, t, interface_flux; color = :seagreen, linewidth = 2, linestyle = :dot,
        label = "∫_Γ v⋅n  (interface flux)")
    axislegend(ax3; position = :lt)

    ax4 = Axis(fig[2, 2]; ylabel = "connected components", xlabel = "time t",
        title = "topology: 4-connected components on the background grid")
    lines!(ax4, t, Float64.(nodes); color = :black, linewidth = 2.5,
        label = "nodes with φ < 0")
    lines!(ax4, t, Float64.(cells); color = :firebrick, linewidth = 2, linestyle = :dash,
        label = "active cells (in + cut)")
    ax4.yticks = 1:max(maximum(nodes), maximum(cells))
    axislegend(ax4; position = :rt)
    Label(fig[0, :], "NutrientDrivenGrowth diagnostics", fontsize = 19)
    save(joinpath(outdir, "fig_metrics.png"), fig)
    return fig
end

"""Read `metrics.csv` into a vector of `Dict{String,Float64}`."""
function read_metrics(path)
    lines = readlines(path)
    header = split(lines[1], ",")
    rows = Vector{Dict{String,Float64}}()
    for line in lines[2:end]
        isempty(strip(line)) && continue
        fields = split(line, ",")
        length(fields) == length(header) || error("malformed row in $path: $line")
        push!(rows, Dict(header[k] => parse(Float64, fields[k]) for k in eachindex(header)))
    end
    return rows
end

function main(args = ARGS)
    outdir = length(args) >= 2 ? args[2] :
             joinpath(@__DIR__, "output", length(args) >= 1 ? args[1] : "n64")
    data = deserialize(joinpath(outdir, "fields.jls"))
    metrics = read_metrics(joinpath(outdir, "metrics.csv"))
    axes = data.axes
    frames = data.frames
    @printf("%s: %d frames, %d metric rows, %d nodes per axis\n", outdir, length(frames),
        length(metrics), length(axes))
    figure_domain(frames, axes, outdir)
    figure_nutrient(frames, axes, outdir)
    figure_pressure_velocity(frames, axes, outdir)
    figure_marker(frames, axes, outdir)
    figure_metrics(metrics, outdir)
    @printf("figures written to %s\n", outdir)
end

if abspath(PROGRAM_FILE) == @__FILE__
    main()
end
