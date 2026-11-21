#!/usr/bin/env julia
"""
Figures for the MeanCurvatureSphere3D benchmark.

Reads only `output/` artifacts written by run.MeanCurvatureSphere3D.jl (timeseries.csv,
snapshots.csv, interface_t*.csv, phi_n*_t*.bin) — no EvolvingDomains, no native library, no
Curvature.jl:

    julia --project=validation/env visualize.MeanCurvatureSphere3D.jl [outputdir]

Figures, all saved as PNG and PDF:

  * `fig_interface3d`  — three-dimensional interface at the first, middle and last snapshot
                         (facet centroids, colored by radius deviation);
  * `fig_slices`       — orthogonal mid-plane slices of phi with the zero contour, the
                         dependable fallback for the 3D view;
  * `fig_radius`       — numerical radius, axis crossings and exact sqrt(R0^2 - 4 mu t);
  * `fig_volume`       — enclosed volume error and radius error over time;
  * `fig_curvature`    — curvature distribution at the final time and mean curvature history.
"""

module MeanCurvatureSphereFigures

using CairoMakie
using Printf

# ---------------------------------------------------------------------------------------
# Readers
# ---------------------------------------------------------------------------------------

"""Read a headered CSV (one leading '#' comment line allowed) into a Dict of column vectors."""
function readTable(path)
    lines = readlines(path)
    first_line = findfirst(l -> !startswith(strip(l), "#") && !isempty(strip(l)), lines)
    isnothing(first_line) && error("empty table: $path")
    header = split(lines[first_line], ",")
    cols = Dict(name => Float64[] for name in header)
    for line in lines[first_line+1:end]
        isempty(strip(line)) && continue
        vals = split(line, ",")
        for (i, name) in enumerate(header)
            push!(cols[name], parse(Float64, vals[i]))
        end
    end
    return cols
end

"""Snapshot list written by the run: tag, time, node count."""
function readSnapshotList(dir)
    tags, times = String[], Float64[]
    n = 0
    for line in readlines(joinpath(dir, "snapshots.csv"))[2:end]
        isempty(strip(line)) && continue
        parts = split(line, ",")
        push!(tags, parts[1])
        push!(times, parse(Float64, parts[2]))
        n = parse(Int, parts[3])
    end
    return (tags=tags, t=times, n=n)
end

function readPhi(dir, n, tag)
    path = joinpath(dir, "phi_n$(n)_t$tag.bin")
    isfile(path) || error("missing level-set snapshot: $path")
    bytes = read(path)
    length(bytes) == 8n^3 || error("$path has $(length(bytes)) bytes, expected $(8n^3)")
    data = reinterpret(Float64, bytes)
    return reshape(data, n, n, n)          # first axis fastest, as written by `write`
end

function readInterface(dir, tag)
    path = joinpath(dir, "interface_t$tag.csv")
    isfile(path) || error("missing interface snapshot: $path")
    return readTable(path)
end

# ---------------------------------------------------------------------------------------
# Figures
# ---------------------------------------------------------------------------------------

function saveFigure(fig, dir, name)
    png = joinpath(dir, name * ".png")
    pdf = joinpath(dir, name * ".pdf")
    save(png, fig)
    save(pdf, fig)
    @info "saved figure" png pdf
    return png
end

"""
    figureInterface3D(dir, tags, rows, h)

Three panels of the interface point cloud (facet centroids) at the first, middle and last
saved time, colored by the deviation of each centroid's radius from the area-weighted mean
radius of that snapshot — the visual counterpart of the sphericity number. A point cloud is
used deliberately: CairoMakie's 3D surface reconstruction is not a dependable guarantee,
whereas scatter of the cut centroids reproduces exactly what the solver produced.
"""
function figureInterface3D(dir, tags, rows, h)
    picks = unique([1, max(1, length(tags) ÷ 2 + 1), length(tags)])
    fig = Figure(; size=(1500, 520))
    limit = h
    for (col, idx) in enumerate(picks)
        tag = tags[idx]
        iface = readInterface(dir, tag)
        x, y, z = iface["x"], iface["y"], iface["z"]
        area = iface["area"]
        r = @. sqrt(x^2 + y^2 + z^2)
        rbar = sum(area .* r) / sum(area)
        ax = Axis3(fig[1, col]; aspect=:data, title=@sprintf("t = %.2f, R = %.4f", rows["t"][idx], rbar),
                   xlabel="x", ylabel="y", zlabel="z",
                   limits=(-0.35, 0.35, -0.35, 0.35, -0.35, 0.35))
        scatter!(ax, x, y, z; color=(r .- rbar) ./ h, colorrange=(-limit, limit),
                 colormap=:balance, markersize=3)
    end
    Colorbar(fig[1, 4]; colormap=:balance, colorrange=(-limit, limit),
             label="(R_p - mean R) / h", width=18)
    return saveFigure(fig, dir, "fig_interface3d")
end

"""
    figureSlices(dir, n, tags)

Mid-plane slices (z = 0 plane, y = 0 plane, x = 0 plane) of the level set for every saved
time, with the zero contour drawn. This is the dependable view of the interface: two-
dimensional `contour!` on a nodal matrix, no reconstructed 3D surface.
"""
function figureSlices(dir, n, tags, rows, h)
    x = range(-0.5, 0.5; length=n)
    mid = n ÷ 2
    levels = [0.0]
    fig = Figure(; size=(1250, 330 * length(tags)))
    for (row, tag) in enumerate(tags)
        phi = readPhi(dir, n, tag)
        planes = ((phi[:, :, mid], "z = 0", 1, 2),
                  (phi[:, mid, :], "y = 0", 1, 3),
                  (phi[mid, :, :], "x = 0", 2, 3))
        for (col, (plane, label, ax1, ax2)) in enumerate(planes)
            ax = Axis(fig[row, col]; aspect=DataAspect(),
                      title=@sprintf("%s slice, t = %.2f", label, rows["t"][row]))
            heatmap!(ax, x, x, plane ./ h; colormap=:balance, colorrange=(-2, 2))
            contour!(ax, x, x, plane; levels=levels, color=:black, linewidth=2)
            hidedecorations!(ax; label=false)
        end
    end
    return saveFigure(fig, dir, "fig_slices")
end

"""    figureRadius(dir, rows, cfg)

Numerical radius against the exact shrinking-radius curve, plus the axis crossings.
"""
function figureRadius(dir, rows)
    t = rows["t"]
    fig = Figure(; size=(1150, 480))
    ax = Axis(fig[1, 1]; xlabel="t", ylabel="radius", title="Numerical vs exact radius")
    lines!(ax, t, rows["exact_radius"]; color=:black, linewidth=2, label="exact sqrt(R0²-4μt)")
    scatterlines!(ax, t, rows["radius"]; color=:crimson, markersize=10,
                  label="mean of cut-surface radii")
    axis_radius = (rows["axis_rx"] .+ rows["axis_ry"] .+ rows["axis_rz"]) ./ 3
    scatter!(ax, t, axis_radius; marker=:utriangle, color=:steelblue, markersize=9,
             label="axis crossings (mean)")
    axislegend(ax; position=:lb)
    lims = (minimum(rows["radius"]) - 0.02, maximum(rows["radius"]) + 0.01)
    ylims!(ax, lims)
    ax2 = Axis(fig[1, 2]; xlabel="t", ylabel="R_num - R_exact", title="Radius error", yscale=log10)
    errors = abs.(rows["radius_err"])
    errors[errors .== 0] .= NaN
    scatterlines!(ax2, t, errors; color=:crimson, markersize=8,
                  label="|cut-surface radius error|")
    axis_error = abs.((rows["axis_rx"] .+ rows["axis_ry"] .+ rows["axis_rz"]) ./ 3 .-
                      rows["exact_radius"])
    axis_error[axis_error .== 0] .= NaN
    scatterlines!(ax2, t, axis_error; color=:steelblue, marker=:utriangle, markersize=8,
                  label="|axis-crossing error|")
    axislegend(ax2; position=:lt)
    return saveFigure(fig, dir, "fig_radius")
end

"""
    figureVolume(dir, rows, h)

Enclosed-volume relative error and radius error (in cells) over time. The radius panel also
carries the area-weighted standard deviation of the interface radius, i.e. the departure
from spherical symmetry.
"""
function figureVolume(dir, rows, h)
    t = rows["t"]
    fig = Figure(; size=(1000, 480))
    ax1 = Axis(fig[1, 1]; xlabel="t", ylabel="relative volume error",
               title="Enclosed volume vs 4πR(t)³/3")
    lines!(ax1, t, 100 .* rows["volume_rel_err"]; color=:darkgreen, linewidth=2)
    hlines!(ax1, [0.0]; color=:gray, linestyle=:dash)
    ax2 = Axis(fig[1, 2]; xlabel="t", ylabel="error / h", title="Radius error and sphericity")
    lines!(ax2, t, rows["radius_err"] ./ h; color=:crimson, linewidth=2,
           label="(R_num - R_exact)/h")
    lines!(ax2, t, rows["radius_std"] ./ h; color=:navy, linewidth=2,
           label="std(R_interface)/h")
    hlines!(ax2, [0.0]; color=:gray, linestyle=:dash)
    axislegend(ax2; position=:lb)
    return saveFigure(fig, dir, "fig_volume")
end

"""
    figureCurvature(dir, rows, tags, h)

Left: distribution of the Curvature.jl facet curvatures at the final time against 2/R.
Right: mean curvature history from both operators against the exact 2/R(t).
"""
function figureCurvature(dir, rows, tags)
    tag = tags[end]
    iface = readInterface(dir, tag)
    kappa = iface["kappa_cv"]
    area = iface["area"]
    r = @. sqrt(iface["x"]^2 + iface["y"]^2 + iface["z"]^2)
    rbar = sum(area .* r) / sum(area)
    kexact = 2 / rbar

    fig = Figure(; size=(1150, 480))
    ax1 = Axis(fig[1, 1]; xlabel="mean curvature κ = κ₁+κ₂",
               ylabel="area-weighted density",
               title=@sprintf("Facet curvature distribution, t = %.2f", rows["t"][end]))
    weights = area ./ sum(area)
    hist!(ax1, kappa; weights=weights, bins=60, color=:steelblue, normalization=:pdf)
    vlines!(ax1, [kexact]; color=:black, linewidth=2, linestyle=:dash,
            label=@sprintf("2/R = %.3f", kexact))
    axislegend(ax1; position=:rt)

    ax2 = Axis(fig[1, 2]; xlabel="t", ylabel="mean curvature",
               title="Average mean curvature vs 2/R(t)")
    lines!(ax2, rows["t"], rows["kappa_exact"]; color=:black, linewidth=2,
           label="2/R_exact")
    scatterlines!(ax2, rows["t"], rows["kappa_cv"]; color=:crimson, markersize=9,
                  label="Curvature.jl (area mean)")
    scatterlines!(ax2, rows["t"], rows["kappa_fd"]; color=:steelblue, marker=:utriangle,
                  markersize=9, label="finite-difference band mean")
    axislegend(ax2; position=:lt)
    return saveFigure(fig, dir, "fig_curvature")
end

"""
    makeFigures(dir)

Load the run's artifacts and save every figure. Called by the run script's entry point and by
this script's own entry point.
"""
function makeFigures(dir)
    rows = readTable(joinpath(dir, "timeseries.csv"))
    snaps = readSnapshotList(dir)
    n, h, tags = snaps.n, 1.0 / (snaps.n - 1), snaps.tags
    figureInterface3D(dir, tags, rows, h)
    figureSlices(dir, n, tags, rows, h)
    figureRadius(dir, rows)
    figureVolume(dir, rows, h)
    figureCurvature(dir, rows, tags)
    return nothing
end

end # module

"""Default output directory: the finest `output/n*` run present, or `output` itself."""
function defaultOutputDir()
    root = joinpath(@__DIR__, "output")
    isdir(root) || return root
    names = filter(n -> occursin(r"^n\d+$", n), readdir(root))
    isempty(names) && return root
    return joinpath(root, last(sort(names; by=n -> parse(Int, n[2:end]))))
end

if abspath(PROGRAM_FILE) == @__FILE__
    dir = isempty(ARGS) ? defaultOutputDir() : abspath(ARGS[1])
    @info "plotting" dir
    MeanCurvatureSphereFigures.makeFigures(dir)
end
