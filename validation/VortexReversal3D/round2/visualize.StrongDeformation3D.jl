# =============================================================================
# StrongDeformation3D — CairoMakie visualization only.
#
# Reads what the round-2 driver saved (`params.txt`, `fields.bin`, the CSV
# tables and the per-snapshot `interface_points_t<S>.csv` files) and draws
# figures. No physics, no transport, no diagnostics: every number shown comes
# from disk, and every interface drawing is the saved point set itself — no
# surface is ever reconstructed here.
#
# Standalone use, on previously produced output:
#   julia --project=validation/env validation/VortexReversal3D/round2/visualize.StrongDeformation3D.jl
# =============================================================================
using CairoMakie

const CMAP = :viridis
const SCHEMES = ((1, "WENO5 + SSP-RK3"), (2, "CIP (Lagrangian)"))

# `interface_points_t<S>.csv` drops a trailing ".0": 0.0 -> "t0", 1.0 -> "t1".
tag(t) = replace(string(t), r"\.0$" => "")

"""Header-keyed CSV columns as `Float64`; empty or textual fields become `NaN`."""
function readCsv(path)
    lines = readlines(path)
    header = String.(split(strip(lines[1]), ','))
    cols = Dict(name => Float64[] for name in header)
    for line in lines[2:end]
        isempty(strip(line)) && continue
        for (name, field) in zip(header, split(line, ','))
            push!(cols[name], something(tryparse(Float64, strip(field)), NaN))
        end
    end
    return cols
end

"""`params.txt`: `key value...` lines, values kept as tokens (`times` is a list)."""
function readParams(dir)
    params = Dict{String,Vector{String}}()
    for line in eachline(joinpath(dir, "params.txt"))
        tokens = split(strip(line))
        isempty(tokens) && continue
        params[tokens[1]] = tokens[2:end]
    end
    return params
end

"""One `nodes³` record of `fields.bin`, column-major, first index = x. Records 1 and 2 are
`phi0`/`c0`; snapshot `k` occupies records `3k` (φ), `3k+1` (c WENO5), `3k+2` (c CIP)."""
function readRecord(path, n, record)
    open(path, "r") do io
        seek(io, (record - 1) * n^3 * sizeof(Float64))
        read!(io, Array{Float64,3}(undef, n, n, n))
    end
end

"""Saved interface vertices of one snapshot, as 3D points."""
function readPoints(path)
    isfile(path) || return Point3f[]
    points = Point3f[]
    for line in Iterators.drop(eachline(path), 1)
        fields = split(strip(line), ',')
        length(fields) == 3 || continue
        push!(points, Point3f(parse.(Float64, fields)...))
    end
    return points
end

nearest(axis, x) = argmin(abs.(axis .- x))
label3(vals) = [string(round(v; sigdigits = 3)) for v in vals]

# ------------------------------------------------------- 1. deforming interface
"""Five-time montage of the saved interface point clouds, shared axes and colour scale."""
function figureMontage(dir, times, path)
    clouds = [readPoints(joinpath(dir, "interface_points_t$(tag(t)).csv")) for t in times]
    x = vcat([[p[1] for p in cloud] for cloud in clouds]...)
    span = isempty(x) ? (0.0, 1.0) : extrema(x)
    fig = Figure(size = (2100, 520))
    for (k, t) in enumerate(times)
        ax = Axis3(fig[1, k], aspect = :data, limits = ((0, 1), (0, 1), (0, 1)),
                   elevation = 0.4, azimuth = 0.45, title = "t = $(tag(t))")
        points = clouds[k]
        isempty(points) && (ax.title = "t = $(tag(t))\n(missing points)"; continue)
        xs = [p[1] for p in points]
        scatter!(ax, xs, [p[2] for p in points], [p[3] for p in points];
                 color = xs, colormap = CMAP, colorrange = span, markersize = 2)
    end
    Colorbar(fig[1, length(times) + 1], colorrange = span, colormap = CMAP, label = "x")
    save(path, fig, px_per_unit = 2)
    println("saved ", path)
end

# ------------------------------------------------------------ 2. orthogonal slices
"""xy/xz/yz slices of φ at the middle snapshot, zero level contoured, symmetric colour range."""
function figureSlices(dir, n, times, path)
    axis = range(0, 1; length = n)
    k = nearest(collect(times), 0.5)
    phi = readRecord(joinpath(dir, "fields.bin"), n, 3k)
    scale = maximum(abs, phi)
    scale = iszero(scale) ? eps() : scale
    planes = ((phi[:, :, nearest(axis, 0.56)], "xy, z = 0.56", "x", "y"),
              (phi[:, nearest(axis, 0.48), :], "xz, y = 0.48", "x", "z"),
              (phi[nearest(axis, 0.41), :, :], "yz, x = 0.41", "y", "z"))
    fig = Figure(size = (1950, 640))
    local heat
    for (col, (matrix, name, xl, yl)) in enumerate(planes)
        ax = Axis(fig[1, col], aspect = DataAspect(), xlabel = xl, ylabel = yl,
                  title = "$name,  t = $(tag(times[k]))")
        heat = heatmap!(ax, axis, axis, matrix; colormap = :balance, colorrange = (-scale, scale))
        contour!(ax, axis, axis, matrix; levels = [0.0], color = :black, linewidth = 1.5)
    end
    Colorbar(fig[1, 4], heat, label = "φ")
    save(path, fig, px_per_unit = 2)
    println("saved ", path)
end

# --------------------------------------------------------- 3. scalar schemes
"""WENO5 against CIP on the same z = 0.56 slice at t = 0, T/2, T, one shared scale."""
function figureScalar(dir, n, times, path)
    axis = range(0, 1; length = n)
    kz = nearest(axis, 0.56)
    snapshots = (1, (length(times) + 1) ÷ 2, length(times))
    values = [readRecord(joinpath(dir, "fields.bin"), n, 3 * k + offset)
              for (offset, _) in SCHEMES, k in snapshots]
    peak = iszero(maximum(maximum, values)) ? eps() : maximum(maximum, values)
    fig = Figure(size = (1950, 1000))
    local heat
    for (row, (_, name)) in enumerate(SCHEMES), (col, k) in enumerate(snapshots)
        ax = Axis(fig[row, col], aspect = DataAspect(), xlabel = "x", ylabel = "y",
                  title = "$name,  t = $(tag(times[k]))")
        heat = heatmap!(ax, axis, axis, values[row, col][:, :, kz];
                        colormap = :magma, colorrange = (0, peak))
    end
    Colorbar(fig[:, length(snapshots) + 1], heat, label = "c")
    save(path, fig, px_per_unit = 2)
    println("saved ", path)
end

# ------------------------------------------------------ 4. symmetric difference / drift
"""Interface symmetric difference and relative volume drift over the snapshots."""
function figureMetrics(dir, path)
    m = readCsv(joinpath(dir, "metrics.csv"))
    t = m["time"]
    fig = Figure(size = (1550, 640))
    ax = Axis(fig[1, 1], title = "interface symmetric difference", xlabel = "t",
              ylabel = "symmetric difference (dimensionless)")
    for (key, colour) in (("dsym_smeared", :steelblue), ("dsym_voxel", :darkorange))
        series = m[key]
        all(isnan, series) && continue
        scatterlines!(ax, t, series; label = replace(key, "_" => " "), color = colour,
                      marker = :circle, markersize = 9)
    end
    axislegend(ax; position = :lt)
    ax2 = Axis(fig[1, 2], title = "relative volume drift", xlabel = "t",
               ylabel = "(V(t) − V(0)) / V(0)  (dimensionless)")
    for (key, colour) in (("voxel", :seagreen), ("smeared", :purple), ("cutcell", :goldenrod))
        series = m[key]
        (all(isnan, series) || !isfinite(series[1]) || iszero(series[1])) && continue
        scatterlines!(ax2, t, (series .- series[1]) ./ series[1]; label = key, color = colour,
                      marker = :circle, markersize = 9)
    end
    axislegend(ax2; position = :lt)
    save(path, fig, px_per_unit = 2)
    println("saved ", path)
end

# ------------------------------------------------------- 5. interface distances
"""Forward/backward interface distance error per snapshot, normalised by `a_min`."""
function figureDistance(dir, path)
    d = readCsv(joinpath(dir, "interface.csv"))
    fig = Figure(size = (1550, 640))
    for (col, (name, meankey, maxkey, colour)) in
        enumerate((("forward", "mean_fwd", "max_fwd", :steelblue),
                   ("backward", "mean_bwd", "max_bwd", :darkorange)))
        ax = Axis(fig[1, col], xlabel = "t", ylabel = "distance / a_min",
                  title = "$name distance, normalised by a_min = 0.09")
        scatterlines!(ax, d["time"], d[meankey]; label = "mean", color = colour,
                      marker = :circle, markersize = 9)
        scatterlines!(ax, d["time"], d[maxkey]; label = "max", color = colour,
                      marker = :rect, linestyle = :dash, markersize = 9)
        axislegend(ax; position = :lt)
    end
    save(path, fig, px_per_unit = 2)
    println("saved ", path)
end

# ------------------------------------------------------------- 6. allocations
"""Per-kernel bytes and seconds per step, horizontal bars sharing the operation order."""
function figureAlloc(dir, path)
    operations, bytes_per_node, seconds = String[], Float64[], Float64[]
    for line in Iterators.drop(eachline(joinpath(dir, "allocs.csv")), 1)
        fields = split(strip(line), ',')
        length(fields) == 5 || continue
        push!(operations, String(fields[1]))
        push!(bytes_per_node, something(tryparse(Float64, fields[3]), 0.0))
        push!(seconds, something(tryparse(Float64, fields[4]), 0.0))
    end
    isempty(operations) && return
    fig = Figure(size = (1850, 700))
    for (col, (name, vals)) in enumerate((("bytes / node / step", bytes_per_node),
                                          ("seconds / step", seconds)))
        positive = max.(vals, eps())
        logged = log10(maximum(positive)) - log10(minimum(positive)) > 2
        ax = Axis(fig[1, col], title = name, xlabel = name, ylabel = "kernel",
                  xscale = logged ? log10 : identity,
                  yticks = (1:length(operations), operations), yreversed = true)
        barplot!(ax, positive; direction = :x, color = :steelblue)
        text!(ax, logged ? positive .* 1.12 : positive .+ 0.02 * maximum(positive),
              1:length(operations); text = label3(vals), align = (:left, :center), fontsize = 12)
        xlims!(ax, logged ? minimum(positive) / 2 : 0,
               maximum(positive) * (logged ? 20 : 1.5))
    end
    save(path, fig, px_per_unit = 2)
    println("saved ", path)
end

# --------------------------------------------------------------- convergence
"""Return error against resolution, when `convergence.csv` holds at least two of them."""
function figureConvergence(csv, path)
    nodes, quantity, l1, l2 = Float64[], String[], Float64[], Float64[]
    for line in Iterators.drop(eachline(csv), 1)
        fields = split(strip(line), ',')
        length(fields) == 4 || continue
        push!(nodes, something(tryparse(Float64, fields[1]), NaN))
        push!(quantity, String(fields[2]))
        push!(l1, something(tryparse(Float64, fields[3]), NaN))
        push!(l2, something(tryparse(Float64, fields[4]), NaN))
    end
    resolutions = sort(unique(nodes))
    length(resolutions) < 2 && (println("convergence.csv holds one resolution — skipping ", path);
                               return)
    fig = Figure(size = (1500, 640))
    for (col, (name, errs)) in enumerate((("L1", l1), ("L2", l2)))
        ax = Axis(fig[1, col], xscale = log10, yscale = log10, xlabel = "nodes per axis",
                  ylabel = "relative return error", title = "return error vs resolution, $name",
                  xticks = (resolutions, string.(round.(Int, resolutions))))
        for (q, colour) in (("phi", :steelblue), ("c_weno5", :darkorange), ("c_cip", :seagreen))
            sel = findall(==(q), quantity)
            isempty(sel) && continue
            order = sortperm(nodes[sel])
            scatterlines!(ax, nodes[sel][order], max.(errs[sel][order], eps());
                          label = q, color = colour, marker = :circle, markersize = 10)
        end
        axislegend(ax; position = :lb)
    end
    save(path, fig, px_per_unit = 2)
    println("saved ", path)
end

"""Draw the six per-resolution figures, plus the convergence figure when it applies."""
function visualizeAll(outputdir::AbstractString)
    isdir(outputdir) || (println("no output directory at $outputdir — nothing to draw"); return)
    dirs = sort([joinpath(outputdir, entry) for entry in readdir(outputdir)
                 if isdir(joinpath(outputdir, entry)) && occursin(r"^n\d+$", entry)],
                by = dir -> parse(Int, basename(dir)[2:end]))
    isempty(dirs) && (println("no n<N> resolution directory under $outputdir — nothing to draw");
                      return)
    for dir in dirs
        n = parse(Int, basename(dir)[2:end])
        times = parse.(Float64, readParams(dir)["times"])
        figureMontage(dir, times, joinpath(outputdir, "montage_interface_n$n.png"))
        figureSlices(dir, n, times, joinpath(outputdir, "slices_n$n.png"))
        figureScalar(dir, n, times, joinpath(outputdir, "scalar_schemes_n$n.png"))
        figureMetrics(dir, joinpath(outputdir, "metrics_n$n.png"))
        figureDistance(dir, joinpath(outputdir, "interface_distance_n$n.png"))
        figureAlloc(dir, joinpath(outputdir, "allocation_n$n.png"))
    end
    csv = joinpath(outputdir, "convergence.csv")
    isfile(csv) && figureConvergence(csv, joinpath(outputdir, "convergence.png"))
    return nothing
end

if abspath(PROGRAM_FILE) == @__FILE__
    visualizeAll(joinpath(@__DIR__, "output"))
end
