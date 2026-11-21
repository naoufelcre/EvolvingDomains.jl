# =============================================================================
# VortexReversal3D — CairoMakie visualization only.
#
# Reads what `run.VortexReversal3D.jl` wrote (`params.txt`, `fields.bin`,
# `metrics.csv`, `errors.csv`) and draws figures. It never transports, never
# modifies a field, and contains no physics: every array it shows comes from disk.
#
# CairoMakie 0.15.17 renders neither 3D isosurfaces nor 3D volumes on the Cairo backend
# (`contour!`/`volume!` on a 3D array build a `Volume` plot that the backend drops with a
# warning), so the "3D" panels are orthographic projections of the zero-level point set —
# the φ = 0 crossings of the mesh edges, which is exactly the geometry the simulation
# represents.
#
# Standalone use, on previously produced output:
#   julia --project=validation/env validation/VortexReversal3D/visualize.VortexReversal3D.jl
# =============================================================================
using CairoMakie

"""
    loadRun(dir)

Read one run directory. `fields.bin` is `Float64` records of `nodes³` in column-major order
(first axis fastest) written in the order documented in `writeRun`; the axes are the
uniform `nodes`-node discretization of the box recorded in `params.txt`.
"""
function loadRun(dir)
    params = Dict{String,String}()
    for line in eachline(joinpath(dir, "params.txt"))
        key, value = split(line; limit = 2)
        params[strip(key)] = strip(value)
    end
    n = parse(Int, params["nodes"])
    lo, hi = parse.(Float64, split(params["box"]))
    axis = collect(range(lo, hi; length = n))

    io = open(joinpath(dir, "fields.bin"), "r")
    record() = read!(io, Array{Float64,3}(undef, n, n, n))
    phi0, c0 = record(), record()
    phi = [record() for _ in 1:3]
    c_weno5 = [record() for _ in 1:3]
    c_cip = [record() for _ in 1:3]
    close(io)

    metrics = readCsv(joinpath(dir, "metrics.csv"))
    err = readCsv(joinpath(dir, "errors.csv"))
    return (; dir, nodes = n, axis, phi0, c0, phi, c_weno5, c_cip,
              times = metrics["time"], metrics,
              errors = (phi = (err["l1"][1], err["l2"][1]),
                        c_weno5 = (err["l1"][2], err["l2"][2]),
                        c_cip = (err["l1"][3], err["l2"][3])))
end

function readCsv(path)
    lines = readlines(path)
    header = String.(split(lines[1], ','))
    columns = Dict(name => Float64[] for name in header)
    for line in lines[2:end]
        isempty(strip(line)) && continue
        for (name, field) in zip(header, split(line, ','))
            push!(columns[name], something(tryparse(Float64, field), NaN))
        end
    end
    return columns
end

nearest(axis, x) = argmin(abs.(axis .- x))

"""
    zeroCrossings(axis, phi)

Points where φ changes sign along a mesh edge, linearly interpolated. Read-only: it shells
the isosurface out of the saved field and touches no state.
"""
function zeroCrossings(axis, phi)
    n = length(axis)
    points = NTuple{3,Float64}[]
    for k in 1:n-1, j in 1:n-1, i in 1:n-1
        a = phi[i, j, k]
        for (di, dj, dk) in ((1, 0, 0), (0, 1, 0), (0, 0, 1))
            b = phi[i + di, j + dj, k + dk]
            (a < 0) == (b < 0) && continue
            t = a / (a - b)
            push!(points, (axis[i] + t * di * (axis[i + 1] - axis[i]),
                           axis[j] + t * dj * (axis[j + 1] - axis[j]),
                           axis[k] + t * dk * (axis[k + 1] - axis[k])))
        end
    end
    return points
end

# ------------------------------------------------------------------------- geometry 3D
"""Zero level set of φ at the three snapshot times, as two oblique projections each."""
function figureGeometry(run, path)
    axis = run.axis
    views = (("u = x + 0.6z,  v = y + 0.35z",
              p -> (p[1] + 0.6p[3], p[2] + 0.35p[3])),
             ("u = -y + 0.6z,  v = x + 0.35z",
              p -> (-p[2] + 0.6p[3], p[1] + 0.35p[3])))
    fig = Figure(size = (1450, 1650))
    for (row, k) in enumerate(eachindex(run.phi))
        points = zeroCrossings(axis, run.phi[k])
        for (col, (label, project)) in enumerate(views)
            uv = project.(points)
            ax = Axis(fig[row, col], aspect = DataAspect(),
                      title = string("φ = 0, t = ", round(run.times[k], digits = 2),
                                     "\n", label))
            scatter!(ax, first.(uv), last.(uv); markersize = 3,
                     color = [p[3] for p in points], colormap = :viridis)
        end
    end
    Colorbar(fig[:, 0], label = "z", colormap = :viridis)
    save(path, fig)
    return path
end

# --------------------------------------------------------------------------- slices
"""Orthogonal xy, xz, yz slices through the box centre with the zero-level contour."""
function figureSlices(run, path)
    axis = run.axis
    c = nearest(axis, 0.5)
    fig = Figure(size = (1250, 1150))
    local heat
    for (row, k) in enumerate(eachindex(run.phi))
        A = run.phi[k]
        panels = (("xy (z = 0.5)", A[:, :, c]), ("xz (y = 0.5)", A[:, c, :]),
                  ("yz (x = 0.5)", A[c, :, :]))
        for (col, (name, matrix)) in enumerate(panels)
            ax = Axis(fig[row, col], aspect = DataAspect(),
                      title = string(name, ",  t = ", round(run.times[k], digits = 2)))
            heat = heatmap!(ax, axis, axis, matrix; colormap = :balance,
                            colorrange = (-0.4, 0.4))
            contour!(ax, axis, axis, matrix; levels = [0.0], color = :black, linewidth = 1.5)
        end
    end
    Colorbar(fig[:, 0], heat; label = "φ")
    save(path, fig)
    return path
end

# ---------------------------------------------------------------------- passive scalar
"""Passive scalar on the xy slice through the bump, WENO5 against CIP, at all three times."""
function figureScalar(run, path)
    axis = run.axis
    k = nearest(axis, 0.49)
    peak = max(maximum(maximum, run.c_weno5), maximum(maximum, run.c_cip))
    fig = Figure(size = (1250, 800))
    local heat
    rows = (("WENO5 + SSP-RK3", run.c_weno5), ("CIP (lagrangian)", run.c_cip))
    for (row, (name, series)) in enumerate(rows)
        for (col, kt) in enumerate(eachindex(series))
            ax = Axis(fig[row, col], aspect = DataAspect(),
                      title = string(name, ",  t = ", round(run.times[kt], digits = 2)))
            heat = heatmap!(ax, axis, axis, series[kt][:, :, k]; colormap = :magma,
                            colorrange = (0.0, peak))
        end
    end
    Colorbar(fig[:, 0], heat; label = "c")
    save(path, fig)
    return path
end

# ---------------------------------------------------------------------------- metrics
"""Volume drift over the three snapshots, and the return errors at t = T."""
function figureMetrics(run, path)
    m = run.metrics
    fig = Figure(size = (1300, 480))
    ax = Axis(fig[1, 1], title = "volume drift, $(run.nodes)³",
              xlabel = "t", ylabel = "(V(t) - V(0)) / V(0)")
    for (label, key, colour) in zip(("voxel (grid)", "smeared (grid)", "cut-cell"),
                                    ("voxel", "smeared", "cutcell"),
                                    (:steelblue, :darkorange, :seagreen))
        series = m[key]
        all(isnan, series) && continue
        scatterlines!(ax, run.times, (series .- series[1]) ./ series[1]; label, color = colour,
               marker = :circle, markersize = 9)
    end
    axislegend(ax; position = :lb)
    ax2 = Axis(fig[1, 2], yscale = log10, title = "relative return error at t = T",
               xticks = (1:3, ["φ", "c (WENO5)", "c (CIP)"]))
    for (label, index) in (("L1", 1), ("L2", 2))
        values = [max(e[index], eps()) for e in
                  (run.errors.phi, run.errors.c_weno5, run.errors.c_cip)]
        scatterlines!(ax2, 1:3, values; label, marker = :circle, markersize = 11)
    end
    axislegend(ax2; position = :rt)
    save(path, fig)
    return path
end

# ------------------------------------------------------------------------ convergence
"""Return errors against resolution for every run directory handed in."""
function figureConvergence(runs, path)
    fig = Figure(size = (820, 540))
    ax = Axis(fig[1, 1], xscale = identity, yscale = log10, xlabel = "nodes per axis",
              xticks = [r.nodes for r in runs], ylabel = "relative return error",
              title = "grid convergence of the return error")
    nodes = [r.nodes for r in runs]
    series = (("φ, L1", e -> e.phi[1], :steelblue, :solid),
              ("φ, L2", e -> e.phi[2], :steelblue, :dash),
              ("c WENO5, L1", e -> e.c_weno5[1], :darkorange, :solid),
              ("c WENO5, L2", e -> e.c_weno5[2], :darkorange, :dash),
              ("c CIP, L1", e -> e.c_cip[1], :seagreen, :solid),
              ("c CIP, L2", e -> e.c_cip[2], :seagreen, :dash))
    for (label, getter, colour, style) in series
        scatterlines!(ax, nodes, [max(getter(r.errors), eps()) for r in runs]; label,
               color = colour, linestyle = style, marker = :circle, markersize = 9)
    end
    axislegend(ax; position = :lb)
    save(path, fig)
    return path
end

"""Draw every figure for each run directory; add the convergence figure when two exist."""
function visualizeAll(dirs, outdir)
    runs = loadRun.(dirs)
    for r in runs
        println("figures in ", r.dir)
        figureGeometry(r, joinpath(r.dir, "geometry3d.png"))
        figureSlices(r, joinpath(r.dir, "slices.png"))
        figureScalar(r, joinpath(r.dir, "scalar_schemes.png"))
        figureMetrics(r, joinpath(r.dir, "metrics.png"))
    end
    length(runs) > 1 && figureConvergence(runs, joinpath(outdir, "convergence.png"))
    return runs
end

if abspath(PROGRAM_FILE) == @__FILE__
    outdir = isempty(ARGS) ? joinpath(@__DIR__, "output") : abspath(ARGS[1])
    dirs = filter(d -> isfile(joinpath(d, "params.txt")),
                  [joinpath(outdir, entry) for entry in readdir(outdir)])
    isempty(dirs) && error("no run directory with params.txt under $outdir")
    visualizeAll(sort(dirs), outdir)
end
