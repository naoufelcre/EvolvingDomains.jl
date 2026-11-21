"""
CurvatureExtensionFigures — figures for the CurvatureExtension3D experiment.

Reads only the files written by `run.CurvatureExtension3D.jl` (`output/n<N>/<geometry>/<method>/`),
so nothing here recomputes the simulation:

    julia --project=../env visualize.CurvatureExtension3D.jl [outputdir]

`makeFigures(root)` is also called at the end of the run script.

Files per run directory: `timeseries.csv`, `snapshots.csv`, `cost.csv`, `facets_<tag>.csv`
(centroid x, y, z, area, kappa from Curvature.jl, kappa from the nodal finite differences),
`phi_<tag>.bin` (Float32, n^3, first axis fastest).
"""
module CurvatureExtensionFigures

using CairoMakie
using Printf
using Statistics: mean, quantile

const COLORS = Dict(:cv => :steelblue, :fd => :darkorange)

"csvLines(path) — non-comment lines of a CSV file, split on commas."
csvLines(path) = [split(l, ',') for l in readlines(path) if !isempty(strip(l)) &&
                                                     !startswith(l, "#")]

"hasLabels(ax) — true when at least one plot in the axis carries a legend label."
hasLabels(ax) = any(p -> p.label[] !== nothing, ax.scene.plots)

"readData(path) — numeric rows of a CSV with a header row; unparsable cells become NaN."
function readData(path)
    isfile(path) || return nothing
    lines = csvLines(path)
    length(lines) < 2 && return nothing
    rows = lines[2:end]
    out = Matrix{Float64}(undef, length(rows), length(rows[1]))
    for (i, r) in enumerate(rows), j in 1:length(rows[1])
        v = j <= length(r) ? tryparse(Float64, r[j]) : nothing
        out[i, j] = v === nothing ? NaN : v
    end
    return out
end

"A column by header name, or `nothing`."
function column(path, name)
    isfile(path) || return nothing
    lines = csvLines(path)
    i = findfirst(==(name), string.(strip.(lines[1])))
    i === nothing && return nothing
    return [something(tryparse(Float64, r[i]), NaN) for r in lines[2:end]]
end

readPhi(path, n) = isfile(path) ?
    reshape(reinterpret(Float32, read(path)), n, n, n) : nothing

readFacets(path) = readData(path)   # x, y, z, area, kappa_cv, kappa_fd

"Run directories present under `root`, as (n, geometry, method, dir)."
function runDirs(root)
    found = Tuple{Int,Symbol,Symbol,String}[]
    for n_dir in readdir(root; join=true)
        m = match(r"^n(\d+)$", basename(n_dir))
        m === nothing && continue
        n = parse(Int, m.captures[1])
        for g in ("A", "B"), meth in ("cv", "fd")
            dir = joinpath(n_dir, g, meth)
            isdir(dir) && push!(found, (n, Symbol(g), Symbol(meth), dir))
        end
    end
    return sort(found, by=t -> (t[1], t[2], t[3]))
end

# ---------------------------------------------------------------------------------------
# 1. Three-dimensional interface evolution
# ---------------------------------------------------------------------------------------
function figInterfaces(root, n)
    dirs = filter(d -> d[1] == n && d[3] == :cv, runDirs(root))
    isempty(dirs) && return
    fig = Figure(size=(1200, 760))
    for (row, (_, g, _, dir)) in enumerate(dirs)
        tags = filter(t -> isfile(joinpath(dir, "facets_$t.csv")),
                      ["t0000", "t0001", "t0002", "t0003", "final"])
        for (col, tag) in enumerate(tags[1:min(3, length(tags))])
            data = readFacets(joinpath(dir, "facets_$tag.csv"))
            data === nothing && continue
            ax = Axis3(fig[row, col], title="$(g) $(tag), Curvature.jl velocity",
                       xlabel="x", ylabel="y", zlabel="z", aspect=:data)
            kappa = data[:, 5]
            kappa = ifelse.(isfinite.(kappa), kappa, NaN)
            lo, hi = quantile(filter(isfinite, kappa), [0.02, 0.98])
            sc = scatter!(ax, data[:, 1], data[:, 2], data[:, 3]; color=kappa,
                          colormap=:balance, colorrange=(lo, hi), markersize=2.0)
            row == 1 && col == 3 && Colorbar(fig[row, col + 1], sc, label="kappa (1/length)")
            limits!(ax, -0.35, 0.35, -0.35, 0.35, -0.35, 0.35)
        end
    end
    Label(fig[0, 1:3], "Interface evolution, facet colour = Curvature.jl mean curvature, n = $n",
          fontsize=18)
    save(joinpath(root, "n$n", "fig_interfaces_3d.png"), fig)
    return fig
end

# ---------------------------------------------------------------------------------------
# 2. Orthogonal cross-sections of the level set
# ---------------------------------------------------------------------------------------
function figCrossSections(root, n)
    fig = Figure(size=(1300, 780))
    for (row, g) in enumerate((:A, :B))
        dir = joinpath(root, "n$n", string(g), "cv")
        phi = readPhi(joinpath(dir, "phi_final.bin"), n)
        phi === nothing && continue
        x = collect(range(-0.5, 0.5; length=n))
        mid = n ÷ 2 + 1
        planes = ((phi[:, :, mid], "z = 0 plane"), (phi[:, mid, :], "y = 0 plane"),
                  (phi[mid, :, :], "x = 0 plane"))
        for (col, (slice, label)) in enumerate(planes)
            ax = Axis(fig[row, col], title="geometry $g: $label", xlabel="axis 1",
                      ylabel="axis 2", aspect=DataAspect())
            hm = heatmap!(ax, x, x, permutedims(slice); colormap=:greys, colorrange=(-0.2, 0.2))
            contour!(ax, x, x, permutedims(slice); levels=[0.0], color=:red, linewidth=2)
            row == 1 && col == 3 && Colorbar(fig[row, col + 1], hm, label="phi")
        end
    end
    Label(fig[0, 1:3], "Orthogonal cross-sections at t_end, Curvature.jl-driven run, n = $n",
          fontsize=18)
    save(joinpath(root, "n$n", "fig_cross_sections.png"), fig)
    return fig
end

# ---------------------------------------------------------------------------------------
# 3. Curvature distributions on the interface
# ---------------------------------------------------------------------------------------
function figCurvature(root, n)
    fig = Figure(size=(1300, 800))
    for (row, g) in enumerate((:A, :B))
        dir = joinpath(root, "n$n", string(g), "cv")
        for (col, tag) in enumerate(("t0000", "final"))
            data = readFacets(joinpath(dir, "facets_$tag.csv"))
            data === nothing && continue
            ax = Axis(fig[row, col], title="$(g) $(tag): interface curvature",
                      xlabel="kappa", ylabel="facet count")
            kcv = filter(isfinite, data[:, 5])
            kfd = filter(isfinite, data[:, 6])
            hist!(ax, kcv; bins=40, color=(COLORS[:cv], 0.55),
                  label="Curvature.jl (area-weighted mean $(round(sum(data[:, 4] .* data[:, 5]) / sum(data[:, 4]); digits=3)))")
            hist!(ax, kfd; bins=40, color=(COLORS[:fd], 0.45), label="nodal FD")
            hasLabels(ax) && axislegend(ax; position=:lt, labelsize=10)
        end
    end
    ax = Axis(fig[3, 1:3], title="Curvature.jl vs nodal FD at the same facet centroids (t_end)",
              xlabel="kappa, Curvature.jl", ylabel="kappa, nodal FD", aspect=DataAspect())
    for g in (:A, :B)
        data = readFacets(joinpath(root, "n$n", string(g), "cv", "facets_final.csv"))
        data === nothing && continue
        scatter!(ax, data[:, 5], data[:, 6]; markersize=4, color=(g === :A ? :black : :purple),
                 label="geometry $g")
    end
    hasLabels(ax) && axislegend(ax; position=:lt)
    save(joinpath(root, "n$n", "fig_curvature_distribution.png"), fig)
    return fig
end

# ---------------------------------------------------------------------------------------
# 4. Sphere radius against the exact solution
# ---------------------------------------------------------------------------------------
function figRadius(root, n)
    fig = Figure(size=(1200, 760))
    ax = Axis(fig[1, 1], title="Geometry A: area-weighted interface radius", xlabel="t",
              ylabel="R")
    ax2 = Axis(fig[2, 1], title="Radius error against the exact solution", xlabel="t",
               ylabel="R_num - R_exact")
    ts = collect(range(0.0, 0.12; length=100))
    lines!(ax, ts, sqrt.(max.(0.23^2 .- 4 * 0.01 .* ts, 0.0)); color=:black,
           linestyle=:dash, linewidth=2, label="exact sqrt(R0^2 - 4 mu t)")
    for meth in (:cv, :fd)
        path = joinpath(root, "n$n", "A", string(meth), "timeseries.csv")
        t = column(path, "t")
        t === nothing && continue
        r = column(path, "radius")
        rex = column(path, "radius_exact")
        lines!(ax, t, r; color=COLORS[meth], linewidth=2,
               label=meth === :cv ? "Curvature.jl + extension" : "nodal finite differences")
        good = isfinite.(rex)
        lines!(ax2, t[good], (r .- rex)[good]; color=COLORS[meth], linewidth=2)
        scatter!(ax2, t[good], (r .- rex)[good]; color=COLORS[meth], markersize=6)
    end
    hasLabels(ax) && axislegend(ax; position=:rt)
    hlines!(ax2, [0.0]; color=:black, linestyle=:dot)
    for r in 1:2
        rowsize!(fig.layout, r, Auto(1.0))
    end
    colsize!(fig.layout, 1, Auto(1.0))
    Label(fig[0, 1], "Sphere control, n = $n", fontsize=18)
    save(joinpath(root, "n$n", "fig_radius_vs_exact.png"), fig)
    return fig
end

# ---------------------------------------------------------------------------------------
# 5. Perturbed-shape anisotropy
# ---------------------------------------------------------------------------------------
function figAnisotropy(root, n)
    fig = Figure(size=(1200, 760))
    ax = Axis(fig[1, 1], title="Geometry B: axis crossings", xlabel="t", ylabel="radius")
    ax2 = Axis(fig[2, 1], title="Anisotropy (R_z - R_equator) and sphericity",
               xlabel="t", ylabel="length / relative")
    for meth in (:cv, :fd)
        path = joinpath(root, "n$n", "B", string(meth), "snapshots.csv")
        t = column(path, "t")
        t === nothing && continue
        cz = column(path, "crossing_z")
        cx = column(path, "crossing_x")
        cy = column(path, "crossing_y")
        lines!(ax, t, cz; color=COLORS[meth], linewidth=2, label="$meth: R along z")
        lines!(ax, t, (cx .+ cy) ./ 2; color=COLORS[meth], linewidth=2, linestyle=:dash,
               label="$meth: R in the equator")
        aniso = cz .- (cx .+ cy) ./ 2
        lines!(ax2, t, aniso; color=COLORS[meth], linewidth=2, label="anisotropy, $meth")
        spher = column(joinpath(root, "n$n", "B", string(meth), "timeseries.csv"),
                       "sphericity")
        tt = column(joinpath(root, "n$n", "B", string(meth), "timeseries.csv"), "t")
        lines!(ax2, tt, spher; color=COLORS[meth], linestyle=:dot, linewidth=2,
               label="sphericity, $meth")
    end
    hasLabels(ax) && axislegend(ax; position=:rt)
    hasLabels(ax2) && axislegend(ax2; position=:rt, labelsize=10)
    for r in 1:2
        rowsize!(fig.layout, r, Auto(1.0))
    end
    colsize!(fig.layout, 1, Auto(1.0))
    Label(fig[0, 1], "Perturbed sphere (P2, eps = 0.15), n = $n", fontsize=18)
    save(joinpath(root, "n$n", "fig_anisotropy.png"), fig)
    return fig
end

# ---------------------------------------------------------------------------------------
# 6. Curvature-driven versus finite-difference-driven evolution
# ---------------------------------------------------------------------------------------
function figCompare(root, n)
    fig = Figure(size=(1200, 760))
    ax = Axis(fig[1, 1], title="Radius: Curvature.jl-driven minus FD-driven", xlabel="t",
              ylabel="R_cv - R_fd")
    ax2 = Axis(fig[2, 1], title="Surface area (both methods)", xlabel="t", ylabel="area")
    for g in (:A, :B)
        cv = readData(joinpath(root, "n$n", string(g), "cv", "timeseries.csv"))
        fd = readData(joinpath(root, "n$n", string(g), "fd", "timeseries.csv"))
        (cv === nothing || fd === nothing) && continue
        m = min(size(cv, 1), size(fd, 1))
        # columns: step, t, nfacets, radius, radius_exact, volume, volume_exact, area, sphericity
        lines!(ax, cv[1:m, 2], cv[1:m, 4] .- fd[1:m, 4]; linewidth=2,
               label="geometry $g", color=(g === :A ? :black : :purple))
        lines!(ax2, cv[1:m, 2], cv[1:m, 8]; linewidth=2, color=(g === :A ? :black : :purple),
               label="geometry $g, cv")
        lines!(ax2, fd[1:m, 2], fd[1:m, 8]; linewidth=2, linestyle=:dash,
               color=(g === :A ? :black : :purple), label="geometry $g, fd")
    end
    hasLabels(ax) && axislegend(ax; position=:lt)
    hasLabels(ax2) && axislegend(ax2; position=:rt, labelsize=9)
    for r in 1:2
        rowsize!(fig.layout, r, Auto(1.0))
    end
    colsize!(fig.layout, 1, Auto(1.0))
    Label(fig[0, 1], "Cross-package versus nodal-finite-difference evolution, n = $n",
          fontsize=18)
    save(joinpath(root, "n$n", "fig_cv_vs_fd.png"), fig)
    return fig
end

# ---------------------------------------------------------------------------------------
# 7. Cost breakdown
# ---------------------------------------------------------------------------------------
function figCost(root, n)
    fig = Figure(size=(1200, 700))
    ax = Axis(fig[1, 1], title="Per-step cost breakdown, n = $n (fraction of the step)",
              xlabel="run (absolute seconds per step below)", ylabel="fraction of step time",
              width=1050, height=430)
    comps = (:cut, :curvature, :extension, :fd, :transport, :reinit)
    colors = [(COLORS[:cv], 0.35), COLORS[:cv], :seagreen, :goldenrod, :grey, :firebrick]
    labels, mats, totals = String[], Vector{Float64}[], Float64[]
    for (g, meth) in ((:A, :cv), (:B, :cv), (:A, :fd), (:B, :fd))
        path = joinpath(root, "n$n", string(g), string(meth), "cost.csv")
        isfile(path) || continue
        per_step = readData(path)[:, 4] ./ 1e3          # ms -> s
        push!(labels, "$g/$meth")
        push!(mats, per_step[1:length(comps)])
        push!(totals, per_step[length(comps) + 1])
    end
    isempty(mats) && return
    # Fractions, stacked with `fillto` (heights are cumulative): the absolute cost of a
    # Curvature.jl step is 2.2-2.6 s against 0.6 s for a finite-difference step, so a
    # linear stack would show one colour only.
    for (r, per_step) in enumerate(mats)
        frac = per_step ./ totals[r]
        cum = 0.0
        for (i, c) in enumerate(comps)
            h = frac[i]
            h <= 0 && continue
            barplot!(ax, [r], [cum + h]; fillto=[cum], color=colors[i], width=0.55,
                     label=r == 1 ? string(c) : nothing)
            cum += h
        end
    end
    ax.xticks = (1:length(labels), [string(l, "\n", round(t; digits=2), " s/step")
                                    for (l, t) in zip(labels, totals)])
    hasLabels(ax) && axislegend(ax; position=:rt, labelsize=10)
    Label(fig[0, 1], "Where the time goes: Curvature.jl dominates", fontsize=18)
    save(joinpath(root, "n$n", "fig_cost.png"), fig)
    return fig
end

# ---------------------------------------------------------------------------------------
# 8. Extension self-test convergence (Stage C)
# ---------------------------------------------------------------------------------------
function figExtensionTest(root)
    path = joinpath(root, "stageC", "extension_test.csv")
    isfile(path) || return
    data = readData(path)
    fig = Figure(size=(1000, 500))
    ax = Axis(fig[1, 1], title="Stage C: extension error on a known surface function",
              xlabel="h", ylabel="error", xscale=log10, yscale=log10)
    h = data[:, 2]
    scatterlines!(ax, h, data[:, 7]; label="|g_x - x*|max", color=:steelblue)
    scatterlines!(ax, h, data[:, 8]; label="|g_x - x*| rms", color=:darkorange)
    scatterlines!(ax, h, data[:, 10]; label="|v - n_analytic|max", color=:black)
    axislegend(ax; position=:rb)
    Label(fig[0, 1], "Closest-point extension: measured, first-order accurate", fontsize=16)
    save(joinpath(root, "stageC", "fig_extension_test.png"), fig)
    return fig
end

"makeFigures(root) — every figure of this experiment, from the saved outputs only."
function makeFigures(root)
    files = runDirs(root)
    isempty(files) && (@warn "no run directories under $root"; return)
    made = String[]
    figExtensionTest(root)
    for n in sort(unique(d[1] for d in files))
        figInterfaces(root, n)
        figCrossSections(root, n)
        figCurvature(root, n)
        figRadius(root, n)
        figAnisotropy(root, n)
        figCompare(root, n)
        figCost(root, n)
        append!(made, filter(f -> endswith(f, ".png"), readdir(joinpath(root, "n$n"); join=true)))
    end
    println("figures written: ", join(basename.(made), ", "))
    return made
end

end # module

if abspath(PROGRAM_FILE) == @__FILE__
    root = isempty(ARGS) ? joinpath(@__DIR__, "output") : ARGS[1]
    CurvatureExtensionFigures.makeFigures(root)
end
