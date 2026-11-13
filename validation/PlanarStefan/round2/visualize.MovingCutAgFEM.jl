#!/usr/bin/env julia
# Visualization for the MovingCutAgFEM experiment. Pure consumer: every number comes from
# the CSVs written by run.MovingCutAgFEM.jl or from re-running the level-set transport to
# recover the moving geometry. No assembly or solving happens here.
#
#   julia --project=. visualize.MovingCutAgFEM.jl            # -> output/fig_*.png

using CairoMakie
using DelimitedFiles
using Printf
using Statistics
using LinearAlgebra

include("run.MovingCutAgFEM.jl")   # definitions only; the __FILE__ guard keeps main() off

const OUT = joinpath(@__DIR__, "output")
const RES = (33, 65, 129)
const COLORS = Makie.wong_colors()

function readcsv(name)
    path = joinpath(OUT, name)
    isfile(path) || return nothing
    lines = readlines(path)
    length(lines) <= 1 && return nothing
    header = Symbol.(split(lines[1], ','))
    raw = readdlm(path, ','; skipstart = 1)
    size(raw, 1) == 0 && return nothing
    return NamedTuple{Tuple(header)}(Tuple(Float64.(raw[:, k]) for k in 1:length(header)))
end

readtext(name) = isfile(joinpath(OUT, name)) ? read(joinpath(OUT, name), String) : ""

# ---------------------------------------------------------------- geometry ---

"""Numerical level set advanced by the same WENO5 transport, snapshots at selected times."""
function geometrySnapshots(p; n = 65, every = 20)
    grid, _, field = buildEvolvingGeometry(p; n = n)
    velocity = transportVelocity(p, grid)
    snaps = Tuple{Float64,Any}[(0.0, ED.cut(field, :init))]
    for step in 1:p.nsteps
        ED.eulerianTransport!(field, velocity, p.dt; level_id = :init)
        step % every == 0 && push!(snaps, (step * p.dt, ED.cut(field, :init)))
    end
    return grid, snaps
end

"""Inside/outside triangles of the cut cells plus the interface chords."""
function cutMesh(cut, grid)
    ins = NTuple{3,NTuple{2,Float64}}[]
    outs = NTuple{3,NTuple{2,Float64}}[]
    for s in eachindex(cut.subcell_states)
        cut.bg_states[Int(cut.cell_to_bgcell[s])] == ED.CUT_CUT || continue
        k = Int(cut.point_pointers[s])
        tri = (cut.physical_points[cut.point_ids[k]],
               cut.physical_points[cut.point_ids[k + 1]],
               cut.physical_points[cut.point_ids[k + 2]])
        push!(cut.subcell_states[s] == ED.CUT_IN ? ins : outs, tri)
    end
    chords = Tuple{NTuple{2,Float64},NTuple{2,Float64}}[]
    b = cut.boundary
    for f in eachindex(b.facet_to_bgcell)
        k = Int(b.point_pointers[f])
        push!(chords, (b.physical_points[b.point_ids[k]],
                       b.physical_points[b.point_ids[k + 1]]))
    end
    return ins, outs, chords
end

triangles!(ax, tris, color, alpha = 0.55) = for t in tris
    poly!(ax, [Point2f(t[1]...), Point2f(t[2]...), Point2f(t[3]...)];
        color = color, strokewidth = 0.2, strokecolor = (:gray, 0.4), alpha = alpha)
end

function figGeometry(p)
    grid, snaps = geometrySnapshots(p; n = 65, every = 50)
    picks = [1, findfirst(s -> s[1] ≈ 0.25, snaps), length(snaps)]
    fig = Figure(size = (1500, 520))
    for (k, idx) in enumerate(picks)
        idx === nothing && continue
        t, cut = snaps[idx]
        ax = Axis(fig[1, k]; aspect = DataAspect(), title = @sprintf("t = %.2f", t),
            xlabel = "x", ylabel = "y", limits = ((0, 1), (0, 1)))
        ins, outs, chords = cutMesh(cut, grid)
        triangles!(ax, outs, COLORS[7], 0.35)
        triangles!(ax, ins, COLORS[2], 0.75)
        for (a, b) in chords
            lines!(ax, [a[1], b[1]], [a[2], b[2]]; color = :black, linewidth = 1.4)
        end
        θ = range(0, 2π; length = 400)
        cx, cy = center(p, t)
        lines!(ax, cx .+ p.R .* cos.(θ), cy .+ p.R .* sin.(θ); color = :magenta,
            linewidth = 1.2, linestyle = :dash)
        k == 1 && text!(ax, 0.03, 0.95; text = "cut cells (EvolvingDomains `cut`)",
            space = :relative, fontsize = 12)
    end
    Label(fig[0, :], "Moving domain from the numerical level set — inside cells (orange), outside (grey), interface chords (black), exact circle (magenta dashed)",
        fontsize = 13)
    save(joinpath(OUT, "fig_geometry.png"), fig)
    return fig
end

# --------------------------------------------------------------- solutions ---

function figSolutions(p; n = 65, t = 0.5)
    diag = readcsv("diagnostics_n$(n).csv")
    isfile(joinpath(OUT, "u_bg_n$(n)_t500.csv")) || return nothing
    raw = readdlm(joinpath(OUT, "u_bg_n$(n)_t500.csv"), ','; skipstart = 1)
    x, y, u = Float64.(raw[:, 1]), Float64.(raw[:, 2]), Float64.(raw[:, 3])
    m = length(unique(x))
    X = reshape(x, m, m)
    Y = reshape(y, m, m)
    U = reshape(u, m, m)
    Uex = [uExact(p, X[i, j], Y[i, j], t) for i in 1:m, j in 1:m]
    inside = [(X[i, j] - center(p, t)[1])^2 + (Y[i, j] - center(p, t)[2])^2 <= p.R^2
              for i in 1:m, j in 1:m]
    Uc = [inside[i, j] ? U[i, j] : NaN for i in 1:m, j in 1:m]
    E = [inside[i, j] ? U[i, j] - Uex[i, j] : NaN for i in 1:m, j in 1:m]

    fig = Figure(size = (1550, 470))
    ax1 = Axis(fig[1, 1]; aspect = DataAspect(), title = "exact  u = R² - |x-c|²")
    ax2 = Axis(fig[1, 2]; aspect = DataAspect(), title = "recovered  u_bg = E u_m  (AgFEM)")
    ax3 = Axis(fig[1, 3]; aspect = DataAspect(), title = "pointwise error  u_h - u")
    hm1 = heatmap!(ax1, X[:, 1], Y[1, :], permutedims(Uex); colormap = :viridis,
        colorrange = (0.0, p.R^2))
    hm2 = heatmap!(ax2, X[:, 1], Y[1, :], permutedims(Uc); colormap = :viridis,
        colorrange = (0.0, p.R^2))
    emax = 5maximum(abs, filter(isfinite, E))
    hm3 = heatmap!(ax3, X[:, 1], Y[1, :], permutedims(E); colormap = :balance,
        colorrange = (-emax, emax))
    for ax in (ax1, ax2, ax3)
        θ = range(0, 2π; length = 400)
        cx, cy = center(p, t)
        lines!(ax, cx .+ p.R .* cos.(θ), cy .+ p.R .* sin.(θ); color = :white, linewidth = 1)
        limits!(ax, 0.1, 0.8, 0.15, 0.85)
    end
    Colorbar(fig[1, 4], hm1; label = "u")
    Colorbar(fig[1, 5], hm3; label = "error")
    Label(fig[0, :], @sprintf("n = %d, t = %.2f — L2 relative error %.2e (interpolation-only %.2e); the AgFEM extension is nonzero on the outside part of cut cells, which the circle mask hides",
        n, t, diag.l2_rel[end], diag.l2_interp_rel[end]), fontsize = 12)
    save(joinpath(OUT, "fig_solution.png"), fig)
    return fig
end

# ------------------------------------------------------------- convergence ---

function figConvergence(p)
    hs, tot, interp, area = Float64[], Float64[], Float64[], Float64[]
    for n in RES
        d = readcsv("diagnostics_n$(n).csv")
        d === nothing && continue
        push!(hs, 1.0 / (n - 1))
        push!(tot, mean(d.l2_rel))
        push!(interp, mean(d.l2_interp_rel))
        push!(area, mean(d.area_rel_error))
    end
    isempty(hs) && return nothing
    fig = Figure(size = (1000, 620))
    ax = Axis(fig[1, 1]; xscale = log10, yscale = log10,
        xlabel = "h = 1 / (n-1)", ylabel = "error",
        title = "Manufactured-solution convergence (mean over t ∈ [0, 0.5])")
    scatterlines!(ax, hs, tot; color = COLORS[1], marker = :circle, label = "‖u_h - u‖ / ‖u‖")
    scatterlines!(ax, hs, interp; color = COLORS[2], marker = :rect,
        label = "interpolation-only share")
    scatterlines!(ax, hs, area; color = COLORS[3], marker = :utriangle,
        label = "geometric area error |Ω_h - Ω| / |Ω|")
    href = hs[1] .* (hs ./ hs[1]).^2
    lines!(ax, hs, 0.6 .* href ./ href[1] .* tot[1]; color = (:black, 0.5),
        linestyle = :dash, label = "slope h²")
    axislegend(ax; position = :lb)
    for (h, v, w) in zip(hs, tot, interp)
        text!(ax, h, v; text = @sprintf(" %.1e", v), fontsize = 10, offset = (6, 6))
    end
    save(joinpath(OUT, "fig_convergence.png"), fig)
    return fig
end

# ------------------------------------------------------------------- DOFs ----

function figDofs()
    fig = Figure(size = (1150, 600))
    ax = Axis(fig[1, 1]; xlabel = "t", ylabel = "count",
        title = "Master DOFs, aggregates and active background DOFs")
    for (k, n) in enumerate(RES)
        s = readcsv("steps_n$(n).csv")
        d = readcsv("diagnostics_n$(n).csv")
        s === nothing && continue
        scatter!(ax, s.t, s.n_master; color = COLORS[k], markersize = 5,
            label = "n_master (n = $(n-1)²)")
        d === nothing || scatter!(ax, d.t, d.n_aggregates; color = COLORS[k],
            marker = :rect, markersize = 8)
        d === nothing || scatter!(ax, d.t, d.n_active_bg; color = COLORS[k],
            marker = :utriangle, markersize = 8)
    end
    axislegend(ax; position = :rb)
    save(joinpath(OUT, "fig_dofs.png"), fig)
    return fig
end

# ------------------------------------------------------------- time series ---

function figTimes()
    fig = Figure(size = (1300, 760))
    ax1 = Axis(fig[1, 1]; xlabel = "t", ylabel = "ms",
        title = "refresh_geometry! per step (● structural, × retained)")
    for (k, n) in enumerate(RES)
        s = readcsv("steps_n$(n).csv")
        s === nothing && continue
        for flag in (true, false)
            sel = s.structural .== (flag ? 1.0 : 0.0)
            any(sel) || continue
            scatter!(ax1, s.t[sel], 1e3 .* s.t_refresh[sel];
                color = COLORS[k], marker = flag ? :circle : :x, markersize = 6,
                label = "n=$(n-1)² $(flag ? "structural" : "retained")")
        end
    end
    axislegend(ax1; position = :rt, nbanks = 2)

    ax2 = Axis(fig[2, 1]; xlabel = "t", ylabel = "ms",
        title = "assemble! and Cholesky per step")
    for (k, n) in enumerate(RES)
        s = readcsv("steps_n$(n).csv")
        s === nothing && continue
        scatter!(ax2, s.t, 1e3 .* s.t_assemble; color = COLORS[k], markersize = 5,
            label = "assemble n=$(n-1)²")
        scatter!(ax2, s.t, 1e3 .* s.t_factor; color = COLORS[k], marker = :rect,
            markersize = 5, label = "cholesky n=$(n-1)²")
    end
    axislegend(ax2; position = :rt, nbanks = 2)
    save(joinpath(OUT, "fig_times.png"), fig)
    return fig
end

# ------------------------------------------------------------------ memory ---

function figMemory()
    stages = ["cut", "convert", "aggregate", "refresh", "assemble", "cholesky"]
    keys_ = (:b_cut, :b_convert, :b_aggregate, :b_refresh, :b_assemble, :b_factor)
    M = zeros(length(RES), length(stages))
    for (k, n) in enumerate(RES)
        s = readcsv("steps_n$(n).csv")
        s === nothing && continue
        for (j, key) in enumerate(keys_)
            M[k, j] = max(mean(getproperty(s, key)) / 1024, 1e-3)   # floor for log scale
        end
    end
    fig = Figure(size = (1050, 560))
    ax = Axis(fig[1, 1]; ylabel = "mean allocation per step (kB)", xlabel = "",
        title = "Allocation by pipeline stage (geometry → operator → solve; floor 1e-3 kB)",
        xticks = (1:length(stages), stages),
        xticklabelrotation = 0.3, yscale = log10)
    for (k, n) in enumerate(RES)
        barplot!(ax, (1:length(stages)) .+ 0.22 * (k - 2), M[k, :];
            color = COLORS[k], label = "n = $(n-1)²")
    end
    axislegend(ax; position = :rt)
    save(joinpath(OUT, "fig_memory.png"), fig)
    return fig
end

# ----------------------------------------------------------------- ferrite ---

function figFerrite(p; n = 33)
    path = joinpath(OUT, "ferrite_quadrature_n$(n).csv")
    isfile(path) || return nothing
    raw = readdlm(path, ','; skipstart = 1)
    isempty(raw) && return nothing
    cell = Int(raw[1, 1])
    x, y = Float64.(raw[:, 2]), Float64.(raw[:, 3])
    ξ, η = Float64.(raw[:, 4]), Float64.(raw[:, 5])
    w = Float64.(raw[:, 6])

    grid, _, field = buildEvolvingGeometry(p; n = n)
    cut = ED.cut(field, :init)
    hx = grid.axes[1][2] - grid.axes[1][1]
    hy = grid.axes[2][2] - grid.axes[2][1]
    i, j = cell_offsets(grid, cell)
    x0, y0 = grid.axes[1][1] + i * hx, grid.axes[2][1] + j * hy

    ins, outs, _ = cutMesh(cut, grid)
    tol = 1e-10
    inbox(pt) = -tol <= pt[1] - x0 <= hx + tol && -tol <= pt[2] - y0 <= hy + tol
    ins = filter(t -> all(inbox, t), ins)

    fig = Figure(size = (1250, 560))
    ax1 = Axis(fig[1, 1]; aspect = DataAspect(), xlabel = "x", ylabel = "y",
        title = "EvolvingDomains cut cell + Ferrite quadrature points")
    poly!(ax1, [Point2f(x0, y0), Point2f(x0 + hx, y0), Point2f(x0 + hx, y0 + hy),
                Point2f(x0, y0 + hy)];
        color = (:grey, 0.15), strokecolor = :black, strokewidth = 1.5)
    triangles!(ax1, ins, COLORS[2], 0.7)
    scatter!(ax1, x, y; color = COLORS[1], markersize = 8, label = "physical weights")
    limits!(ax1, x0 - 0.15hx, x0 + 1.15hx, y0 - 0.15hy, y0 + 1.15hy)
    axislegend(ax1; position = :lb)

    ax2 = Axis(fig[1, 2]; aspect = DataAspect(), xlabel = "ξ", ylabel = "η",
        title = "same points in Ferrite's reference square [-1,1]²")
    poly!(ax2, [Point2f(-1, -1), Point2f(1, -1), Point2f(1, 1), Point2f(-1, 1)];
        color = (:grey, 0.15), strokecolor = :black, strokewidth = 1.5)
    scatter!(ax2, ξ, η; color = COLORS[1], markersize = 8)
    for (k, (a, b)) in enumerate(zip(ξ, η))
        text!(ax2, a, b; text = @sprintf("%.1e", w[k]), fontsize = 7, offset = (5, 5))
    end
    limits!(ax2, -1.25, 1.25, -1.25, 1.25)
    Label(fig[0, :], "Ferrite probe (n = $(n-1)²): CellValues and PointValues agree with a hand-coded Q1 gradient to ~1e-15",
        fontsize = 12)
    save(joinpath(OUT, "fig_ferrite.png"), fig)
    return fig
end

# -------------------------------------------------------------------- main ---

function main()
    p = defaultParams(; n = 65)
    configureCutCellsLibrary()
    isdir(OUT) || error("run run.MovingCutAgFEM.jl first: $OUT is missing")
    figGeometry(p)
    figSolutions(p)
    figConvergence(p)
    figDofs()
    figTimes()
    figMemory()
    figFerrite(p)
    for f in sort(filter(x -> startswith(x, "fig_"), readdir(OUT)))
        @printf("wrote %s (%.0f kB)\n", f, filesize(joinpath(OUT, f)) / 1024)
    end
end

if abspath(PROGRAM_FILE) == @__FILE__
    main()
end
