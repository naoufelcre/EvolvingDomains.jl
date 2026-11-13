# Visualisation of the one-phase planar Stefan melting run (standalone: reads the
# CSV files written by run.PlanarStefan.jl, no Gridap / EvolvingDomains needed).
using CairoMakie
using DelimitedFiles
using Printf

CairoMakie.activate!()

"Read a CSV with a header line into `Dict{Symbol,Vector{Float64}}`."
function read_csv(path::AbstractString)
    raw = readdlm(path, ',', skipstart = 1)
    header = strip.(split(readline(path), ','))
    cols = Dict{Symbol,Vector{Float64}}()
    for (k, name) in enumerate(header)
        v = vec(raw[:, k])
        cols[Symbol(name)] = Float64.(v)
    end
    return cols
end

"CairoMakie keeps no window open; drop the plots so the figure can be collected."
release!(fig) = (empty!(fig.scene.plots); nothing)

uniq_sorted(v) = sort!(unique(round.(Float64.(v), digits = 12)))

"Snapshots for a single time -> (xs, ys, T[i,j], phi[i,j])."
function snapshot_grid(snap, t)
    idx = findall(x -> isapprox(x, t; atol = 1e-12), snap[:t])
    isempty(idx) && error("no snapshot at t = $t")
    nx = round(Int, maximum(snap[:i][idx]))
    ny = round(Int, maximum(snap[:j][idx]))
    x = zeros(nx); y = zeros(ny)
    T = fill(NaN, nx, ny); P = fill(NaN, nx, ny)
    for k in idx
        i = round(Int, snap[:i][k]); j = round(Int, snap[:j][k])
        x[i] = snap[:x][k]; y[j] = snap[:y][k]
        T[i, j] = snap[:T][k]; P[i, j] = snap[:phi][k]
    end
    return x, y, T, P
end

function fig_heatmaps(snap, hist, times, path)
    fig = Figure(size = (1100, 900))
    hm = nothing
    for (n, t) in enumerate(times)
        x, y, T, P = snapshot_grid(snap, t)
        ax = Axis(fig[n > 2 ? 2 : 1, n > 2 ? n - 2 : n];
                  title = @sprintf("t = %.2f", t), aspect = DataAspect(),
                  xlabel = "x [-]", ylabel = "y [-]")
        solid = P .> 0
        Tp = copy(T); Tp[solid] .= NaN          # blank the solid phase
        hm = heatmap!(ax, x, y, Tp; colormap = :thermal, colorrange = (0, 1),
                      nan_color = (:gray35, 0.85))
        contour!(ax, x, y, P; levels = [0.0], color = :white, linewidth = 1.4)
        k = argmin(abs.(hist[:t] .- t))
        vlines!(ax, [hist[:s_exact][k]]; color = :black, linestyle = :dash,
                linewidth = 1.2, label = "exact interface")
        n == 1 && axislegend(ax; position = :rt, framevisible = false,
                             labelsize = 9)
    end
    Colorbar(fig[1:2, 3], hm; label = "T [-]", width = 14)
    Label(fig[0, 1:2], "temperature, interface (phi = 0, white) and exact interface (dashed); gray = solid",
          fontsize = 12)
    save(path, fig); release!(fig)
    return path
end

function fig_history(h, path)
    t = h[:t]
    fig = Figure(size = (1100, 800))
    ax1 = Axis(fig[1, 1]; ylabel = "s [L]", xlabel = "t [s]",
               title = "Interface position")
    lines!(ax1, t, h[:s_exact]; color = :black, linestyle = :dash,
           label = "s_exact (similarity)")
    lines!(ax1, t, h[:s_num]; color = :crimson, label = "s_num (cut FEM)")
    axislegend(ax1; position = :rb, framevisible = false)

    ax2 = Axis(fig[1, 2]; ylabel = "s_err [L]", xlabel = "t [s]",
               title = "Interface-position error")
    lines!(ax2, t, h[:s_err]; color = :navy)
    hlines!(ax2, [0.0]; color = :gray60, linestyle = :dot)

    ax3 = Axis(fig[2, 1]; ylabel = "q [K/L]", xlabel = "t [s]",
               title = "Interface heat flux")
    lines!(ax3, t, h[:q_exact]; color = :black, linestyle = :dash,
           label = "q_exact")
    lines!(ax3, t, h[:q_num]; color = :seagreen, label = "q_num")
    axislegend(ax3; position = :rb, framevisible = false)

    ax4 = Axis(fig[2, 2]; ylabel = "|r_stefan| [-]", xlabel = "t [s]",
               title = "Stefan-condition residual", yscale = log10)
    r = abs.(h[:r_stefan])
    r[r .== 0] .= NaN
    lines!(ax4, t, r; color = :darkorange)

    Label(fig[0, 1:2], "history.csv — n = $(length(t)) steps", fontsize = 12)
    save(path, fig); release!(fig)
    return path
end

function fig_profiles(prof, hist, times, path)
    fig = Figure(size = (900, 620))
    ax = Axis(fig[1, 1]; xlabel = "x [-]", ylabel = "T [-]", limits = ((-0.02, 1.02), (-0.03, 1.06)),
              title = "Temperature profiles at mid-height",
              subtitle = "markers = T_num, lines = T_exact, dash-dotted = interface (exact, num)")
    colors = cgrad(:viridis, length(times); categorical = true)
    for (n, t) in enumerate(times)
        idx = findall(x -> isapprox(x, t; atol = 1e-12), prof[:t])
        x, Te, Tn, phi = prof[:x][idx], prof[:T_exact][idx], prof[:T_num][idx], prof[:phi][idx]
        liquid = phi .< 0                     # exact solution is only meaningful in the liquid
        lines!(ax, x[liquid], Te[liquid]; color = colors[n], linewidth = 2,
               label = @sprintf("t = %.2f", t))
        scatter!(ax, x, Tn; color = colors[n], markersize = 5, marker = :circle)
        k = argmin(abs.(hist[:t] .- t))
        vlines!(ax, [hist[:s_exact][k], hist[:s_num][k]];
                color = (colors[n], 0.55), linestyle = :dashdot, linewidth = 1.1)
    end
    axislegend(ax; position = :rt, framevisible = false, labelsize = 9)
    save(path, fig); release!(fig)
    return path
end

"""
    visualizePlanarStefan(outdir; nheat = 4, nprof = 5)

Read `history.csv`, `profiles.csv`, `snapshots.csv` from `outdir` and write
`fig_heatmaps.png`, `fig_history.png`, `fig_profiles.png`. Returns the paths.
"""
function visualizePlanarStefan(outdir::AbstractString; nheat::Int = 4, nprof::Int = 5)
    hist = read_csv(joinpath(outdir, "history.csv"))
    prof = read_csv(joinpath(outdir, "profiles.csv"))
    snap = read_csv(joinpath(outdir, "snapshots.csv"))

    snakes = uniq_sorted(snap[:t])
    ptimes = uniq_sorted(prof[:t])
    heat_times = snakes[round.(Int, range(1, length(snakes); length = min(nheat, length(snakes))))]
    prof_times = ptimes[round.(Int, range(1, length(ptimes); length = min(nprof, length(ptimes))))]

    out = String[]
    push!(out, fig_heatmaps(snap, hist, heat_times, joinpath(outdir, "fig_heatmaps.png")))
    push!(out, fig_history(hist, joinpath(outdir, "fig_history.png")))
    push!(out, fig_profiles(prof, hist, prof_times, joinpath(outdir, "fig_profiles.png")))
    return out
end

if abspath(PROGRAM_FILE) == @__FILE__
    for p in visualizePlanarStefan(isempty(ARGS) ? joinpath(@__DIR__, "output") : ARGS[1])
        println("wrote ", p)
    end
end
