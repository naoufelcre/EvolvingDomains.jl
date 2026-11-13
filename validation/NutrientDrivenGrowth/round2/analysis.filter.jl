#!/usr/bin/env julia
# Read-only analysis of the sub-resolution filter on the NutrientGrowthRefined experiment.
# Consumes output/<tag>/ (metrics.csv, frames.jls) and writes output/figures/filter_*.png.
using Serialization, Printf, Statistics
using CairoMakie
CairoMakie.activate!()

const OUT = joinpath(@__DIR__, "output")
const FIGD = joinpath(OUT, "figures")
mkpath(FIGD)

const ED = try
    @eval using EvolvingDomains
    EvolvingDomains
catch err
    error("load EvolvingDomains in the project environment: $err")
end

# ------------------------------------------------------------------ CSV helpers
function readcsv(tag)
    path = joinpath(OUT, tag, "metrics.csv")
    isfile(path) || return nothing
    lines = [l for l in readlines(path) if !isempty(strip(l))]
    hdr = Symbol.(split(lines[1], ","))
    rows = [split(l, ",") for l in lines[2:end]]
    return Dict(h => [parse(Float64, r[k]) for r in rows] for (k, h) in enumerate(hdr))
end

# ------------------------------------------------------------------ cut measure
function subcell_measure(snapshot, subcell)
    ids = Int(snapshot.point_pointers[subcell]):(Int(snapshot.point_pointers[subcell + 1]) - 1)
    a, b, c = Int.(snapshot.point_ids[ids])
    p, q, r = snapshot.physical_points[[a, b, c]]
    return abs((q[1] - p[1]) * (r[2] - p[2]) - (q[2] - p[2]) * (r[1] - p[1])) / 2
end

inside_area(snapshot) = sum(subcell_measure(snapshot, i) for i in eachindex(snapshot.subcell_states)
                            if snapshot.subcell_states[i] == ED.CUT_IN; init = 0.0)

function first_merged(metrics, column)
    (metrics === nothing || !haskey(metrics, column)) && return NaN
    v = metrics[column]
    i = findfirst(==(1.0), v)
    return i === nothing ? NaN : metrics[:t_end][i]
end

window_mean(x, t, lo, hi) = (m = (t .>= lo) .& (t .<= hi); isempty(x[m]) ? NaN : sum(x[m]) / sum(m))
window_max(x, t, lo, hi) = (m = (t .>= lo) .& (t .<= hi); isempty(x[m]) ? NaN : maximum(x[m]))

function summary(tag)
    m = readcsv(tag)
    m === nothing && return nothing
    t = m[:t_end]
    area_gain = m[:area_end][end] - m[:area_start][1]
    growth = sum(m[:integral_growth] .* m[:dt])                 # ∫ ∫_Ω g dt
    flux = sum(m[:interface_flux] .* m[:dt])
    filterΔ = haskey(m, :filter_area_change) ? sum(m[:filter_area_change]) : 0.0
    active = haskey(m, :filter_nodes_changed) ? count(>(0), m[:filter_nodes_changed]) : 0
    nodes = haskey(m, :filter_nodes_changed) ? sum(m[:filter_nodes_changed]) : 0.0
    eps_scale = 1.0e-3 * 1.5 * pi * 0.15^2
    debiased = haskey(m, :filter_area_change) ?
        [abs(m[:d_area_dt][i] - m[:filter_area_change][i] / m[:dt][i] - m[:integral_growth][i]) /
         max(abs(m[:integral_growth][i]), eps_scale) for i in eachindex(t)] :
        m[:balance_residual]
    return (; tag, ncells = m[:area_start][1] > 0 ? NaN : NaN,
        bal_debiased_mean = mean(debiased), bal_debiased_max = maximum(debiased),
        bal_debiased_pre = window_mean(debiased, t, 0.0, 0.34),
        bal_debiased_post = window_mean(debiased, t, 0.35, 0.50),
        debiased,
        merged = first_merged(m, :components_geometric),
        merged_cells = first_merged(m, :components_cells),
        area_final = m[:area_end][end], perimeter_final = m[:perimeter][end],
        area_gain, growth, flux, filterΔ, active, nodes,
        filter_nodes_max = haskey(m, :filter_nodes_changed) ? maximum(m[:filter_nodes_changed]) : 0.0,
        filterΔ_max = haskey(m, :filter_area_change) ? maximum(abs.(m[:filter_area_change])) : 0.0,
        bal_mean = mean(m[:balance_residual]), bal_max = maximum(m[:balance_residual]),
        bal_pre = window_mean(m[:balance_residual], t, 0.0, 0.34),
        bal_post = window_mean(m[:balance_residual], t, 0.35, 0.50),
        flux_bal_mean = mean(m[:flux_balance_residual]),
        ncut_mean = mean(m[:n_cut_cells]), ncut_final = m[:n_cut_cells][end],
        minfraction_min = minimum(m[:min_cut_fraction]),
        n_cut_min = minimum(m[:n_cut_cells]),
        t = t, area = m[:area_end], balance = m[:balance_residual],
        area_at = [m[:area_end][argmin(abs.(t .- s))] for s in (0.10, 0.20, 0.30, 0.40, 0.50)],
        filterΔ_series = haskey(m, :filter_area_change) ? m[:filter_area_change] : zeros(length(t)),
        ncut = m[:n_cut_cells])
end

const TAGS = ["n64_A", "n64_A_filt_c1", "n64_A_filt_c2", "n64_A_filt_o1",
              "n128_A", "n128_A_filt_c1", "n128_A_filt_o1"]

rows = Dict(tag => summary(tag) for tag in TAGS)
S = [rows[t] for t in TAGS if rows[t] !== nothing]

println("\n================ run summary ================")
@printf("%-16s %8s %8s %10s %10s %10s %8s %8s %8s %8s\n",
    "tag", "merg_g", "merg_c", "area_fin", "A@0.20", "A@0.50", "bal_mean", "bal_max", "bal_post", "ncut_end")
for r in S
    @printf("%-16s %8.3f %8.3f %10.6f %10.6f %10.6f %8.4f %8.4f %8.4f %8.0f\n",
        r.tag, r.merged, r.merged_cells, r.area_final, r.area_at[2], r.area_at[5],
        r.bal_mean, r.bal_max, r.bal_post, r.ncut_final)
end

println("\n--- growth-balance residual with the filter's own area injection removed ---")
@printf("%-16s %10s %10s %10s %10s %10s\n",
    "tag", "raw mean", "debiased", "deb pre", "deb post", "max deb")
for r in S
    @printf("%-16s %10.4f %10.4f %10.4f %10.4f %10.4f\n",
        r.tag, r.bal_mean, r.bal_debiased_mean, r.bal_debiased_pre, r.bal_debiased_post,
        r.bal_debiased_max)
end

println("\n--- area growth vs growth integral (cumulative) ---")
@printf("%-16s %10s %10s %10s %10s %10s %8s\n",
    "tag", "ΔA_total", "∫∫g", "ΔA-∫g", "(ΔA-∫g)/∫g", "Σ filterΔA", "% of ΔA")
for r in S
    @printf("%-16s %10.6f %10.6f %+10.6f %+10.2f%% %10.6f %7.2f%%\n",
        r.tag, r.area_gain, r.growth, r.area_gain - r.growth,
        100 * (r.area_gain - r.growth) / r.growth, r.filterΔ,
        100 * r.filterΔ / r.area_gain)
end

println("\n--- filter activity ---")
@printf("%-16s %8s %8s %8s %10s %10s %8s %8s %10s\n",
    "tag", "steps", "nodes", "nodes/st", "maxΔA/step", "ΣΔA", "ncut mean", "ncut min", "min_cut_frac")
for r in S
    @printf("%-16s %8d %8.0f %8.2f %10.2e %10.6f %8.1f %8.0f %10.2e\n",
        r.tag, r.active, r.nodes, r.active > 0 ? r.nodes / r.active : 0.0,
        r.filterΔ_max, r.filterΔ, r.ncut_mean, r.n_cut_min, r.minfraction_min)
end

# ------------------------------------------- offline probe on the baseline frames
# The frames carry the reinitialised φ, so |φ|/h at a node is its distance to the interface in
# cells. Filtering a copy of every saved frame therefore measures, per step, how far from the
# interface the operation reaches.
function probe(tag; mode = :closing, radius = 1)
    path = joinpath(OUT, tag, "frames.jls")
    isfile(path) || return nothing
    fr = deserialize(path)
    n = size(fr.frames[1].phi, 1)
    h = 1 / (n - 1)
    grid = ED.Grid((fr.axes, fr.axes))
    out = NamedTuple[]
    for f in fr.frames
        mask = f.phi .<= 0
        phi = copy(f.phi)
        field = ED.Field(ED.Discretization(Dict(:init => grid)), Dict(:init => phi))
        before = ED.cut(field, :init)
        stats = ED.filterSubResolutionFeatures!(field, :init; phase = :negative, mode, radius)
        after = ED.cut(field, :init)
        flipped = [I for I in CartesianIndices(phi) if (phi[I] <= 0) != mask[I]]
        depths = [abs(f.phi[I]) / h for I in flipped]           # distance to interface, cells
        push!(out, (; t = f.t, flipped = length(flipped),
            max_depth = isempty(depths) ? 0.0 : maximum(depths),
            mean_depth = isempty(depths) ? 0.0 : mean(depths),
            deep = count(>(2.0), depths),
            area_change = inside_area(after) - inside_area(before),
            ncut_change = length(ED.cut_cell_ids(after)) - length(ED.cut_cell_ids(before)),
            components_before = stats.n_components_before,
            components_after = stats.n_components_after,
            xlo = isempty(flipped) ? NaN : minimum(fr.axes[I[1]] for I in flipped),
            xhi = isempty(flipped) ? NaN : maximum(fr.axes[I[1]] for I in flipped),
            ylo = isempty(flipped) ? NaN : minimum(fr.axes[I[2]] for I in flipped),
            yhi = isempty(flipped) ? NaN : maximum(fr.axes[I[2]] for I in flipped)))
    end
    return out
end

println("\n--- offline probe: closing radius 1 applied to every saved baseline frame ---")
println("depths are |φ|/h of the flipped nodes, i.e. their distance from the interface in cells")
probes = Dict(tag => probe(tag) for tag in ("n64_A", "n128_A"))
for (tag, rows) in probes
    rows === nothing && continue
    h = tag == "n64_A" ? 1 / 64 : 1 / 128
    flipped = sum(r.flipped for r in rows)
    deep = sum(r.deep for r in rows)
    @printf("\n%s (%d frames, h=%.4g): %d nodes flipped in total, %d of them deeper than 2 cells\n",
        tag, length(rows), h, flipped, deep)
    @printf("%8s %8s %10s %10s %8s %12s %10s   %s\n",
        "t", "flipped", "mean depth", "max depth", ">2 cells", "Δarea", "Δcut cells", "modified x-range")
    for r in rows
        @printf("%8.3f %8d %10.2f %10.2f %8d %12.2e %10d   %s\n",
            r.t, r.flipped, r.mean_depth, r.max_depth, r.deep, r.area_change, r.ncut_change,
            r.flipped == 0 ? "-" : @sprintf("[%.3f,%.3f]x[%.3f,%.3f]", r.xlo, r.xhi, r.ylo, r.yhi))
    end
end

# ---------------------------------------------------------------------- figures
function frame_at(frames, t)
    i = argmin([abs(f.t - t) for f in frames])
    return frames[i]
end

window = (0.40, 0.62, 0.40, 0.60)      # xlo, xhi, ylo, yhi
tagsplot = [("n64_A", "baseline"), ("n64_A_filt_c1", "closing r=1")]
ts = [0.10, 0.35, 0.50]

fig = Figure(size = (1100, 1000))
hm = nothing
for (col, (tag, label)) in enumerate(tagsplot)
    fr = deserialize(joinpath(OUT, tag, "frames.jls"))
    for (row, t) in enumerate(ts)
        f = frame_at(fr.frames, t)
        ax = Axis(fig[row, col]; title = @sprintf("%s  t=%.2f  A=%.4f", label, f.t, f.area),
            xlabel = col == 1 && row == 3 ? "x" : "", ylabel = col == 1 ? "y" : "",
            aspect = DataAspect(), limits = (window[1], window[2], window[3], window[4]))
        global hm = heatmap!(ax, fr.axes, fr.axes, f.phi; colormap = :balance, colorrange = (-0.03, 0.03))
        contour!(ax, fr.axes, fr.axes, f.phi; levels = [0.0], color = :black, linewidth = 2)
        scatter!(ax, f.facet_mx, f.facet_my; markersize = 4, color = :white, strokewidth = 0.6,
            strokecolor = :black)
        vlines!(ax, [0.5], color = (:gray50, 0.4), linestyle = :dash)
        hlines!(ax, [0.5], color = (:gray50, 0.4), linestyle = :dash)
    end
end
Colorbar(fig[1:3, 3], hm; label = "φ", width = 14)
save(joinpath(FIGD, "filter_merger_zoom.png"), fig)
println("\nwrote ", joinpath(FIGD, "filter_merger_zoom.png"))

# --- tight magnification: is the void channel gone from the *reconstructed* geometry? ---
function window_cells(tag, t; xlo = 0.47, xhi = 0.53, ylo = 0.47, yhi = 0.53)
    fr = deserialize(joinpath(OUT, tag, "frames.jls"))
    f = frame_at(fr.frames, t)
    n = size(f.phi, 1)
    h = 1 / (n - 1)
    grid = ED.Grid((fr.axes, fr.axes))
    field = ED.Field(ED.Discretization(Dict(:init => grid)), Dict(:init => copy(f.phi)))
    snapshot = ED.cut(field, :init)
    nc = n - 1
    i0, i1 = max(1, floor(Int, xlo / h) + 1), min(nc, ceil(Int, xhi / h))
    j0, j1 = max(1, floor(Int, ylo / h) + 1), min(nc, ceil(Int, yhi / h))
    cells = [i + nc * (j - 1) for i in i0:i1 for j in j0:j1]
    states = snapshot.bg_states[cells]
    return (; t = f.t, t_neg = count(==(ED.CUT_IN), states), t_pos = count(==(ED.CUT_OUT), states),
        t_cut = count(==(ED.CUT_CUT), states), n = length(cells))
end

println("\n--- central window x,y ∈ [0.47,0.53]^2: reconstructed cell states ---")
println("(t_neg = entirely inside, t_pos = entirely in the void, t_cut = straddling)")
@printf("%6s  %-16s %6s %6s %6s\n", "t", "tag", "neg", "cut", "pos")
for t in (0.05, 0.10, 0.15, 0.20, 0.35, 0.50), tag in ("n64_A", "n64_A_filt_c1")
    w = window_cells(tag, t)
    @printf("%6.2f  %-16s %6d %6d %6d\n", w.t, tag, w.t_neg, w.t_cut, w.t_pos)
end

fig_micro = Figure(size = (900, 560))
hm2 = nothing
for (col, (tag, label)) in enumerate(tagsplot)
    fr = deserialize(joinpath(OUT, tag, "frames.jls"))
    for (row, t) in enumerate((0.10, 0.35))
        f = frame_at(fr.frames, t)
        ax = Axis(fig_micro[row, col];
            title = @sprintf("%s  t=%.2f  A=%.4f", label, f.t, f.area),
            aspect = DataAspect(), limits = (0.455, 0.545, 0.465, 0.535),
            xticks = 0.46:0.02:0.54, yticks = 0.47:0.02:0.53)
        global hm2 = heatmap!(ax, fr.axes, fr.axes, f.phi; colormap = :balance,
            colorrange = (-0.01, 0.01))
        contour!(ax, fr.axes, fr.axes, f.phi; levels = [0.0], color = :black, linewidth = 2.5)
        scatter!(ax, f.facet_mx, f.facet_my; markersize = 7, color = :white,
            strokewidth = 0.8, strokecolor = :black)
    end
end
Colorbar(fig_micro[1:2, 3], hm2; label = "φ", width = 14)
save(joinpath(FIGD, "filter_merger_micro.png"), fig_micro)
println("wrote ", joinpath(FIGD, "filter_merger_micro.png"))

fig2 = Figure(size = (1100, 750))
ax1 = Axis(fig2[1, 1]; xlabel = "t", ylabel = "area", title = "area evolution")
ax2 = Axis(fig2[1, 2]; xlabel = "t", ylabel = "balance residual", title = "growth balance |dA/dt − ∫g|/max(|∫g|,ε)")
ax3 = Axis(fig2[2, 1]; xlabel = "t", ylabel = "filter Δarea per step", title = "area injected by filtering")
ax4 = Axis(fig2[2, 2]; xlabel = "t", ylabel = "straddling cells (CUT)", title = "reconstructed cut complexity")
clr = Dict("n64_A" => :dodgerblue, "n64_A_filt_c1" => :orangered, "n64_A_filt_c2" => :gold, "n64_A_filt_o1" => :seagreen,
           "n128_A" => :navy, "n128_A_filt_c1" => :darkred, "n128_A_filt_o1" => :darkgreen)
for r in S
    lw = startswith(r.tag, "n128") ? 1.4 : 2.0
    lines!(ax1, r.t, r.area; color = clr[r.tag], linewidth = lw, label = r.tag)
    lines!(ax2, r.t, r.balance; color = clr[r.tag], linewidth = lw, label = r.tag)
    startswith(r.tag, "n64") && lines!(ax3, r.t, r.filterΔ_series; color = clr[r.tag], linewidth = lw,
        label = r.tag)
    lines!(ax4, r.t, r.ncut; color = clr[r.tag], linewidth = lw, label = r.tag)
end
axislegend(ax1; position = :lt, framevisible = false)
axislegend(ax2; position = :lt, framevisible = false)
save(joinpath(FIGD, "filter_histories.png"), fig2)
println("wrote ", joinpath(FIGD, "filter_histories.png"))
