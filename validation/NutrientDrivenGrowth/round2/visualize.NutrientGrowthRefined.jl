#!/usr/bin/env julia
# Plots for NutrientGrowthRefined. Read-only on the data produced by run.NutrientGrowthRefined.jl.
# Usage: julia --project=... visualize.NutrientGrowthRefined.jl [tag ...]
using CairoMakie, Serialization, Printf, Statistics
CairoMakie.activate!()

const OUT = joinpath(@__DIR__, "output")
const FIGD = joinpath(OUT, "figures")
const CLR = Dict("A" => :dodgerblue, "B" => :orangered, "" => :gray40)
warn(msg) = println("WARNING: ", msg)

# ---------------------------------------------------------------- input helpers
function readcsv(path)
    isfile(path) || return nothing
    lines = [strip(l) for l in readlines(path) if !isempty(strip(l))]
    length(lines) < 2 && return nothing
    hdr = Symbol.(split(lines[1], ','))
    rows = [split(l, ',') for l in lines[2:end]]
    return Dict(h => [r[k] for r in rows] for (k, h) in enumerate(hdr) if all(length(r) >= k for r in rows))
end
function num(c, name)
    key = Symbol(name)
    (c !== nothing && haskey(c, key)) || (warn("column :$key absent"); return nothing)
    return [something(tryparse(Float64, strip(s)), NaN) for s in c[key]]
end
strcol(c, name) = (key = Symbol(name); haskey(c, key) ? c[key] : (warn("column :$key absent"); nothing))
maybelegend(ax; kw...) = try axislegend(ax; kw...) catch; warn("no labelled plots, legend skipped") end
pickidx(n, k) = isempty(n) ? Int[] : unique(round.(Int, range(1, length(n), length = min(k, length(n)))))
selectframes(fr; k = 3) = (i = findall(f -> f.synchronized, fr.frames);
                           isempty(i) && (i = collect(eachindex(fr.frames))); fr.frames[pickidx(i, k)])
loadframes(dir) = isfile(joinpath(dir, "frames.jls")) ? deserialize(joinpath(dir, "frames.jls")) :
                  (warn("no frames.jls in $dir"); nothing)
function taginfo(tag)
    dir = joinpath(OUT, tag)
    fr = loadframes(dir)
    m = match(r"n(\d+)_([AB])", tag)
    nc = m === nothing ? (fr === nothing ? 0 : fr.ncells) : parse(Int, m.captures[1])
    meth = m === nothing ? (fr === nothing ? "" : (occursin("centroid", lowercase(fr.method)) ? "A" : "B")) : m.captures[2]
    return (; tag, dir, ncells = nc, method = meth, fr, metrics = readcsv(joinpath(dir, "metrics.csv")),
              profile = readcsv(joinpath(dir, "profile.csv")))
end
# one representative tag per method letter: the 64² primary grid when present, else the
# first grid that has it (REPORT.md designates 64² as the primary controlled comparison).
const PRIMARY_N = 64
function reps(infos)
    d = Dict{String,Any}()
    isempty(infos) && return d
    for ti in infos; ti.ncells == PRIMARY_N && get!(d, ti.method, ti); end
    n0 = first(infos).ncells
    for ti in infos; ti.ncells == n0 && get!(d, ti.method, ti); end
    for ti in infos; get!(d, ti.method, ti); end
    return d
end
# the montage/field figures use the 64² interface-aware run when it exists (the baseline of
# the report), then 64² centroid, then whatever has frames
function primary_index(infos)
    for want in ((PRIMARY_N, "B"), (PRIMARY_N, "A"))
        i = findfirst(ti -> ti.fr !== nothing && !isempty(ti.fr.frames) &&
                           ti.ncells == want[1] && ti.method == want[2], infos)
        i === nothing || return i
    end
    return findfirst(ti -> ti.fr !== nothing && !isempty(ti.fr.frames), infos)
end
function clims(mats)
    v = [x for M in mats for x in M if isfinite(x)]
    isempty(v) && return (0.0, 1.0)
    lo, hi = extrema(v)
    return lo == hi ? (lo, lo + 1) : (lo, hi)
end
function domain_limits(fr, frames; pad = 2)
    nx = length(fr.axes)
    i1, i2, j1, j2 = nx, 1, nx, 1
    for f in frames
        m = f.phi .< 0
        any(m) || continue
        ii = findall(m)
        i1 = min(i1, minimum(I -> I[1], ii)); i2 = max(i2, maximum(I -> I[1], ii))
        j1 = min(j1, minimum(I -> I[2], ii)); j2 = max(j2, maximum(I -> I[2], ii))
    end
    i1 > i2 && return ((0.0, 1.0), (0.0, 1.0))
    return ((fr.axes[max(i1 - pad, 1)], fr.axes[min(i2 + pad, nx)]),
            (fr.axes[max(j1 - pad, 1)], fr.axes[min(j2 + pad, nx)]))
end
# The interface is drawn from the cut-facet midpoints, not from a φ = 0 contour: the two-disc
# initial level set is a min() of two distances, so it has a crease on the medial line where the
# *bilinear* interpolant of φ takes a spurious zero crossing in the void. That crossing produces
# no cut facet and never enters the computation, but it does show up in a contour plot.
function draw_domain!(ax, fr, f)
    heatmap!(ax, fr.axes, fr.axes, ifelse.(f.phi .< 0, 1.0, NaN); colormap = cgrad([:steelblue, :steelblue]),
             colorrange = (0, 1))
    scatter!(ax, f.facet_mx, f.facet_my; color = :black, markersize = 3)
end
files = String[]
function emit(name, fig)
    path = joinpath(FIGD, name); save(path, fig); push!(files, path)
    @printf("  wrote %-28s %8.1f kB\n", name, filesize(path) / 1024)
end
# ---------------------------------------------------------------------- figures
function fig_domain_evolution(ti)
    ti.fr === nothing && return warn("tag $(ti.tag): no frames.jls")
    fr = ti.fr
    sel = fr.frames[pickidx(eachindex(fr.frames), 8)]
    isempty(sel) && return warn("tag $(ti.tag): no frames")
    xl, yl = domain_limits(fr, sel)
    nc = min(4, length(sel)); nr = cld(length(sel), nc)
    fig = Figure(size = (330 * nc, 250 * nr))
    for (k, f) in enumerate(sel)
        c = f.components
        ax = Axis(fig[cld(k, nc), mod1(k, nc)]; xlabel = "x [–]", ylabel = "y [–]", aspect = DataAspect(),
                  limits = (xl, yl), titlesize = 10, title = @sprintf("t = %.4f  |Ω| = %.4f", f.t, f.area),
                  subtitle = "components n/c/g = $(c.nodes)/$(c.cells)/$(c.geometric)")
        draw_domain!(ax, fr, f)
    end
    Label(fig[0, :], "material region Ω(t), interface φ = 0 — tag $(ti.tag) (n = $(fr.ncells), method $(ti.method))",
          tellwidth = false)
    emit("fig_domain_evolution.png", fig)
end

function fig_fields(ti)
    ti.fr === nothing && return warn("tag $(ti.tag): no frames.jls")
    fr = ti.fr; sel = selectframes(fr; k = 3)
    isempty(sel) && return warn("tag $(ti.tag): no frames")
    xl, yl = domain_limits(fr, sel)
    fig = Figure(size = (330 * length(sel) + 90, 560))
    for (row, key, title, lbl, cmap) in ((1, :c, "nutrient c", "c [–]", :viridis), (2, :p, "pressure p", "p [–]", :balance))
        lim = clims([getproperty(f, key) for f in sel])
        hm = nothing
        for (k, f) in enumerate(sel)
            ax = Axis(fig[row, k]; xlabel = "x [–]", ylabel = "y [–]", aspect = DataAspect(), limits = (xl, yl),
                      titlesize = 10, title = @sprintf("%s, t = %.4f", title, f.t))
            hm = heatmap!(ax, fr.axes, fr.axes, getproperty(f, key); colormap = cmap, colorrange = lim)
            contour!(ax, fr.axes, fr.axes, f.phi; levels = [0.0], color = :black, linewidth = 1.2)
        end
        hm === nothing || Colorbar(fig[row, length(sel) + 1], hm; label = lbl)
    end
    Label(fig[0, :], "Fields on Ω(t) — tag $(ti.tag)", tellwidth = false)
    emit("fig_fields.png", fig)
end

function fig_velocity(ti)
    ti.fr === nothing && return warn("tag $(ti.tag): no frames.jls")
    fr = ti.fr; f = fr.frames[end]
    h = fr.axes[2] - fr.axes[1]
    xl, yl = domain_limits(fr, [f])
    nc = size(f.cell_vx, 1)
    xc = fr.axes[1:nc] .+ h / 2
    spd = hypot.(f.cell_vx, f.cell_vy)
    vmax = maximum(filter(isfinite, spd)); sc = iszero(vmax) ? 0.0 : 2.5h / vmax
    fig = Figure(size = (1160, 430))
    ax = Axis(fig[1, 1]; xlabel = "x [–]", ylabel = "y [–]", aspect = DataAspect(), limits = (xl, yl), titlesize = 11,
              title = @sprintf("cell speed |v| and (vx, vy) every 4 nodes, t = %.4f", f.t))
    hm = heatmap!(ax, xc, xc, spd; colormap = :magma, colorrange = (0, vmax))
    Colorbar(fig[1, 2], hm; label = "|v| [–]")
    contour!(ax, fr.axes, fr.axes, f.phi; levels = [0.0], color = :black, linewidth = 1.3)
    pts = Point2f[]; dirs = Vec2f[]
    for i in 1:4:length(fr.axes), j in 1:4:length(fr.axes)
        (isfinite(f.vx[i, j]) && isfinite(f.vy[i, j]) && f.phi[i, j] < 2.5h) || continue
        push!(pts, Point2f(fr.axes[i], fr.axes[j])); push!(dirs, Vec2f(f.vx[i, j] * sc, f.vy[i, j] * sc))
    end
    arrows2d!(ax, pts, dirs; color = (:white, 0.9), tiplength = 6, shaftwidth = 1.2)
    ax2 = Axis(fig[1, 3]; xlabel = "x [–]", ylabel = "y [–]", aspect = DataAspect(), limits = (xl, yl),
               title = "interface normal velocity Vₙ", titlesize = 11)
    contour!(ax2, fr.axes, fr.axes, f.phi; levels = [0.0], color = (:gray60, 0.7), linewidth = 1.0)
    sc2 = scatter!(ax2, f.facet_mx, f.facet_my; color = f.facet_vn, colormap = :viridis, markersize = 9)
    Colorbar(fig[1, 4], sc2; label = "Vₙ [–]")
    Label(fig[0, :], "Velocity field — tag $(ti.tag), final frame", tellwidth = false)
    emit("fig_velocity.png", fig)
end

function fig_marker(ti)
    ti.fr === nothing && return warn("tag $(ti.tag): no frames.jls")
    fr = ti.fr; sel = selectframes(fr; k = 3)
    isempty(sel) && return warn("tag $(ti.tag): no frames")
    xl, yl = domain_limits(fr, sel)
    lim = clims([f.q for f in sel])
    fig = Figure(size = (330 * length(sel) + 90, 280))
    hm = nothing
    for (k, f) in enumerate(sel)
        ax = Axis(fig[1, k]; xlabel = "x [–]", ylabel = "y [–]", aspect = DataAspect(), limits = (xl, yl),
                  titlesize = 10, title = @sprintf("t = %.4f", f.t))
        hm = heatmap!(ax, fr.axes, fr.axes, f.q; colormap = :plasma, colorrange = lim)
        contour!(ax, fr.axes, fr.axes, f.phi; levels = [0.0], color = :black, linewidth = 1.3)
    end
    hm === nothing || Colorbar(fig[1, length(sel) + 1], hm; label = "marker q [–]")
    Label(fig[0, :], "Passive marker q with interface φ = 0 — tag $(ti.tag)", tellwidth = false)
    emit("fig_marker.png", fig)
end

function fig_balance(R)
    isempty(R) && return warn("fig_balance: no tags")
    fig = Figure(size = (1200, 450))
    ax = Axis(fig[1, 1]; xlabel = "t [–]", ylabel = "area rate / integral [–]", titlesize = 11,
              title = "Area balance terms")
    ax2 = Axis(fig[1, 2]; xlabel = "t [–]", ylabel = "relative residual [–]", yscale = log10, titlesize = 11,
               title = "Balance and closure residuals")
    for meth in sort(collect(keys(R)))
        ti = R[meth]; m = ti.metrics; m === nothing && (warn("tag $(ti.tag): no metrics.csv"); continue)
        t = num(m, "t_end"); (t === nothing || isempty(t)) && continue
        col = CLR[meth]
        for (name, lbl) in (("d_area_dt", "d|Ω|/dt"), ("integral_growth", "∫ growth"), ("interface_flux", "flux"))
            y = num(m, name); y === nothing || lines!(ax, t, y; color = col, linewidth = 2, label = "$meth: $lbl")
        end
        y = num(m, "balance_residual")
        y === nothing || lines!(ax2, t, abs.(y); color = col, linewidth = 2, label = "$meth: |area residual|")
        ifl = num(m, "interface_flux"); ig = num(m, "integral_growth")
        if ifl !== nothing && ig !== nothing
            lines!(ax2, t, abs.(ifl .- ig) ./ ifelse.(abs.(ig) .== 0, NaN, abs.(ig)); color = col, linewidth = 2,
                   linestyle = :dash, label = "$meth: |flux−growth|/|growth|")
        end
    end
    maybelegend(ax; position = :rb, fontsize = 9); maybelegend(ax2; position = :lb, fontsize = 9)
    emit("fig_balance.png", fig)
end

function fig_area(infos)
    fig = Figure(size = (900, 450))
    ax = Axis(fig[1, 1]; xlabel = "t [–]", ylabel = "area [–]", title = "Physical area versus time", titlesize = 11)
    got = false
    for ti in infos
        ti.metrics === nothing && continue
        t = num(ti.metrics, "t_end"); a = num(ti.metrics, "area_end")
        (t === nothing || a === nothing) && continue
        col = CLR[ti.method]
        lines!(ax, t, a; color = col, linewidth = 1.8, label = "$(ti.tag) (ED)")
        af = num(ti.metrics, "area_end_fem")
        af === nothing || lines!(ax, t, af; color = col, linewidth = 1.4, linestyle = :dash, label = "$(ti.tag) (FE)")
        got = true
    end
    got || return warn("fig_area: no area data")
    maybelegend(ax; position = :rb, fontsize = 8)
    emit("fig_area.png", fig)
end

function fig_topology(R)
    isempty(R) && return warn("fig_topology: no tags")
    fig = Figure(size = (1150, 380))
    for (k, name) in enumerate(("components_nodes", "components_cells", "components_geometric"))
        ax = Axis(fig[1, k]; xlabel = "t [–]", ylabel = "count [–]", titlesize = 11,
                  title = replace(name, "components_" => ""))
        for meth in sort(collect(keys(R)))
            ti = R[meth]; ti.metrics === nothing && continue
            t = num(ti.metrics, "t_end"); y = num(ti.metrics, name)
            (t === nothing || y === nothing) && continue
            lines!(ax, t, y; color = CLR[meth], linewidth = 2, label = "method $meth")
            scatter!(ax, t, y; color = CLR[meth], marker = :circle, markersize = 6)
        end
        k == 1 && maybelegend(ax; position = :lb)
    end
    emit("fig_topology.png", fig)
end

# a single point on a log axis otherwise gives absurd auto-limits, so pad instead of drawing a slope line
function refline!(ax, pts; slope = 1.0, lbl = "∝ h")
    isempty(pts) && return
    hs = first.(pts); ys = last.(pts)
    (all(>(0), hs) && all(>(0), ys)) || return
    lo, hi = extrema(hs)
    if lo == hi
        xlims!(ax, lo / 3, lo * 3); ylims!(ax, first(ys) / 3, first(ys) * 3); return
    end
    h0 = exp(mean(log.(hs))); y0 = exp(mean(log.(ys)))
    lines!(ax, [lo, hi], y0 .* ([lo, hi] ./ h0) .^ slope; color = :gray45, linestyle = :dash, linewidth = 1.6, label = lbl)
end

function fig_convergence(infos)
    fig = Figure(size = (1150, 450))
    ax = Axis(fig[1, 1]; xlabel = "h = 1/n", ylabel = "mean |balance residual| [–]", xscale = log10, yscale = log10,
              titlesize = 11, title = "Area-balance residual versus h")
    allpts = Point2f[]
    for meth in sort(unique(ti.method for ti in infos))
        hs = Float64[]; ys = Float64[]
        for ti in infos
            (ti.method == meth && ti.metrics !== nothing && ti.ncells > 0) || continue
            y = num(ti.metrics, "balance_residual"); y === nothing && continue
            v = mean(abs, filter(isfinite, y))
            (isfinite(v) && v > 0) || continue
            push!(hs, 1.0 / ti.ncells); push!(ys, v)
        end
        isempty(hs) && continue
        o = sortperm(hs)
        lines!(ax, hs[o], ys[o]; color = CLR[meth], linewidth = 2, label = "method $meth")
        scatter!(ax, hs[o], ys[o]; color = CLR[meth], marker = :circle, markersize = 8)
        append!(allpts, Point2f.(hs, ys))
    end
    refline!(ax, allpts); maybelegend(ax; position = :lb, fontsize = 9)

    ax2 = Axis(fig[1, 2]; xlabel = "h = 1/n", ylabel = "recovery error [–]", xscale = log10, yscale = log10,
               titlesize = 11, title = "Verification recovery errors versus h")
    v = readcsv(joinpath(OUT, "verification.csv"))
    v === nothing && warn("no verification.csv")
    nn = v === nothing ? nothing : num(v, "verification_ncells")
    if nn !== nothing
        o = sortperm(nn); h = 1.0 ./ nn[o]
        vpts = Point2f[]
        for (name, lbl, ls) in (("verification_recovery_A_mean_abs_error", "A: mean |error|", :solid),
                                ("verification_recovery_B_mean_abs_error", "B: mean |error|", :solid),
                                ("verification_recovery_A_flux_relative_error", "A: flux rel. error", :dash),
                                ("verification_recovery_B_flux_relative_error", "B: flux rel. error", :dash))
            y = num(v, name); y === nothing && continue
            y = y[o]; col = CLR[startswith(lbl, "A") ? "A" : "B"]
            lines!(ax2, h, y; color = col, linewidth = 1.8, linestyle = ls, label = lbl)
            scatter!(ax2, h, y; color = col, marker = ls == :solid ? :circle : :utriangle, markersize = 6)
            append!(vpts, Point2f.(h, y))
        end
        refline!(ax2, vpts; slope = 2.0, lbl = "∝ h²")
    end
    maybelegend(ax2; position = :lb, fontsize = 8)
    Label(fig[0, :], "Convergence: balance residual and recovery error", tellwidth = false)
    emit("fig_convergence.png", fig)
end

function fig_performance(infos)
    agg = Dict{String,Tuple{Dict{String,Tuple{Float64,Float64}},String}}()  # tag -> (stage -> (s, bytes), method)
    for ti in infos
        ti.profile === nothing && (warn("tag $(ti.tag): no profile.csv"); continue)
        st = strcol(ti.profile, "stage"); se = num(ti.profile, "seconds"); by = num(ti.profile, "bytes")
        sp = num(ti.profile, "step")
        (st === nothing || se === nothing || by === nothing || sp === nothing) && continue
        keep = sp .>= 3
        any(keep) || (keep .= true; warn("tag $(ti.tag): no profile steps >= 3, using all steps"))
        d = Dict{String,Tuple{Float64,Float64}}()
        for s in unique(st)
            idx = [i for i in eachindex(st) if st[i] == s && keep[i]]
            isempty(idx) || (d[s] = (mean(se[idx]), mean(by[idx])))
        end
        agg[ti.tag] = (d, ti.method)
    end
    isempty(agg) && return warn("fig_performance: no profile data")
    sv(d, s, k) = haskey(d, s) ? d[s][k] : NaN
    stages = unique(reduce(vcat, [collect(keys(x[1])) for x in values(agg)]))
    sort!(stages; by = s -> -maximum(sv(x[1], s, 1) for x in values(agg)))
    tgs = sort(collect(keys(agg))); K = length(tgs); bw = 0.8 / K
    fig = Figure(size = (1250, max(430, 24 * length(stages) + 150)))
    axs = [Axis(fig[1, 1]; xlabel = "mean seconds per step [s]", yticklabelsize = 8, titlesize = 11,
                yticks = (collect(1:length(stages)), stages), title = "Time per stage"),
           Axis(fig[1, 2]; xlabel = "mean bytes per step [B]", yticklabelsize = 8, titlesize = 11, xscale = log10,
                yticks = (collect(1:length(stages)), fill("", length(stages))), title = "Allocations per stage")]
    lo = Float64[]
    for (j, tag) in enumerate(tgs)
        d, meth = agg[tag]; ys = collect(1:length(stages)) .+ (j - (K + 1) / 2) * bw
        ts = [sv(d, s, 1) for s in stages]
        barplot!(axs[1], ys, [isfinite(x) ? max(x, 0.0) : 0.0 for x in ts]; direction = :x, width = bw,
                 color = CLR[meth], label = tag)
        by = [sv(d, s, 2) for s in stages]; ok = isfinite.(by) .& (by .> 0)
        if any(ok)
            barplot!(axs[2], ys[ok], by[ok]; direction = :x, width = bw, color = CLR[meth], label = tag)
            append!(lo, by[ok])
        else
            warn("tag $tag: no positive byte counts for the log-scale panel")
        end
    end
    !isempty(lo) && xlims!(axs[2], minimum(lo) / 3, maximum(lo) * 3)
    maybelegend(axs[1]; position = :rb, fontsize = 9)
    Label(fig[0, :], "Runtime and allocation profile (steps ≥ 3, mean per stage)", tellwidth = false)
    emit("fig_performance.png", fig)
end

# ------------------------------------------------------------------------ driver
mkpath(FIGD)
tags = isempty(ARGS) ? sort([d for d in readdir(OUT) if isdir(joinpath(OUT, d)) &&
                             (isfile(joinpath(OUT, d, "metrics.csv")) || isfile(joinpath(OUT, d, "frames.jls")))]) :
       copy(ARGS)
isempty(tags) && warn("no tag directories with metrics.csv or frames.jls under $OUT")
infos = [taginfo(t) for t in tags if isdir(joinpath(OUT, t))]
for ti in infos
    (ti.fr === nothing && ti.metrics === nothing) && warn("tag $(ti.tag) has neither frames.jls nor metrics.csv")
end
primary = primary_index(infos)
R = reps(infos)
println("tags: ", isempty(infos) ? "(none)" : join([ti.tag for ti in infos], ", "),
        " | methods: ", isempty(R) ? "(none)" : join(sort(collect(keys(R))), ","))

ptag(fn) = () -> (primary === nothing ? warn("no primary tag: skipping $(nameof(fn))") : fn(infos[primary]))
for (name, fn) in (("fig_domain_evolution", ptag(fig_domain_evolution)), ("fig_fields", ptag(fig_fields)),
                   ("fig_velocity", ptag(fig_velocity)), ("fig_marker", ptag(fig_marker)),
                   ("fig_balance", () -> fig_balance(R)), ("fig_area", () -> fig_area(infos)),
                   ("fig_topology", () -> fig_topology(R)), ("fig_convergence", () -> fig_convergence(infos)),
                   ("fig_performance", () -> fig_performance(infos)))
    try
        fn()
    catch e
        warn("$name failed: $(sprint(showerror, e))")
    end
end

println("\n", length(files), " file(s) written to $FIGD:")
foreach(p -> @printf("  %-34s %8.1f kB\n", basename(p), filesize(p) / 1024), files)
