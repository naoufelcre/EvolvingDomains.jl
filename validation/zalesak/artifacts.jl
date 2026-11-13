# Figures and run records for the Zalesak case.
#
# Plain functions with explicit arguments: a variant builds the data, calls one figure
# function, saves it with `save_figure`, and writes its numbers with `report`.
using CairoMakie
using SHA: sha256
using Libdl: dlext
import EvolvingDomains: Field, Discretization, Grid, cut

_fmt(v) = string(round(v; sigdigits=5))

_refuse(paths) = isempty(filter(ispath, paths)) ||
    throw(ArgumentError("refusing to overwrite $paths; use a new output directory"))

"""
    native_cut_library()

Path to the native cut library, or `nothing`. Uses `EVOLVINGDOMAINS_CUTCELLS_LIBRARY` when
set, otherwise the local developer build from `scripts/build_cxx_bindings.sh`.
"""
function native_cut_library()
    override = get(ENV, "EVOLVINGDOMAINS_CUTCELLS_LIBRARY", "")
    isempty(override) || return isfile(override) ? override : nothing
    candidate = joinpath(tempdir(), "evolvingdomains-cutcells-install", "lib",
                         "libevolvingdomains_cut_cells." * dlext)
    isfile(candidate) ? candidate : nothing
end

"""
    save_figure_or_skip(build, dir, name)

Build and save a figure, or skip it when the 3D cut library is unavailable. The reported
numbers are the validation; the picture is evidence, and a long solve must not fail on it.
"""
function save_figure_or_skip(build, dir, name)
    library = native_cut_library()
    if library === nothing
        println("skipping the $name figure: no native cut library found. Build it with ",
                "scripts/build_cxx_bindings.sh or set EVOLVINGDOMAINS_CUTCELLS_LIBRARY")
        return nothing
    end
    ENV["EVOLVINGDOMAINS_CUTCELLS_LIBRARY"] = library
    save_figure(build(), dir, name)
end

"""Write PNG and PDF, verify both are real files, and confirm the overwrite guard works."""
function save_figure(fig, dir, name)
    paths = (png=joinpath(dir, "$name.png"), pdf=joinpath(dir, "$name.pdf"))
    _refuse(collect(paths))
    mkpath(dir)
    save(paths.png, fig)
    save(paths.pdf, fig)
    open(io -> read(io, 8), paths.png) == UInt8[0x89, 0x50, 0x4e, 0x47, 0x0d, 0x0a, 0x1a, 0x0a] ||
        error("not a PNG: $(paths.png)")
    open(io -> String(read(io, 5)), paths.pdf) == "%PDF-" || error("not a PDF: $(paths.pdf)")
    filesize(paths.png) > 1024 && filesize(paths.pdf) > 1024 ||
        error("empty figure: $paths")
    try
        _refuse(collect(paths))
        error("the overwrite guard accepted an existing file")
    catch e
        e isa ArgumentError || rethrow()
    end
    println("saved ", paths.png)
    paths
end

"""Write `run.txt`: what was run, the numbers it produced, and the source hashes."""
function report(dir, title, metrics, files)
    open(joinpath(dir, "run.txt"), "w") do io
        println(io, title)
        println(io, "Julia=", VERSION, "; threads=", Threads.nthreads())
        for (key, value) in pairs(metrics)
            println(io, key, "=", value)
        end
        println(io, "sampling: midpoint cells of the multilinear reconstruction, not a ",
                "conservative projection; interface distances are nearest edge samples, ",
                "not exact Hausdorff")
        println(io, "timing: solver loop only; setup, warmup, snapshots and plotting excluded")
        println(io, "sha256 (files are mutable in the checkout):")
        for file in files
            isfile(file) && println(io, basename(file), "=", bytes2hex(sha256(read(file))))
        end
    end
end

_metric_text(metrics) = join(["$k=$(_fmt(v))" for (k, v) in pairs(metrics)], "   |   ")

"""
    figure_2d(axis, times, snapshots, exact, metrics)

Four contour panels plus a final error heatmap and the numbers, for a 2D run.
"""
function figure_2d(axis, times, snapshots, exact, metrics)
    fig = Figure(size=(1200, 920))
    Label(fig[0, 1:4], "Field WENO5 / SSP-RK3"; fontsize=20, font=:bold)
    selected = round.(Int, range(1, length(times); length=5))
    for (panel, k) in enumerate(selected)
        row, col = divrem(panel - 1, 3) .+ (1, 1)
        ax = Axis(fig[row, col], aspect=DataAspect(), title="t = $(round(times[k]; sigdigits=4))",
                  xlabel="x", ylabel="y", xgridvisible=false, ygridvisible=false)
        limits!(ax, first(axis), last(axis), first(axis), last(axis))
        contour!(ax, axis, axis, snapshots[k]; levels=[0.0], color=:dodgerblue3, linewidth=2,
                 label="numerical")
        contour!(ax, axis, axis, exact[k]; levels=[0.0], color=:orangered, linewidth=1.5,
                 linestyle=:dash)
    end
    Legend(fig[1, 4], [LineElement(color=:dodgerblue3, linewidth=2),
                       LineElement(color=:orangered, linewidth=1.5, linestyle=:dash)],
           ["numerical", "analytic (sampled)"]; framevisible=false)
    error_ax = Axis(fig[2, 3], aspect=DataAspect(),
                    title="final |numerical − analytic (sampled)|",
                    subtitle="clamped-edge boundary artifacts are expected",
                    xlabel="x", ylabel="y", xgridvisible=false, ygridvisible=false)
    limits!(error_ax, first(axis), last(axis), first(axis), last(axis))
    heatmap!(error_ax, axis, axis, abs.(snapshots[end] .- exact[end]); colormap=:magma)
    Label(fig[3, 1:4], _metric_text(metrics); tellwidth=false, justification=:left, fontsize=14)
    fig
end

"""Three surfaces (initial cut, analytic final, numerical final), one slice, shape history."""
function figure_3d(axes, center, times, snapshots, exact, metrics, title)
    surfaces = [_surface(axes, u) for u in (first(snapshots), last(exact), last(snapshots))]
    fig = Figure(size=(1450, 1050))
    Label(fig[0, 1:3], title; fontsize=22, font=:bold)
    axis = axes[1]
    zoom = 0.35*(last(axis) - first(axis))
    for (k, name) in enumerate(("initial native cut", "analytic final (sampled cut)",
                                "final numerical"))
        ax = Axis3(fig[1, k]; title=name, xlabel="x", ylabel="y", zlabel="z",
                   aspect=:data, azimuth=π/5, elevation=π/6)
        mesh!(ax, surfaces[k].positions, surfaces[k].faces;
              color=k == 3 ? :dodgerblue3 : :darkorange)
        limits!(ax, center[1] .- zoom, center[1] .+ zoom, center[2] .- zoom,
                center[2] .+ zoom, center[3] .- zoom, center[3] .+ zoom)
    end
    mid = argmin(abs.(axis .- center[3]))
    slice_ax = Axis(fig[2, 1], title="final z ≈ $(round(axis[mid]; digits=3))",
                    xlabel="x", ylabel="y", aspect=DataAspect())
    contour!(slice_ax, axis, axis, @view last(snapshots)[:, :, mid]; levels=[0.0],
             color=:dodgerblue3, linewidth=2)
    contour!(slice_ax, axis, axis, @view last(exact)[:, :, mid]; levels=[0.0],
             color=:darkorange, linestyle=:dash)
    Label(fig[2, 2:3], _metric_text(metrics) *
          "\nnative tetrahedral volume: initial=$(_fmt(surfaces[1].volume)), " *
          "analytic=$(_fmt(surfaces[2].volume)), final=$(_fmt(surfaces[3].volume))",
          fontsize=13, tellwidth=false, justification=:left)
    fig
end

"""Initial/final surfaces from owned Cuts, with their native tetrahedral volume."""
function _surface(axes, u)
    geometry = cut(Field(Discretization(Dict(:init => Grid(axes))), Dict(:init => copy(u))),
                   :init)
    boundary = geometry.boundary
    positions = [Point3f(p...) for p in boundary.physical_points]
    faces = [CairoMakie.TriangleFace(
        Int.(boundary.point_ids[boundary.point_pointers[f]:boundary.point_pointers[f+1]-1])...)
        for f in eachindex(boundary.facet_to_bgcell)]
    (; positions, faces, volume=_tetra_volume(geometry))
end

function _tetra_volume(geometry)
    total = 0.0
    for s in eachindex(geometry.cell_to_bgcell)
        geometry.subcell_states[s] == -1 || continue
        ids = geometry.point_ids[geometry.point_pointers[s]:geometry.point_pointers[s+1]-1]
        p, a, b, c = geometry.physical_points[ids]
        a = a .- p; b = b .- p; c = c .- p
        total += abs(a[1]*(b[2]*c[3]-b[3]*c[2]) - a[2]*(b[1]*c[3]-b[3]*c[1]) +
                     a[3]*(b[1]*c[2]-b[2]*c[1]))/6
    end
    total
end

"""Overlay two solvers and the analytic contour on the common grid, at matching times."""
function figure_comparison_2d(axis, times, snapshots_a, snapshots_b, exact, metrics, title)
    fig = Figure(size=(1200, 920))
    Label(fig[0, 1:4], title; fontsize=20, font=:bold)
    selected = round.(Int, range(1, length(times); length=5))
    for (panel, k) in enumerate(selected)
        row, col = divrem(panel - 1, 3) .+ (1, 1)
        ax = Axis(fig[row, col], aspect=DataAspect(), title="t = $(round(times[k]; sigdigits=4))",
                  xlabel="x", ylabel="y", xgridvisible=false, ygridvisible=false)
        limits!(ax, first(axis), last(axis), first(axis), last(axis))
        contour!(ax, axis, axis, snapshots_a[k]; levels=[0.0], color=:dodgerblue3, linewidth=2)
        contour!(ax, axis, axis, snapshots_b[k]; levels=[0.0], color=:seagreen, linewidth=2)
        contour!(ax, axis, axis, exact[k]; levels=[0.0], color=:orangered, linewidth=1.5,
                 linestyle=:dash)
    end
    Legend(fig[1, 4], [LineElement(color=:dodgerblue3, linewidth=2),
                       LineElement(color=:seagreen, linewidth=2),
                       LineElement(color=:orangered, linewidth=1.5, linestyle=:dash)],
           ["Field WENO5", "Trixi DGSEM", "analytic"]; framevisible=false)
    Label(fig[2, 1:4], _metric_text(metrics); tellwidth=false, justification=:left, fontsize=13)
    fig
end

"""Two solvers' final surfaces against the analytic one."""
function figure_comparison_3d(axes, center, snapshots_a, snapshots_b, exact, metrics, title)
    surfaces = [_surface(axes, last(snapshots_a)), _surface(axes, last(snapshots_b)),
                _surface(axes, last(exact))]
    fig = Figure(size=(1450, 1050))
    Label(fig[0, 1:3], title; fontsize=22, font=:bold)
    axis = axes[1]
    zoom = 0.35*(last(axis) - first(axis))
    for (k, name) in enumerate(("final Field WENO5", "final Trixi DGSEM", "final analytic"))
        ax = Axis3(fig[1, k]; title=name, xlabel="x", ylabel="y", zlabel="z",
                   aspect=:data, azimuth=π/5, elevation=π/6)
        mesh!(ax, surfaces[k].positions, surfaces[k].faces;
              color=k == 3 ? :darkorange : :dodgerblue3)
        limits!(ax, center[1] .- zoom, center[1] .+ zoom, center[2] .- zoom,
                center[2] .+ zoom, center[3] .- zoom, center[3] .+ zoom)
    end
    Label(fig[2, 1:3], _metric_text(metrics) *
          "\nnative tetrahedral volume: ED=$(_fmt(surfaces[1].volume)), " *
          "Trixi=$(_fmt(surfaces[2].volume)), analytic=$(_fmt(surfaces[3].volume))",
          fontsize=13, tellwidth=false, justification=:left)
    fig
end
