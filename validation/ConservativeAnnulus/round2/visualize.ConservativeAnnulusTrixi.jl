# ConservativeAnnulusTrixi — figures. Reads output/results.jls; solves nothing.
#
# Run:  julia validation/ConservativeAnnulus/round2/visualize.ConservativeAnnulusTrixi.jl

using Pkg
Pkg.activate(joinpath(@__DIR__, "..", "..", "env"); io=devnull)
using CairoMakie
using Serialization
using Printf

const OUTDIR = joinpath(@__DIR__, "output")
const METHODS = ("WENO5", "CIP", "Map", "Trixi DG3")
const COLORS = (:steelblue, :seagreen, :firebrick, :darkorange)

data = deserialize(joinpath(OUTDIR, "results.jls"))
exact_area, exact_mass = data.exact_area, data.exact_mass
cases = data.cases

byname = Dict(c.label => c for c in cases)
fine = byname[:ncells128]
coarse = byname[:ncells64]
finec = byname[:ncells128_const]
coarsec = byname[:ncells64_const]

rhoExact(x, y, t) = 1 + 0.5 * sum(exp(-((x - 0.41 - 0.24t + mx)^2 + (y - 0.48 - 0.10t + my)^2) / 0.055^2)
                                  for mx in -1:1, my in -1:1)
phiExact(x, y, t) = (cx = x - 0.35 - 0.24t; cy = y - 0.45 - 0.10t; r = hypot(cx, cy);
                     max(0.12 - r, r - 0.24))

snapAt(c, t) = c.snaps[argmin([abs(s.t - t) for s in c.snaps])]

# ---------------------------------------------------------------------------------------
# 1. Material-domain density, every method, three times.
# ---------------------------------------------------------------------------------------
function figure1()
    c = fine
    times = (0.0, 0.5, 1.0)
    fig = Figure(; size=(1150, 1300))
    for (row, ti) in enumerate(times)
        s = snapAt(c, ti)
        n = c.n
        ax_axis = range(0.0, 1.0; length=n)
        for (col, key) in enumerate((:weno, :cip, :map, :trixi))
            rho = getproperty(s, key)
            masked = [s.phi[i, j] < 0 ? rho[i, j] : NaN for i in 1:n, j in 1:n]
            ax = Axis(fig[row, col]; aspect=DataAspect(),
                      title=col == 1 ? @sprintf("%s  t = %.1f", METHODS[1], ti) :
                                           @sprintf("%s  t = %.1f", METHODS[col], ti),
                      titlesize=11, xlabelsize=9, ylabelsize=9)
            hm = heatmap!(ax, ax_axis, ax_axis, masked; colorrange=(0.72, 1.26), colormap=:viridis)
            contour!(ax, ax_axis, ax_axis, s.phi; levels=[0.0], color=:white, linewidth=1.6)
            if col == 4
                Colorbar(fig[row, end + 1], hm; label="ρ")
            end
            hidespines!(ax); hidedecorations!(ax)
        end
    end
    Label(fig[0, :], "Density on the numerically evolved material domain (white: φ = 0 from the level-set Field)",
          fontsize=13)
    save(joinpath(OUTDIR, "1_density_maps.png"), fig)
end

# ---------------------------------------------------------------------------------------
# 2. Density error at t = 1 on the same nodes.
# ---------------------------------------------------------------------------------------
function figure2()
    c = fine
    s = snapAt(c, 1.0)
    n = c.n
    ax_axis = range(0.0, 1.0; length=n)
    f = Figure(; size=(1150, 420))
    for (col, key) in enumerate((:weno, :cip, :map, :trixi))
        rho = getproperty(s, key)
        err = [s.phi[i, j] < 0 ? max(abs(rho[i, j] - rhoExact(ax_axis[i], ax_axis[j], 1.0)), 1e-7) : NaN
               for i in 1:n, j in 1:n]
        ax = Axis(f[1, col]; aspect=DataAspect(), title=METHODS[col], titlesize=12)
        hm = heatmap!(ax, ax_axis, ax_axis, err; colorrange=(1e-7, 1e-1), colorscale=log10,
                      colormap=:magma, lowclip=:white)
        contour!(ax, ax_axis, ax_axis, s.phi; levels=[0.0], color=:cyan, linewidth=1.2)
        Colorbar(f[2, col], hm; vertical=false, label="|ρ - ρ_exact|")
        hidespines!(ax); hidedecorations!(ax)
    end
    Label(f[0, :], "|Density error| at t = 1, 128² cells / 32² elements (log colour scale, 1e-7 to 1e-1)", fontsize=13)
    save(joinpath(OUTDIR, "2_error_maps.png"), f)
end

# ---------------------------------------------------------------------------------------
# 3. Physical mass drift.
# ---------------------------------------------------------------------------------------
function figure3()
    f = Figure(; size=(1000, 400))
    for (col, (c, ref, ttl)) in enumerate(((coarse, exact_mass, "64² cells / 16² elements"),
                                           (fine, exact_mass, "128² cells / 32² elements")))
        ax = Axis(f[1, col]; xlabel="t", ylabel="M_h(t)/M_exact - 1",
                  title=ttl, yscale=log10)
        for k in 1:4
            lines!(ax, c.times, [abs(m[k] / ref - 1) + 1e-16 for m in c.masses];
                   color=COLORS[k], label=METHODS[k], linewidth=1.8)
        end
        col == 1 && axislegend(ax; position=:rb, framevisible=false, labelsize=10)
    end
    Label(f[0, :], "Physical mass drift, |M_h/M_exact - 1| (exact material transport ⇒ 0)", fontsize=13)
    save(joinpath(OUTDIR, "3_mass_drift.png"), f)
end

# ---------------------------------------------------------------------------------------
# 4. Accuracy against DOF count and 5. against solver time.
# ---------------------------------------------------------------------------------------
function figure45()
    summ = [s for s in data.summary if s.label in (:ncells64, :ncells128)]
    f = Figure(; size=(1000, 400))
    ax = Axis(f[1, 1]; xlabel="degrees of freedom", ylabel="relative L1 error at t = 1",
              xscale=log10, yscale=log10, title="accuracy vs resolution")
    ax2 = Axis(f[1, 2]; xlabel="solver time for the whole run (s)", ylabel="relative L1 error at t = 1",
               xscale=log10, yscale=log10, title="accuracy vs cost")
    for k in 1:4
        dof = [k == 4 ? s.ndof_trixi : s.ndof_ed for s in summ]
        times = [s.time[k + 1] for s in summ]
        err = [s.l1[k] for s in summ]
        scatterlines!(ax, dof, err; color=COLORS[k], markersize=8, linewidth=1.8, label=METHODS[k])
        scatterlines!(ax2, times, err; color=COLORS[k], markersize=8, linewidth=1.8, label=METHODS[k])
    end
    axislegend(ax; position=:lb, framevisible=false, labelsize=10)
    axislegend(ax2; position=:lb, framevisible=false, labelsize=10)
    save(joinpath(OUTDIR, "4_accuracy_vs_dof.png"), f)
end

# ---------------------------------------------------------------------------------------
# 6. Numerical boundary layer: interface zoom and error against φ/h.
# ---------------------------------------------------------------------------------------
function figure6()
    c = fine
    s = snapAt(c, 1.0)
    n = c.n
    ax_axis = range(0.0, 1.0; length=n)
    f = Figure(; size=(1150, 900))

    # (a) zoom on the inner interface at t = 1
    xlim = (0.42, 0.56)
    ylim = (0.50, 0.64)
    for (row, key) in enumerate((:weno, :cip, :map, :trixi))
        rho = getproperty(s, key)
        ax = Axis(f[row, 1]; aspect=DataAspect(), title=METHODS[row], titlesize=11)
        hm = heatmap!(ax, ax_axis, ax_axis, rho; colorrange=(1.00, 1.13), colormap=:viridis)
        contour!(ax, ax_axis, ax_axis, s.phi; levels=[0.0], color=:white, linewidth=2.0)
        th = range(0, 2π; length=400)
        lines!(ax, 0.59 .+ 0.12 .* cos.(th), 0.55 .+ 0.12 .* sin.(th); color=:red, linestyle=:dash)
        limits!(ax, xlim..., ylim...)
        row == 4 && Colorbar(f[4, 1], hm; vertical=false, label="ρ")
        hidespines!(ax); hidedecorations!(ax)
    end

    # (b) mean nodal error against signed distance in cells
    centers = (-5.0, -3.0, -1.5, -0.5, 0.5, 1.5, 3.0)
    ax = Axis(f[1:3, 2:3]; xlabel="φ/h  (negative inside the material)",
              ylabel="mean |ρ_h - ρ_exact| per node", yscale=log10)
    for k in 1:4
        ys = [s.bins[k][b][3] + 1e-18 for b in 1:7]
        scatterlines!(ax, collect(centers), ys; color=COLORS[k], label=METHODS[k], linewidth=1.8)
    end
    axislegend(ax; position=:rt, framevisible=false, labelsize=10)
    Label(f[0, :], "Numerical boundary layer at t = 1: interface zoom (left, dashed = analytic circle) and nodal error by signed distance",
          fontsize=12)
    save(joinpath(OUTDIR, "5_boundary_layer.png"), f)
end

# ---------------------------------------------------------------------------------------
# 7. Constant-density test: partition of unity and support-boundary defects.
# ---------------------------------------------------------------------------------------
function figure7()
    c = finec
    s = snapAt(c, 1.0)
    n = c.n
    ax_axis = range(0.0, 1.0; length=n)
    f = Figure(; size=(1150, 420))
    for (col, key) in enumerate((:weno, :cip, :map, :trixi))
        rho = getproperty(s, key)
        err = [s.phi[i, j] < 0 ? max(abs(rho[i, j] - 1.0), 1e-17) : NaN for i in 1:n, j in 1:n]
        ax = Axis(f[1, col]; aspect=DataAspect(), title=METHODS[col], titlesize=12)
        hm = heatmap!(ax, ax_axis, ax_axis, err; colorrange=(1e-17, 1e-1), colorscale=log10,
                      colormap=:magma, lowclip=:white)
        contour!(ax, ax_axis, ax_axis, s.phi; levels=[0.0], color=:cyan, linewidth=1.2)
        Colorbar(f[2, col], hm; vertical=false, label="|ρ_h - 1|")
        hidespines!(ax); hidedecorations!(ax)
    end
    Label(f[0, :], "Constant state ρ ≡ 1 after transport to t = 1 (128²/32²): partition of unity on the moving domain",
          fontsize=13)
    save(joinpath(OUTDIR, "6_constant_state.png"), f)
end

# ---------------------------------------------------------------------------------------
# 8. Maps as a linear operator: the two identities of §5, measured per step.
# ---------------------------------------------------------------------------------------
function figure8()
    c = fine
    f = Figure(; size=(1000, 400))
    ax = Axis(f[1, 1]; xlabel="t", ylabel="relative change of the node-set sum",
              title="unweighted sum over the fixed node set", yscale=log10)
    for (k, key) in enumerate((1, 2, 3))
        lines!(ax, c.times, [abs(m[key] / c.masses[1][key] - 1) + 1e-17 for m in c.masses];
               color=COLORS[key], label=METHODS[key], linewidth=1.8)
    end
    axislegend(ax; framevisible=false, labelsize=10)
    ax2 = Axis(f[1, 2]; xlabel="t", ylabel="|Σq_after/Σq_before - 1| per step",
               title="TransportMap invariant on its own support", yscale=log10)
    lines!(ax2, c.times[2:end], [abs(d) + 1e-17 for d in c.mapdrift]; color=:firebrick, linewidth=1.8)
    hlines!(ax2, [eps()]; color=:black, linestyle=:dash, label="eps()")
    axislegend(ax2; framevisible=false, labelsize=10)
    Label(f[0, :], "Two different conservation statements (128²/32²; φ(t) is the same for every method)", fontsize=13)
    save(joinpath(OUTDIR, "7_conservation_split.png"), f)
end

figure1(); figure2(); figure3(); figure45(); figure6(); figure7(); figure8()
println("figures written to $OUTDIR")

# --- compact numeric digest used by REPORT.md
for s in data.summary
    @printf("%-16s ndof=%-6d L1=%.3e L2=%.3e Linf=%.3e mass=%.3e time=%.3f s\n",
            String(s.label), s.ndof_trixi, s.l1[4], s.l2[4], s.linf[4], s.mass_drift[4], s.time[5])
end
