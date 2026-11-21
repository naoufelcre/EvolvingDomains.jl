# Canonical 3D Zalesak sphere: Field WENO5 and Trixi DGSEM at the same physical times.
# Runs with any project: the pinned validation environment is activated below.
using Pkg
Pkg.activate(joinpath(@__DIR__, "..", "env"); io=devnull)
using EvolvingDomains
using Trixi
using StaticArrays: SVector
include(joinpath(@__DIR__, "zalesak.jl"))
include(joinpath(@__DIR__, "measure.jl"))
include(joinpath(@__DIR__, "artifacts.jl"))
include(joinpath(@__DIR__, "TrixiZalesak3D.jl"))

n = 50                                   # cells per axis, Field
rotations = 1
cfl = 0.4
dg_cells = 16                            # elements per axis, Trixi
polydeg = 3

axis = collect(range(-1.0, 1.0; length=n+1))
times = collect(range(0.0, rotations*PERIOD_3D; length=4rotations+1))
u0 = [sphere(x, y, z) for x in axis, y in axis, z in axis]
velocity = vec([sphere_velocity(x, y, z) for x in axis, y in axis, z in axis])

# Field WENO5: πn is the CFL bound used for v = π(y,-x,0)
rate = π*n
steps_per_revolution = 4 * ceil(Int, PERIOD_3D*rate/(4cfl))
dt = PERIOD_3D/steps_per_revolution
field = Field(Discretization(Dict(:init => Grid((axis, axis, axis)))),
              Dict(:init => copy(u0)))
cache = WENO5Cache()
eulerianTransport!(Field(field.discretization, Dict(:init => copy(u0))), velocity, dt;
                   level_id=:init, cache)
snapshots, ed_times = [copy(u0)], [0.0]
ed_seconds = @elapsed for step in 1:rotations*steps_per_revolution
    eulerianTransport!(field, velocity, dt; level_id=:init, cache)
    if step % (steps_per_revolution ÷ 4) == 0
        push!(ed_times, step*dt)
        push!(snapshots, copy(field.data[:init]))
    end
end

# Trixi DGSEM on the same grid, sampled at the same times
dg = TrixiZalesak3D.run_trixi3d(Grid((axis, axis, axis)), sphere, times;
                                cells=dg_cells, polydeg, cfl=0.25)
exact = [[sphere_rotated(x, y, z, t) for x in axis, y in axis, z in axis] for t in times]

# measurements for each solver, and of their difference on the common grid
ed_ref = phase_measure(first(snapshots), (axis, axis, axis), sphere_rotated, 0.0; samples=2)
ed = phase_measure(last(snapshots), (axis, axis, axis), sphere_rotated, last(ed_times);
                   samples=2)
dg_ref = phase_measure(first(dg.snapshots), (axis, axis, axis), sphere_rotated, 0.0;
                       samples=2)
dg_final = phase_measure(last(dg.snapshots), (axis, axis, axis), sphere_rotated, last(times);
                         samples=2)
ed_norms = field_norms(last(snapshots), last(exact), (axis, axis, axis), CENTER_3D, 0.65)
dg_norms = field_norms(last(dg.snapshots), last(exact), (axis, axis, axis), CENTER_3D, 0.65)
paired = field_norms(last(snapshots), last(dg.snapshots), (axis, axis, axis), CENTER_3D, 0.65)
metrics = (cells=n, dg_cells=dg_cells, polydeg=polydeg, rotations=rotations,
           ed_seconds=ed_seconds, dg_seconds=dg.elapsed,
           ed_measure_drift=(ed.measure - ed_ref.measure)/ed_ref.measure,
           dg_measure_drift=(dg_final.measure - dg_ref.measure)/dg_ref.measure,
           ed_symmetric_difference=ed.symmetric_difference,
           dg_symmetric_difference=dg_final.symmetric_difference,
           ed_l1=ed_norms.l1, ed_linf=ed_norms.linf, dg_l1=dg_norms.l1,
           dg_linf=dg_norms.linf, solver_difference_l1=paired.l1)

# acceptance gates: each solver against the analytic solution, plus their difference
within(metrics.ed_measure_drift, "ed_measure_drift", -0.05, 0.05)
within(metrics.ed_symmetric_difference, "ed_symdiff", 0.0, 0.15EXACT_VOLUME_3D)
within(metrics.dg_measure_drift, "dg_measure_drift", -0.03, 0.03)
within(metrics.dg_symmetric_difference, "dg_symdiff", 0.0, 0.04)
within(metrics.dg_l1, "dg_l1", 0.0, 0.005)
within(metrics.dg_linf, "dg_linf", 0.0, 0.05)
within(metrics.solver_difference_l1, "solver_difference_l1", 0.0, 0.02)

dir = length(ARGS) >= 1 ? abspath(ARGS[1]) : mktempdir(mkpath(joinpath(@__DIR__, "output"));
                                                       prefix="run-", cleanup=false)
save_figure_or_skip(() -> figure_comparison_3d((axis, axis, axis), CENTER_3D, snapshots,
                                               dg.snapshots, exact, metrics,
                                               "Field WENO5 vs Trixi DGSEM — canonical 3D Zalesak"),
                    dir, "trixi_3d")
report(dir, "3D Zalesak: Field $(n)³ vs Trixi $(dg_cells)³ degree-$(polydeg), " *
            "$(rotations) revolution", metrics,
       (joinpath(@__DIR__, "zalesak.jl"), joinpath(@__DIR__, "measure.jl"),
        joinpath(@__DIR__, "TrixiZalesak3D.jl"), @__FILE__))
