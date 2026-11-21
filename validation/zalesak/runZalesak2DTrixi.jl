# 2D Zalesak disk: Field WENO5 and Trixi DGSEM at the same physical times.
# Runs with any project: the pinned validation environment is activated below.
using Pkg
Pkg.activate(joinpath(@__DIR__, "..", "env"); io=devnull)
using EvolvingDomains
using Trixi
using StaticArrays: SVector
include(joinpath(@__DIR__, "zalesak.jl"))
include(joinpath(@__DIR__, "measure.jl"))
include(joinpath(@__DIR__, "artifacts.jl"))
include(joinpath(@__DIR__, "TrixiZalesakSolver.jl"))

n = 100                                  # cells per axis, Field
rotations = 10
cfl = 0.4
dg_cells = 32                            # elements per axis, Trixi
polydeg = 3

axis = collect(range(0.0, 1.0; length=n+1))
times = collect(range(0.0, rotations*PERIOD_2D; length=4rotations+1))
u0 = [disk(x, y) for x in axis, y in axis]
velocity = vec([disk_velocity(x, y) for x in axis, y in axis])

# Field WENO5: CFL from the nodal velocity bound
rate = n * maximum(v -> abs(v[1]) + abs(v[2]), velocity)
steps_per_period = 4 * ceil(Int, PERIOD_2D*rate/(4cfl))
dt = PERIOD_2D/steps_per_period
field = Field(Discretization(Dict(:init => Grid((axis, axis)))), Dict(:init => copy(u0)))
cache = WENO5Cache()
eulerianTransport!(Field(field.discretization, Dict(:init => copy(u0))), velocity, dt;
                   level_id=:init, cache)
snapshots, ed_times = [copy(u0)], [0.0]
ed_seconds = @elapsed for step in 1:rotations*steps_per_period
    eulerianTransport!(field, velocity, dt; level_id=:init, cache)
    if step % (steps_per_period ÷ 4) == 0
        push!(ed_times, step*dt)
        push!(snapshots, copy(field.data[:init]))
    end
end

# Trixi DGSEM on the same grid, sampled at the same times
dg = TrixiZalesakSolver.run_trixi(Grid((axis, axis)), disk, times;
                                  cells=dg_cells, polydeg, cfl=0.3)
exact = [[disk_rotated(x, y, t) for x in axis, y in axis] for t in times]

# measurements for each solver, and of their difference on the common grid
ed_ref = phase_measure(first(snapshots), (axis, axis), disk_rotated, 0.0)
ed = phase_measure(last(snapshots), (axis, axis), disk_rotated, last(ed_times))
dg_ref = phase_measure(first(dg.snapshots), (axis, axis), disk_rotated, 0.0)
dg_final = phase_measure(last(dg.snapshots), (axis, axis), disk_rotated, last(times))
ed_norms = field_norms(last(snapshots), last(exact), (axis, axis), (0.5, 0.5), 0.45)
dg_norms = field_norms(last(dg.snapshots), last(exact), (axis, axis), (0.5, 0.5), 0.45)
paired = field_norms(last(snapshots), last(dg.snapshots), (axis, axis), (0.5, 0.5), 0.45)
metrics = (cells=n, dg_cells=dg_cells, polydeg=polydeg, rotations=rotations,
           ed_seconds=ed_seconds, dg_seconds=dg.elapsed,
           ed_measure_drift=(ed.measure - ed_ref.measure)/ed_ref.measure,
           dg_measure_drift=(dg_final.measure - dg_ref.measure)/dg_ref.measure,
           ed_symmetric_difference=ed.symmetric_difference,
           dg_symmetric_difference=dg_final.symmetric_difference,
           ed_l1=ed_norms.l1, dg_l1=dg_norms.l1, solver_difference_l1=paired.l1)

# acceptance gates: each solver against the analytic solution, plus their difference
within(metrics.ed_measure_drift, "ed_measure_drift", -0.04, 0.04)
within(metrics.ed_symmetric_difference, "ed_symdiff", 0.0, 0.008)
within(metrics.dg_measure_drift, "dg_measure_drift", -0.3, 0.3)
within(metrics.dg_symmetric_difference, "dg_symdiff", 0.0, EXACT_AREA_2D/2)
within(metrics.dg_l1, "dg_l1", 0.0, 0.02)
within(metrics.solver_difference_l1, "solver_difference_l1", 0.0, 0.02)

dir = length(ARGS) >= 1 ? abspath(ARGS[1]) : mktempdir(mkpath(joinpath(@__DIR__, "output"));
                                                       prefix="run-", cleanup=false)
save_figure(figure_comparison_2d(axis, times, snapshots, dg.snapshots, exact, metrics,
                                 "Field WENO5 vs Trixi DGSEM — $(rotations) rotations"),
            dir, "trixi_2d")
report(dir, "2D Zalesak: Field $(n)² vs Trixi $(dg_cells)² degree-$(polydeg), " *
            "$(rotations) rotations", metrics,
       (joinpath(@__DIR__, "zalesak.jl"), joinpath(@__DIR__, "measure.jl"),
        joinpath(@__DIR__, "TrixiZalesakSolver.jl"), @__FILE__))
