# 2D Zalesak disk: Field WENO5 / SSP-RK3, ten rotations at 100 cells per axis.
# Runs with any project: the pinned validation environment is activated below.
using Pkg
Pkg.activate(joinpath(@__DIR__, "..", "env"); io=devnull)
using EvolvingDomains
using StaticArrays: SVector
include(joinpath(@__DIR__, "zalesak.jl"))
include(joinpath(@__DIR__, "measure.jl"))
include(joinpath(@__DIR__, "artifacts.jl"))

n = 100                                  # cells per axis
rotations = 10
cfl = 0.4

axis = collect(range(0.0, 1.0; length=n+1))

# initial level set and the prescribed rigid rotation
u0 = [disk(x, y) for x in axis, y in axis]
velocity = vec([disk_velocity(x, y) for x in axis, y in axis])

# CFL: one rotation takes PERIOD_2D, split into equal steps
rate = n * maximum(v -> abs(v[1]) + abs(v[2]), velocity)
steps_per_period = 4 * ceil(Int, PERIOD_2D*rate/(4cfl))
dt = PERIOD_2D/steps_per_period
steps = rotations * steps_per_period

# discretization and time loop, snapshotting every quarter turn
field = Field(Discretization(Dict(:init => Grid((axis, axis)))), Dict(:init => copy(u0)))
cache = WENO5Cache()
eulerianTransport!(Field(field.discretization, Dict(:init => copy(u0))), velocity, dt;
                   level_id=:init, cache)          # warm up, excluded from timing
snapshots, times = [copy(u0)], [0.0]
elapsed = @elapsed for step in 1:steps
    eulerianTransport!(field, velocity, dt; level_id=:init, cache)
    if step % (steps_per_period ÷ 4) == 0
        push!(times, step*dt)
        push!(snapshots, copy(field.data[:init]))
    end
end
exact = [[disk_rotated(x, y, t) for x in axis, y in axis] for t in times]

# measurements against the analytic solution
initial = phase_measure(first(snapshots), (axis, axis), disk_rotated, 0.0)
final = phase_measure(last(snapshots), (axis, axis), disk_rotated, last(times))
norms = field_norms(last(snapshots), last(exact), (axis, axis), (0.5, 0.5), 0.45)
distance = interface_distance(interface_points(last(snapshots), axis),
                              interface_points(last(exact), axis))
metrics = (cells=n, rotations=rotations, steps=steps, dt=dt, seconds=elapsed,
           measure_drift=(final.measure - initial.measure)/initial.measure,
           symmetric_difference=final.symmetric_difference, l1=norms.l1, l2=norms.l2,
           linf=norms.linf, interface_distance=distance, exact_area=EXACT_AREA_2D)

# acceptance gates for this variant
within(metrics.measure_drift, "measure_drift", -0.04, 0.04)
within(metrics.symmetric_difference, "symmetric_difference", 0.0, 0.008)
within(metrics.l1, "l1", 0.0, 0.002)
within(metrics.l2, "l2", 0.0, 0.005)
within(metrics.linf, "linf", 0.0, 0.04)
within(metrics.interface_distance, "interface_distance", 0.0, 0.04)

dir = length(ARGS) >= 1 ? abspath(ARGS[1]) : mktempdir(mkpath(joinpath(@__DIR__, "output"));
                                                       prefix="run-", cleanup=false)
save_figure(figure_2d(axis, times, snapshots, exact, metrics), dir, "weno_2d")
report(dir, "2D Zalesak: $(n)² cells, $(rotations) rotations, WENO5 + SSP-RK3", metrics,
       (joinpath(@__DIR__, "zalesak.jl"), joinpath(@__DIR__, "measure.jl"),
        @__FILE__))
