# Tilted-axis 3D Zalesak sphere: ten rotations at 96³, the stress variant of the same case.
# Runs with any project: the pinned validation environment is activated below.
using Pkg
Pkg.activate(joinpath(@__DIR__, "..", "env"); io=devnull)
using EvolvingDomains
using StaticArrays: SVector
include(joinpath(@__DIR__, "zalesak.jl"))
include(joinpath(@__DIR__, "measure.jl"))
include(joinpath(@__DIR__, "artifacts.jl"))

n = 96                                   # cells per axis; the slot is ~4.8 cells wide
rotations = 10
cfl = 0.4

axis = collect(range(0.0, 1.0; length=n+1))

# initial level set and the prescribed tilted rigid rotation about the box center
u0 = [tilted_sphere(x, y, z) for x in axis, y in axis, z in axis]
velocity = vec([tilted_velocity(x, y, z) for x in axis, y in axis, z in axis])

# CFL from the measured |v|₁ bound on this grid
rate = n * maximum(v -> abs(v[1]) + abs(v[2]) + abs(v[3]), velocity)
steps_per_period = 4 * ceil(Int, PERIOD_TILT*rate/(4cfl))
dt = PERIOD_TILT/steps_per_period
steps = rotations * steps_per_period

# discretization and time loop, snapshotting every quarter turn
field = Field(Discretization(Dict(:init => Grid((axis, axis, axis)))),
              Dict(:init => copy(u0)))
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
exact = [[tilted_rotated(x, y, z, t) for x in axis, y in axis, z in axis] for t in times]

# measurements against the analytic solution
initial = phase_measure(first(snapshots), (axis, axis, axis), tilted_rotated, 0.0; samples=2)
final = phase_measure(last(snapshots), (axis, axis, axis), tilted_rotated, last(times);
                      samples=2)
norms = field_norms(last(snapshots), last(exact), (axis, axis, axis), (0.5, 0.5, 0.5), 0.45)
metrics = (cells=n, rotations=rotations, steps=steps, dt=dt, seconds=elapsed,
           measure_drift=(final.measure - initial.measure)/initial.measure,
           symmetric_difference=final.symmetric_difference, l1=norms.l1, linf=norms.linf,
           exact_volume=EXACT_VOLUME_TILT)

# acceptance gates for this variant; the 48³ arm was deliberately dropped
within(metrics.measure_drift, "measure_drift", -0.05, 0.05)
within(metrics.symmetric_difference, "symmetric_difference", 0.0, 0.01)
within(metrics.l1, "l1", 0.0, 0.002)
within(metrics.linf, "linf", 0.0, 0.04)

dir = length(ARGS) >= 1 ? abspath(ARGS[1]) : mktempdir(mkpath(joinpath(@__DIR__, "output"));
                                                       prefix="run-", cleanup=false)
save_figure_or_skip(() -> figure_3d((axis, axis, axis), (0.5, 0.5, 0.5), times, snapshots,
                                    exact, metrics, "3D Zalesak stress test — $(n)³, " *
                                    "$(rotations) rotations, tilted axis"),
                    dir, "weno_3d_ten")
report(dir, "3D Zalesak stress: $(n)³ cells, $(rotations) rotations, tilted axis", metrics,
       (joinpath(@__DIR__, "zalesak.jl"), joinpath(@__DIR__, "measure.jl"),
        @__FILE__))
