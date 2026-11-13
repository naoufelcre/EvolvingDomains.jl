# Canonical 3D Zalesak sphere: Field WENO5 / SSP-RK3, one revolution at 50³.
# Runs with any project: the pinned validation environment is activated below.
using Pkg
Pkg.activate(joinpath(@__DIR__, "..", "env"); io=devnull)
using EvolvingDomains
using StaticArrays: SVector
include(joinpath(@__DIR__, "zalesak.jl"))
include(joinpath(@__DIR__, "measure.jl"))
include(joinpath(@__DIR__, "artifacts.jl"))

n = 50                                   # cells per axis
rotations = 1
cfl = 0.4

axis = collect(range(-1.0, 1.0; length=n+1))

# initial level set and the prescribed rigid rotation
u0 = [sphere(x, y, z) for x in axis, y in axis, z in axis]
velocity = vec([sphere_velocity(x, y, z) for x in axis, y in axis, z in axis])

# CFL: v = π(y,-x,0) gives |v₁|+|v₂| ≤ 2π on the box; πn is the bound used here
rate = π*n
steps_per_revolution = 4 * ceil(Int, PERIOD_3D*rate/(4cfl))
dt = PERIOD_3D/steps_per_revolution
steps = rotations * steps_per_revolution

# discretization and time loop, snapshotting every quarter revolution
field = Field(Discretization(Dict(:init => Grid((axis, axis, axis)))),
              Dict(:init => copy(u0)))
cache = WENO5Cache()
eulerianTransport!(Field(field.discretization, Dict(:init => copy(u0))), velocity, dt;
                   level_id=:init, cache)          # warm up, excluded from timing
snapshots, times = [copy(u0)], [0.0]
elapsed = @elapsed for step in 1:steps
    eulerianTransport!(field, velocity, dt; level_id=:init, cache)
    if step % (steps_per_revolution ÷ 4) == 0
        push!(times, step*dt)
        push!(snapshots, copy(field.data[:init]))
    end
end
exact = [[sphere_rotated(x, y, z, t) for x in axis, y in axis, z in axis] for t in times]

# measurements against the analytic solution
initial = phase_measure(first(snapshots), (axis, axis, axis), sphere_rotated, 0.0; samples=2)
final = phase_measure(last(snapshots), (axis, axis, axis), sphere_rotated, last(times);
                      samples=2)
norms = field_norms(last(snapshots), last(exact), (axis, axis, axis), CENTER_3D, 0.65)
metrics = (cells=n, rotations=rotations, steps=steps, dt=dt, seconds=elapsed,
           measure_drift=(final.measure - initial.measure)/initial.measure,
           symmetric_difference=final.symmetric_difference,
           initial_measure_error=abs(initial.measure - EXACT_VOLUME_3D),
           l1=norms.l1, linf=norms.linf, exact_volume=EXACT_VOLUME_3D)

# acceptance gates for this variant
within(metrics.measure_drift, "measure_drift", -0.05, 0.05)
within(metrics.symmetric_difference, "symmetric_difference", 0.0, 0.15EXACT_VOLUME_3D)
within(metrics.initial_measure_error, "initial_measure_error", 0.0, 0.03EXACT_VOLUME_3D)
within(metrics.l1, "l1", 0.0, 0.01)
within(metrics.linf, "linf", 0.0, 0.09)

dir = length(ARGS) >= 1 ? abspath(ARGS[1]) : mktempdir(mkpath(joinpath(@__DIR__, "output"));
                                                       prefix="run-", cleanup=false)
save_figure_or_skip(() -> figure_3d((axis, axis, axis), CENTER_3D, times, snapshots, exact,
                                    metrics, "Canonical 3D Zalesak — $(n)³, $(rotations) turn"),
                    dir, "weno_3d")
report(dir, "3D Zalesak: $(n)³ cells, $(rotations) revolution, WENO5 + SSP-RK3", metrics,
       (joinpath(@__DIR__, "zalesak.jl"), joinpath(@__DIR__, "measure.jl"),
        @__FILE__))
