module TestGeometryEvolution

using EvolvingDomains
using EvolvingDomains.Geometric
using EvolvingDomains.Kinematic
using Gridap
using Gridap.TensorValues

"""
    run_zalesak(; n=100, nTime=10000, plot_stride=10, on_frame=(args...; kwargs...) -> nothing)

Advect the Zalesak disk by solid-body rotation. Headless unless a frame
callback is supplied; `on_frame(geom; label=...)` fires every `plot_stride`
steps. Returns the evolved geometry.
"""
function run_zalesak(; n=100, nTime=10000, plot_stride=10,
    on_frame=(args...; kwargs...) -> nothing)
    # Define Grid (n x n on [0,1]x[0,1])
    grid = CartesianDiscreteModel((0, 1, 0, 1), (n, n))
    center = VectorValue(0.5, 0.5) # Center of the Grid

    # Create empty geometry container
    geom = EvolvingDiscreteGeometry(grid)

    # Define Geometry (Zalesak Disk)
    disk = EvolvingDomains.Geometric.Circle(VectorValue(0.5, 0.75), 0.15)
    slot = EvolvingDomains.Geometric.Rectangle(VectorValue(0.475, 0.60), VectorValue(0.525, 0.75))
    zalesak_geo = setdiff(disk, slot)

    # Initialize Level Set Field
    pts = Gridap.Geometry.get_node_coordinates(grid)
    phi_values = map(x -> zalesak_geo(x), pts)
    set_levelset!(geom, vec(collect(phi_values)))

    # Define Velocity as Solid Body Rotation.
    ω = 20π # Angular velocity
    function rigid_rotation(x)
        rx = x[1] - center[1]
        ry = x[2] - center[2]
        VectorValue(-ω * ry, ω * rx)
    end
    #Set up
    vel = StaticFunctionVelocity(rigid_rotation)

    # ====== Time Loop ======
    timeHorizon = 1
    Δt = timeHorizon / nTime

    on_frame(geom; label="Iteration 0 / $nTime")

    info = grid_info(grid)
    v_field = sample_velocity(vel, info, 0)
    started = time_ns()

    for i in 1:nTime
        advance!(geom, v_field, Δt)
        i % plot_stride == 0 && on_frame(geom; label="Iteration $i / $nTime   mean = $(round((time_ns() - started) / (1e6 * i), digits=2)) ms/iteration")
    end

    return geom
end

end # module

if abspath(PROGRAM_FILE) == abspath(@__FILE__)
    TestGeometryEvolution.run_zalesak()
end
