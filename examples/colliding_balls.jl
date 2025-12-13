# Colliding Balls Visualization Example
# ======================================
using EvolvingDomains
using Gridap
using GridapEmbedded
using LevelSetMethods
using GLMakie
using FFMPEG

# =============================================================================
# Setup
# =============================================================================

domain = (0.0, 2.0, 0.0, 2.0)
partition = (80, 80)
model = CartesianDiscreteModel(domain, partition)

center_top = (1.0, 1.3)
center_bottom = (1.0, 0.7)
radius = 0.2

function ϕ0(x)
    d_top = sqrt((x[1] - center_top[1])^2 + (x[2] - center_top[2])^2) - radius
    d_bottom = sqrt((x[1] - center_bottom[1])^2 + (x[2] - center_bottom[2])^2) - radius
    return min(d_top, d_bottom)
end

# Constant velocity toward y=1.0
velocity(x) = x[2] > 1.0 ? (0.0, -0.5) : (0.0, 0.5)

evolver = LevelSetMethodsEvolver(;
    bg_model = model,
    initial_ls = ϕ0,
    velocity = velocity,
    spatial_scheme = :WENO5,
    integrator = :RK3,
    bc = :Neumann
)

eg = EvolvingDiscreteGeometry(evolver, model)

info = grid_info(eg)
nx, ny = info.dims
x_range = range(info.origin[1], step=info.spacing[1], length=nx)
y_range = range(info.origin[2], step=info.spacing[2], length=ny)

# =============================================================================
# Create Animation Frames
# =============================================================================

# Fresh evolver for animation
evolver2 = LevelSetMethodsEvolver(;
    bg_model = model,
    initial_ls = ϕ0,
    velocity = velocity,
    spatial_scheme = :WENO5,
    integrator = :RK3,
    bc = :Neumann
)
eg2 = EvolvingDiscreteGeometry(evolver2, model)

# Create frames directory
frames_dir = joinpath(@__DIR__, "colliding_frames")
rm(frames_dir, force=true, recursive=true)
mkpath(frames_dir)

Δt = 0.01
substeps = 2
nframes = 30

for frame in 1:nframes
    t = EvolvingDomains.current_time(eg2)

    # Create figure
    fig = Figure(size = (700, 600))
    ax = Axis(fig[1, 1],
              title = string("Colliding Balls (t = ", round(t, digits=2), ")"),
              xlabel = "x", ylabel = "y",
              aspect = DataAspect())

    ϕ = current_levelset(eg2)
    ϕ_2d = reshape(ϕ, (nx, ny))

    heatmap!(ax, x_range, y_range, ϕ_2d; colormap = :RdBu, colorrange = (-0.4, 0.4))
    contour!(ax, x_range, y_range, ϕ_2d; levels = [0.0], linewidth = 4, color = :black)

    # Save frame
    filename = joinpath(frames_dir, string("frame_", lpad(frame, 3, "0"), ".png"))
    save(filename, fig)

    if frame % 10 == 1
        println("  Frame ", frame, "/", nframes, " (t=", round(t, digits=2), ")")
    end

    # Advance
    for _ in 1:substeps
        EvolvingDomains.advance!(eg2, Δt)
    end
end

println("✓ Frames saved to: ", frames_dir)

# =============================================================================
# Combine frames into video using FFMPEG.jl
# =============================================================================
output_video = joinpath(@__DIR__, "colliding_balls.mp4")
input_pattern = joinpath(frames_dir, "frame_%03d.png")
FFMPEG.exe(`-y -framerate 15 -i $input_pattern -c:v libx264 -pix_fmt yuv420p $output_video`)

println("Video saved: ", output_video)

# Evolve original eg to see final state
for _ in 1:60
    advance!(eg, Δt)
end

center_idx = div(nx, 2) + 1 + (div(ny, 2)) * nx
final_phi = current_levelset(eg)[center_idx]

fig2 = plot_levelset(eg; title=string("Final (t=", round(EvolvingDomains.current_time(eg), digits=1), ")"))
display(fig2)
