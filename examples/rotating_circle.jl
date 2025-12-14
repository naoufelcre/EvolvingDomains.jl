# Rotating Circle Example
# =======================
# A circle advecting across the domain with constant velocity

using EvolvingDomains
using Gridap, GridapEmbedded
using LevelSetMethods
using GLMakie
using FFMPEG

# =============================================================================
# Setup
# =============================================================================
domain = (-1.0, 1.0, -1.0, 1.0)
partition = (100, 100)
model = CartesianDiscreteModel(domain, partition)

# Initial circle at left side
R = 0.2
center = (0.3, 0.5)
ϕ₀(x) = R - sqrt((x[1]-center[1])^2 + (x[2]-center[2])^2)

ω = 2π  # Angular velocity (one revolution in T=1)
u(x) = (-ω * x[2], ω * x[1])

# Create evolver
evolver = LevelSetMethodsEvolver(;
    bg_model = model,
    initial_ls = ϕ₀,
    velocity = u,
    spatial_scheme = :WENO5,
    bc = :Neumann
)
eg = EvolvingDiscreteGeometry(evolver, model)

info = EvolvingDomains.grid_info(eg)
nx, ny = info.dims
x_range = range(info.origin[1], step=info.spacing[1], length=nx)
y_range = range(info.origin[2], step=info.spacing[2], length=ny)

# =============================================================================
# Initial State
# =============================================================================
println("Rotating Circle Example")
println("=" ^ 40)
fig_initial = plot_levelset(eg; title = "Initial Configuration")
display(fig_initial)

# =============================================================================
# Create Animation Frames
# =============================================================================
frames_dir = joinpath(@__DIR__, "rotating_circle_frames")
rm(frames_dir, force=true, recursive=true)
mkpath(frames_dir)

Δt = 0.01
substeps = 2
nframes = 50

for frame in 1:nframes
    t = EvolvingDomains.current_time(eg)

    # Create figure
    fig = Figure(size = (800, 400))
    ax = Axis(fig[1, 1],
              title = string("Advecting Circle (t = ", round(t, digits=2), ")"),
              xlabel = "x", ylabel = "y",
              aspect = DataAspect())

    ϕ = EvolvingDomains.current_levelset(eg)
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
        EvolvingDomains.advance!(eg, Δt)
    end

    # Reinitialize periodically
    if frame % 5 == 0
        EvolvingDomains.reinitialize!(eg)
    end
end

println("✓ Frames saved to: ", frames_dir)

# =============================================================================
# Combine frames into video using FFMPEG.jl
# =============================================================================
output_video = joinpath(@__DIR__, "rotating_circle.mp4")
input_pattern = joinpath(frames_dir, "frame_%03d.png")
FFMPEG.exe(`-y -framerate 15 -i $input_pattern -c:v libx264 -pix_fmt yuv420p $output_video`)

println("✓ Video saved: ", output_video)

# =============================================================================
# Final State
# =============================================================================
println("\nFinal time: ", EvolvingDomains.current_time(eg))
fig_final = plot_levelset(eg; title = "Advecting Circle (Final)")
display(fig_final)
