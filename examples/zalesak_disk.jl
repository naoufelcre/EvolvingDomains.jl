# Zalesak's Rotating Disk
# =======================
# Classic benchmark: slotted disk under rigid body rotation

using EvolvingDomains
using Gridap
using GridapEmbedded
using LevelSetMethods
using GLMakie
using FFMPEG

# =============================================================================
# Domain Setup
# =============================================================================
domain = (-1.5, 1.5, -1.5, 1.5)
partition = (100, 100)
model = CartesianDiscreteModel(domain, partition)

# =============================================================================
# Zalesak Disk Level Set
# =============================================================================
function zalesak_disk(x)
    center = (-0.75, 0.0)
    radius = 0.5

    # Circle SDF: positive inside, negative outside
    d_circle = radius - sqrt((x[1] - center[1])^2 + (x[2] - center[2])^2)

    # Rectangular slot: centered on circle center, extends upward
    h = 1.0   # slot height
    w = 0.2   # slot width
    
    # Slot bounds (centered horizontally on circle center, bottom at circle center)
    xmin = center[1] - w/2
    xmax = center[1] + w/2
    ymin = center[2]
    ymax = center[2] + h
    
    # Slot SDF: positive inside slot, negative outside
    dx = min(x[1] - xmin, xmax - x[1])
    dy = min(x[2] - ymin, ymax - x[2])
    d_slot = min(dx, dy)  # positive inside, negative outside
    
    # Zalesak disk = circle - slot
    # ϕ > 0 inside the slotted disk (positive inside convention)
    # Inside disk AND outside slot: min(d_circle, -d_slot)
    return min(d_circle, -d_slot)
end

# =============================================================================
# Rigid Body Rotation Velocity
# =============================================================================
ω = 2π  # Angular velocity (one revolution in T=1)
velocity(x) = (-ω * x[2], ω * x[1])

# =============================================================================
# Create Evolver
# =============================================================================
evolver = LevelSetMethodsEvolver(;
    bg_model = model,
    initial_ls = zalesak_disk,
    velocity = velocity,
    spatial_scheme = :WENO5,
    integrator = :RK3,
    bc = :Neumann
)
eg = EvolvingDiscreteGeometry(evolver, model)

info = EvolvingDomains.grid_info(eg)
nx, ny = info.dims
x_range = range(info.origin[1], step=info.spacing[1], length=nx)
y_range = range(info.origin[2], step=info.spacing[2], length=ny)

# Store initial level set for comparison
ϕ_initial = copy(EvolvingDomains.current_levelset(eg))

# =============================================================================
# Initial State
# =============================================================================
fig_initial = plot_levelset(eg; title = "Initial Configuration")
display(fig_initial)

# =============================================================================
# Create Animation Frames (One Complete Revolution)
# =============================================================================
frames_dir = joinpath(@__DIR__, "zalesak_frames")
rm(frames_dir, force=true, recursive=true)
mkpath(frames_dir)

Δt = 0.005
substeps = 4
nframes = 50
reinit_freq = 5

println("\nGenerating frames...")

for frame in 1:nframes
    t = EvolvingDomains.current_time(eg)

    # Create figure
    fig = Figure(size = (700, 700))
    ax = Axis(fig[1, 1],
              title = string("Zalesak Disk (t = ", round(t, digits=3), ", ", round(100*t, digits=0), "% revolution)"),
              xlabel = "x", ylabel = "y",
              aspect = DataAspect())

    ϕ = EvolvingDomains.current_levelset(eg)
    ϕ_2d = reshape(ϕ, (nx, ny))

    heatmap!(ax, x_range, y_range, ϕ_2d; colormap = :RdBu, colorrange = (-0.3, 0.3))
    contour!(ax, x_range, y_range, ϕ_2d; levels = [0.0], linewidth = 4, color = :black)

    # Save frame
    filename = joinpath(frames_dir, string("frame_", lpad(frame, 3, "0"), ".png"))
    save(filename, fig)

    if frame % 10 == 1
        println("  Frame ", frame, "/", nframes, " (t=", round(t, digits=3), ")")
    end

    # Advance
    for _ in 1:substeps
        EvolvingDomains.advance!(eg, Δt)
    end

    # Reinitialize periodically
    if frame % reinit_freq == 0
        EvolvingDomains.reinitialize!(eg)
    end
end

println("✓ Frames saved to: ", frames_dir)

# =============================================================================
# Combine frames into video using FFMPEG.jl
# =============================================================================
output_video = joinpath(@__DIR__, "zalesak_disk.mp4")
input_pattern = joinpath(frames_dir, "frame_%03d.png")
FFMPEG.exe(`-y -framerate 15 -i $input_pattern -c:v libx264 -pix_fmt yuv420p $output_video`)

println("✓ Video saved: ", output_video)

# =============================================================================
# Error Analysis
# =============================================================================
ϕ_final = EvolvingDomains.current_levelset(eg)

error_L2 = sqrt(sum((ϕ_final .- ϕ_initial).^2) / length(ϕ_initial))
error_Linf = maximum(abs.(ϕ_final .- ϕ_initial))

area_initial = sum(ϕ_initial .< 0)
area_final = sum(ϕ_final .< 0)
area_change = (area_final - area_initial) / area_initial * 100

println("\n" * "=" ^ 40)
println("Results after one revolution:")
println("-" ^ 40)
println("  L² error:   $(round(error_L2, sigdigits=4))")
println("  L∞ error:   $(round(error_Linf, sigdigits=4))")
println("  Area change: $(round(area_change, digits=2))%")
println("=" ^ 40)

if error_L2 < 0.01
    println("✓ PASS: L² error within tolerance")
else
    println("⚠ WARNING: L² error higher than expected")
end

# =============================================================================
# Final State
# =============================================================================
println("\nFinal time: ", EvolvingDomains.current_time(eg))
fig_final = plot_levelset(eg; title = "Zalesak Disk (After One Revolution)")
display(fig_final)
