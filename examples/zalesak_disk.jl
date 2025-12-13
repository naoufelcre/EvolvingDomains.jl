# Zalesak's Rotating Disk
# =======================

using EvolvingDomains
using Gridap
using GridapEmbedded
using LevelSetMethods

# =============================================================================
# Domain Setup
# =============================================================================
# Centered domain for symmetric rotation about origin
domain = (-0.5, 0.5, -0.5, 0.5)
partition = (100, 100)
model = CartesianDiscreteModel(domain, partition)

println("Zalesak's Rotating Disk Benchmark")
println("=" ^ 40)
println("Grid: $(partition[1])×$(partition[2])")

# =============================================================================
# Zalesak Disk Level Set
# =============================================================================

function zalesak_disk(x)
    # Disk parameters
    cx, cy = 0.0, 0.25          # Disk center
    radius = 0.15               # Disk radius

    # Slot parameters
    slot_width = 0.05           # Slot width
    slot_height = 0.25          # Slot height (extends into disk)

    # Signed distance to circle (negative inside)
    d_circle = sqrt((x[1] - cx)^2 + (x[2] - cy)^2) - radius

    # Slot bounds
    xmin = cx - slot_width/2
    xmax = cx + slot_width/2
    ymin = cy - radius
    ymax = cy - radius + slot_height

    # Signed distance to slot rectangle
    dx = max(xmin - x[1], x[1] - xmax, 0.0)
    dy = max(ymin - x[2], x[2] - ymax, 0.0)
    inside_slot = (xmin ≤ x[1] ≤ xmax) && (ymin ≤ x[2] ≤ ymax)

    if inside_slot
        # Inside slot: distance to nearest slot wall (negative)
        d_slot = -min(x[1] - xmin, xmax - x[1], x[2] - ymin, ymax - x[2])
    else
        # Outside slot: distance to slot boundary
        d_slot = sqrt(dx^2 + dy^2)
    end

    # CSG: disk minus slot = max(-d_circle, d_slot)
    # We want: inside disk AND outside slot
    return max(-d_circle, d_slot)
end

# =============================================================================
# Rigid Body Rotation Velocity
# =============================================================================

ω = 2π  # Angular velocity
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

# Store initial level set for comparison
ϕ_initial = copy(current_levelset(eg))

# =============================================================================
# Time Evolution: One Complete Revolution
# =============================================================================
Δt = 0.005                   # Small time step for accuracy
T_period = 1.0               # One revolution
nsteps = Int(T_period / Δt)
reinit_freq = 20

println("\nTime stepping:")
println("  Δt = $Δt")
println("  Steps: $nsteps")
println("-" ^ 40)

for step in 1:nsteps
    advance!(eg, Δt)

    if step % reinit_freq == 0
        reinitialize!(eg)
    end

    # Progress report every 20%
    if step % (nsteps ÷ 5) == 0
        t = EvolvingDomains.current_time(eg)
        revolutions = t / T_period
        println("  t = $(round(t, digits=3)) ($(round(100*revolutions, digits=0))% revolution)")
    end
end

# =============================================================================
# Error Analysis
# =============================================================================
ϕ_final = current_levelset(eg)

# L2 error in level set values
error_L2 = sqrt(sum((ϕ_final .- ϕ_initial).^2) / length(ϕ_initial))

# L∞ error
error_Linf = maximum(abs.(ϕ_final .- ϕ_initial))

# Mass (area) conservation: count nodes inside
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
