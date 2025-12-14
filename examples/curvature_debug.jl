# Curvature Visualization Debug Example
# =====================================
# Demonstrates the narrow-band curvature computation and visualization

using Pkg
Pkg.activate(dirname(@__DIR__))

using EvolvingDomains
using Gridap, GridapEmbedded
using LevelSetMethods
using GLMakie

# =============================================================================
# Setup: Circle with known curvature κ = 1/R
# =============================================================================
println("=" ^ 50)
println("Curvature Visualization Debug Example")
println("=" ^ 50)

domain = (0.0, 1.0, 0.0, 1.0)
partition = (60, 60)
model = CartesianDiscreteModel(domain, partition)

R = 0.25
center = (0.5, 0.5)
ϕ(x) = sqrt((x[1]-center[1])^2 + (x[2]-center[2])^2) - R

println("\nCircle Parameters:")
println("  Radius R = $R")
println("  Expected curvature κ = 1/R = $(round(1/R, digits=3))")

evolver = LevelSetMethodsEvolver(;
    bg_model = model,
    initial_ls = ϕ,
    velocity = x -> (0.0, 0.0),
    spatial_scheme = :WENO5,
    bc = :Neumann
)
eg = EvolvingDiscreteGeometry(evolver, model)

# =============================================================================
# Compute and Display Curvature Statistics
# =============================================================================
κ_values, band_nodes = curvature_at_band(eg)
info = grid_info(eg)
ϕ_vals = current_levelset(eg)

# Near-interface curvature
near_interface = findall(i -> abs(ϕ_vals[i]) < 0.03 && κ_values[i] != 0, 1:length(κ_values))
κ_near = κ_values[near_interface]

println("\n--- Narrow Band Statistics ---")
println("  Total nodes: $(length(κ_values))")
println("  Nodes in band: $(length(band_nodes)) ($(round(100*length(band_nodes)/length(κ_values), digits=1))%)")
println("  Nodes near interface: $(length(near_interface))")

println("\n--- Curvature Accuracy ---")
println("  Mean κ at interface: $(round(mean(κ_near), digits=4))")
println("  Expected κ = 1/R: $(round(1/R, digits=4))")
println("  Relative error: $(round(100*abs(mean(κ_near) - 1/R)/(1/R), digits=2))%")

# =============================================================================
# Side-by-Side Visualization: Level Set vs Curvature
# =============================================================================
println("\n--- Creating Visualization ---")

fig = Figure(size = (1200, 500))

# Left: Level Set
ax1 = Axis(fig[1, 1], 
           title = "Level Set ϕ",
           xlabel = "x", ylabel = "y",
           aspect = DataAspect())

x = range(info.origin[1], step=info.spacing[1], length=info.dims[1])
y = range(info.origin[2], step=info.spacing[2], length=info.dims[2])
ϕ_2d = reshape(ϕ_vals, info.dims)
κ_2d = reshape(κ_values, info.dims)

hm1 = heatmap!(ax1, x, y, ϕ_2d; colormap = :RdBu, colorrange = (-0.3, 0.3))
contour!(ax1, x, y, ϕ_2d; levels = [0.0], linewidth = 3, color = :black)
Colorbar(fig[1, 2], hm1, label = "ϕ")

# Right: Curvature in Narrow Band
ax2 = Axis(fig[1, 3], 
           title = "Curvature κ (narrow band only)\nExpected: κ = $(round(1/R, digits=2))",
           xlabel = "x", ylabel = "y",
           aspect = DataAspect())

κ_max = maximum(abs, κ_near)
hm2 = heatmap!(ax2, x, y, κ_2d; colormap = :viridis, colorrange = (0, κ_max * 1.2))
contour!(ax2, x, y, ϕ_2d; levels = [0.0], linewidth = 3, color = :white)
Colorbar(fig[1, 4], hm2, label = "κ")

# Save figure
output_path = joinpath(@__DIR__, "curvature_debug.png")
save(output_path, fig)
println("\n✓ Saved: $output_path")

display(fig)
println("\nDone!")
