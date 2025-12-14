# =============================================================================
# POC: Grid Cell Classification Visualization
# =============================================================================
# Minimal proof-of-concept for cell-based domain visualization.
# Colors cells: white (outside), red (cut by boundary), yellow (inside)
#
# Usage:
#   julia --project=. examples/poc_cell_classification.jl
#
# Expected: Window shows grid with colored cells based on domain status

using Pkg
Pkg.activate(dirname(@__DIR__))

using EvolvingDomains
using Gridap, GridapEmbedded
using LevelSetMethods
using GLMakie

# =============================================================================
# Cell Classification Logic
# =============================================================================
"""
    classify_cells(ϕ_nodes::Matrix) -> Matrix{Int}

Classify each cell based on corner node signs:
  0 = outside (all ϕ > 0)  → white
  1 = cut     (mixed signs) → red  
  2 = inside  (all ϕ < 0)  → yellow
"""
function classify_cells(ϕ_nodes::Matrix{Float64})
    nx, ny = size(ϕ_nodes)
    # Cells are between nodes: (nx-1) × (ny-1) cells
    cell_class = zeros(Int, nx-1, ny-1)
    
    @inbounds for j in 1:ny-1, i in 1:nx-1
        # Get 4 corner values
        c1 = ϕ_nodes[i,   j]
        c2 = ϕ_nodes[i+1, j]
        c3 = ϕ_nodes[i,   j+1]
        c4 = ϕ_nodes[i+1, j+1]
        
        # Count how many are inside (ϕ < 0)
        n_inside = (c1 < 0) + (c2 < 0) + (c3 < 0) + (c4 < 0)
        
        cell_class[i,j] = if n_inside == 0
            0  # All outside
        elseif n_inside == 4
            2  # All inside
        else
            1  # Cut cell
        end
    end
    return cell_class
end

# =============================================================================
# Setup
# =============================================================================
domain = (-1.0, 1.0, -1.0, 1.0)
partition = (30, 30)  # Coarse grid to see individual cells
model = CartesianDiscreteModel(domain, partition)

# Circle (half inside, half outside for interesting cut cells)
R = 0.5
ϕ₀(x) = R - sqrt(x[1]^2 + x[2]^2)

evolver = LevelSetMethodsEvolver(;
    bg_model = model,
    initial_ls = ϕ₀,
    velocity = (0.0, 0.0),  # Static for clarity
    spatial_scheme = :WENO5,
    bc = :Neumann
)
eg = EvolvingDiscreteGeometry(evolver, model)

# =============================================================================
# POC: Cell Classification Plot
# =============================================================================
println("=" ^ 50)
println("POC: Grid Cell Classification")
println("=" ^ 50)

info = EvolvingDomains.grid_info(eg)
nx, ny = info.dims
ncx, ncy = info.cells

# Node coordinates (for contour overlay)
x_nodes = range(info.origin[1], step=info.spacing[1], length=nx)
y_nodes = range(info.origin[2], step=info.spacing[2], length=ny)

# Cell center coordinates (for heatmap)
x_cells = range(info.origin[1] + info.spacing[1]/2, step=info.spacing[1], length=ncx)
y_cells = range(info.origin[2] + info.spacing[2]/2, step=info.spacing[2], length=ncy)

# Get level set and classify
ϕ = EvolvingDomains.current_levelset(eg)
ϕ_2d = reshape(ϕ, (nx, ny))
cell_class = classify_cells(ϕ_2d)

# Count cells by type
n_outside = count(==(0), cell_class)
n_cut = count(==(1), cell_class)
n_inside = count(==(2), cell_class)
println("Cell counts: outside=$n_outside, cut=$n_cut, inside=$n_inside")

# Create figure
fig = Figure(size = (700, 600))
ax = Axis(fig[1, 1],
          title = "Grid Cell Classification",
          xlabel = "x", ylabel = "y",
          aspect = DataAspect())

# Custom discrete colormap: white (0), red (1), yellow (2)
colors = [:white, :red, :yellow]
cmap = cgrad(colors, 3, categorical=true)

# Plot cell classification as heatmap
hm = heatmap!(ax, x_cells, y_cells, cell_class;
              colormap = cmap, 
              colorrange = (-0.5, 2.5))  # Centers discrete values

# Overlay zero contour (domain boundary)
contour!(ax, x_nodes, y_nodes, ϕ_2d;
         levels = [0.0], linewidth = 3, color = :black)

# Optional: show grid lines
for x in x_nodes
    vlines!(ax, x; color = :gray70, linewidth = 0.5)
end
for y in y_nodes
    hlines!(ax, y; color = :gray70, linewidth = 0.5)
end

# Legend
Legend(fig[1, 2],
       [PolyElement(color = :white, strokecolor = :black),
        PolyElement(color = :red),
        PolyElement(color = :yellow)],
       ["Outside (ϕ > 0)", "Cut (boundary)", "Inside (ϕ < 0)"])

display(fig)

println("=" ^ 50)
println("POC Complete - check the displayed figure")
println("=" ^ 50)

# Keep window open
println("Press Enter to close...")
readline()
