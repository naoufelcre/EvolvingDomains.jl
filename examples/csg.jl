using EvolvingDomains
using Gridap
using CairoMakie
using LevelSetMethods

"""
    tutorial_01_csg.jl

This tutorial demonstrates the Constructive Solid Geometry (CSG) capabilities
of EvolvingDomains.jl. We will construct complex shapes by combining
simple primitives like Circles and Rectangles.

Key Concepts:
- Primitives: `Circle`, `Rectangle`
- Operations: `union` (∪), `intersect` (∩), `setdiff` (setminus)
- Visualization via `plot_levelset`
"""

# 1. Setup Grid
domain = (-2.0, 2.0, -2.0, 2.0)
n = 100
model = CartesianDiscreteModel(domain, (n, n))

# 2. Define Primitives
# A circle centered at origin
c = EvolvingDomains.Circle(VectorValue(0.0, 0.0), 1.0)

# A rectangle to cut out a slot (Zalesak style)
# Slot width = 0.2, Depth = 1.0
r_slot = EvolvingDomains.Rectangle(VectorValue(-0.1, 0.0), VectorValue(0.1, 1.2))

# A box for intersection (clipping)
r_clip = EvolvingDomains.Rectangle(VectorValue(-0.8, -0.8), VectorValue(0.8, 0.8))

# 3. CSG Operations
# Example A: The "Pacman" / Zalesak Disk (Circle minus Slot)
geo_zalesak = setdiff(c, r_slot)

# Example B: A clipped integration (Circle intersection Box)
geo_clipped = intersect(c, r_clip)

# Example C: Complex Union (Two circles)
c2 = EvolvingDomains.Circle(VectorValue(1.0, 1.0), 0.5)
geo_union = union(c, c2)

# 4. Initialize Geometry Objects
# We wrap them in EvolvingDiscreteGeometry to use the visualization tools
# Note: We use static geometries here (velocity=nothing)
# CRITICAL FIX: Ensure x is converted to VectorValue for geometry functions
geom_zalesak = EvolvingDiscreteGeometry(model, x -> geo_zalesak(VectorValue(x...)))
geom_clipped = EvolvingDiscreteGeometry(model, x -> geo_clipped(VectorValue(x...)))
geom_union   = EvolvingDiscreteGeometry(model, x -> geo_union(VectorValue(x...)))

# 5. Visualize
println("Generating visualization...")

fig = Figure(size=(1200, 400))

# Plot A: Zalesak
ax1 = Axis(fig[1, 1], title="A. Difference (Zalesak)", aspect=1)
plot_levelset!(ax1, geom_zalesak)

# Plot B: Clipped
ax2 = Axis(fig[1, 2], title="B. Intersection (Clipped)", aspect=1)
plot_levelset!(ax2, geom_clipped)

# Plot C: Union
ax3 = Axis(fig[1, 3], title="C. Union (Metaballs)", aspect=1)
plot_levelset!(ax3, geom_union)

output_path = joinpath(@__DIR__, "csg_overview.png")
save(output_path, fig)
println("Saved $output_path")

#md # ![CSG Overview](csg_overview.png)
