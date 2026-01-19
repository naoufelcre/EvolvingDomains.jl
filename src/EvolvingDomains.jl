# =============================================================================
# EvolvingDomains.jl
# =============================================================================
#
# V2 Minimal Rebuild — December 2025
#
# This package provides:
# - Level set evolution using LevelSetMethods.jl (WENO5 + reinitialization)
# - Integration with GridapEmbedded for CutFEM simulations
# - Simple API: advance!, current_cut, set_velocity!
#
# =============================================================================

module EvolvingDomains

# Narrow imports for explicit dependencies
using Gridap: CartesianDiscreteModel, DiscreteModel, FESpace, FEFunction, ReferenceFE, lagrangian, interpolate
using Gridap.Geometry: get_cartesian_descriptor
using GridapEmbedded: cut
using GridapEmbedded.LevelSetCutters: DiscreteGeometry
using StaticArrays: SVector
using LevelSetMethods  # Pinned to #main in Project.toml

# Include modules
include("GridInfo.jl")
include("InterpolationUtils.jl")  # [NEW]
using .InterpolationUtils         # [NEW]



# Velocity Sources (Depends on Quadtree for Guided Lookup)
include("VelocitySource.jl")

include("GeometryDesign.jl")
include("EvolvingGeometry.jl")

# =============================================================================
# Exports — Core API
# =============================================================================

# Geometry Container
export EvolvingDiscreteGeometry
export advance!, current_geometry, current_cut, current_levelset, current_time
export set_levelset!, reinitialize!, set_velocity!

# Grid Info (External Solver Integration)
export CartesianGridInfo, grid_info
export domain_mask, narrow_band_mask

# Velocity Sources
export AbstractVelocitySource
export StaticFunctionVelocity, TimeDependentVelocity, FEVelocitySource, GuidedVelocitySource
export sample_velocity, is_time_dependent, update_velocity!, locate_cell

# Geometry Design (CSG)
export AbstractGeometry, Circle, Rectangle
export Translate
export signed_distance



# =============================================================================
# Visualization Stubs (implemented by extension when CairoMakie is loaded)
# =============================================================================

"""
    plot_levelset(geom::EvolvingDiscreteGeometry; kwargs...)

Create a contour plot of the current level set.
Requires CairoMakie: `using CairoMakie` before calling.
"""
function plot_levelset end

"""
    plot_levelset!(ax, geom::EvolvingDiscreteGeometry; kwargs...)

Add level set visualization to an existing Makie axis.
Requires CairoMakie: `using CairoMakie` before calling.
"""
function plot_levelset! end

export plot_levelset, plot_levelset!

end # module EvolvingDomains
