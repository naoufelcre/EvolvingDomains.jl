# =============================================================================
# EvolvingDomains.jl
# =============================================================================


module EvolvingDomains

# Narrow imports for explicit dependencies
using Gridap: CartesianDiscreteModel, DiscreteModel, FESpace, FEFunction, ReferenceFE, lagrangian, interpolate
using Gridap.Geometry: get_cartesian_descriptor
using GridapEmbedded: cut
using GridapEmbedded.LevelSetCutters: DiscreteGeometry
using StaticArrays: SVector
using LevelSetMethods  # Pinned to main because we are not using a registered version
using Interpolations: Interpolations
using NearestNeighbors

# Include modules
include("GridInfo.jl")
include("InterpolationUtils.jl")
using .InterpolationUtils



# Velocity Sources
include("VelocitySource.jl")

# Transfer Operator
include("GridMeshTransfer.jl")

# Extension Operator
include("ExtensionOperator.jl")

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
export StaticFunctionVelocity, TimeDependentVelocity
export sample_velocity, is_time_dependent

# Transfer
export GridMeshTransfer, setup_transfer, get_transfer_op, update_transfer_cache!
export prolong, restrict

# Extension
export ClosestPointExtension, get_extension_op, extend, update_extension_cache!

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
