module EvolvingDomains

# External dependencies (re-exported or used by submodules)
using Gridap: CartesianDiscreteModel, DiscreteModel, FESpace, FEFunction, ReferenceFE, lagrangian, interpolate
using Gridap.Geometry: get_cartesian_descriptor, get_node_coordinates
using GridapEmbedded: cut
using GridapEmbedded.LevelSetCutters: DiscreteGeometry
using StaticArrays: SVector
using Interpolations: Interpolations

# Submodules
include("Geometric/Geometric.jl")
using .Geometric

include("Kinematic/Kinematic.jl")
using .Kinematic

include("Transfer/Transfer.jl")
using .Transfer

# =============================================================================
# Exports

# Kinematic
export AbstractVelocitySource
export StaticFunctionVelocity, TimeDependentVelocity
export sample_velocity, is_time_dependent
export advance!

# Geometric
export EvolvingDiscreteGeometry
export current_cut, current_levelset, ensure_cut!
export set_levelset!, reinitialize!, tangential_smooth!
export InterfaceSamples, interface_samples, interface_curvature, get_curvature, curvature_at
export CartesianGridInfo, grid_info
export AbstractGeometry, Circle, Rectangle, Translate, signed_distance
export WENO5Cache

# Transfer
export GridMeshTransfer, setup_transfer, get_transfer_op, update_transfer_cache!
export prolong, restrict
export ClosestPointExtension, get_extension_op, extend, update_extension_cache!

# Visualization Stubs
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
