module EvolvingDomains

using Gridap
using GridapEmbedded
using Gridap.Geometry
using Gridap.CellData
using LinearAlgebra
using StaticArrays

# Abstract interface (backend-agnostic)
include("AbstractEvolver.jl")

# Velocity sources (must be before backends that use them)
include("VelocitySource.jl")

# Main struct
include("EvolvingGeometry.jl")

# Utilities
include("utils.jl")

# Grid info for external solvers
include("GridInfo.jl")

# Simulation snapshotting
include("SimulationSnapshot.jl")

# Geometric quantities (curvature, normals for CutFEM)
include("GeometricQuantities.jl")

# Backends (loaded conditionally or explicitly)
include("backends/LevelSetMethodsBackend.jl")

# Exports - Core
export AbstractLevelSetEvolver
export EvolvingDiscreteGeometry
export advance!, current_geometry, current_time, current_levelset, reinitialize!
export current_cut, invalidate_cache!
export LevelSetMethodsEvolver
export evolve!, current_values  # Low-level evolver interface

# Exports - Velocity Sources
export AbstractVelocitySource
export StaticFunctionVelocity, TimeDependentVelocity, FEVelocitySource
export sample_velocity, is_time_dependent, update_velocity!

# Exports - Velocity Extension
export AbstractVelocityExtension, extend_velocity
export NoExtension, ZeroExtension, ConstantExtension, NarrowBandExtension
export extend_velocity_narrow_band!, update_levelset!

# Exports - Evolver capabilities
export supports_velocity_update

# Exports - External Solver Integration
export CartesianGridInfo, grid_info
export domain_mask, narrow_band_mask
export set_levelset!, set_values!

# Exports - Simulation Snapshotting
export SimulationFrame, SimulationResult, snapshot

# Exports - Geometric Quantities (for CutFEM surface tension)
export interface_curvature, interface_normal, curvature_at_band
export get_cut_cells, expand_to_band, cells_to_nodes
export plot_curvature, plot_curvature!


# =============================================================================
# Visualization Stubs (implemented by EvolvingDomainsMakieExt when GLMakie loaded)
# =============================================================================
"""
    plot_levelset(eg::EvolvingDiscreteGeometry; kwargs...)

Create an interactive contour plot of the current level set.
Requires GLMakie: `using GLMakie` before calling.

See also: [`plot_levelset!`](@ref)
"""
function plot_levelset end

"""
    plot_levelset!(ax, eg::EvolvingDiscreteGeometry; kwargs...)

Add level set visualization to an existing Makie axis.
Requires GLMakie: `using GLMakie` before calling.
"""
function plot_levelset! end

"""
    viewer(result::SimulationResult; wait=false, colormap=:RdBu) -> Figure

Open an interactive viewer with a time slider for cached simulation frames.
Requires GLMakie: `using GLMakie` before calling.

# Arguments
- `result`: Cached simulation frames from `snapshot()` calls
- `wait`: If true, block until window is closed
- `colormap`: Colormap for level set visualization

# Example
```julia
frames = [snapshot(eg) for _ in 1:100 if (advance!(eg, Δt); true)]
result = SimulationResult(grid_info(eg), frames)
viewer(result)
```
"""
function viewer end

"""
    view_live!(result::SimulationResult; fps=30, loop=false)

Play back cached simulation frames as an animation.
Requires GLMakie: `using GLMakie` before calling.

# Arguments
- `result`: Cached simulation frames
- `fps`: Frames per second for playback
- `loop`: If true, loop animation continuously

This function blocks until playback completes (or Ctrl+C).
"""
function view_live! end

export plot_levelset, plot_levelset!, viewer, view_live!

end # module
