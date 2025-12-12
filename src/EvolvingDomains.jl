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

# Backends (loaded conditionally or explicitly)
include("backends/LevelSetMethodsBackend.jl")

# Test visualization (optional helper)
include("TestVisualization.jl")
using .TestVisualization

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

# Exports - Test Visualization (thin wrapper for EmbeddedViz)
export TestVisualization
export setup_test_viz!, show_geometry, close_test_viz
export start_recording!, save_animation

end # module
