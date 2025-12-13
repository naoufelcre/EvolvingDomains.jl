# =============================================================================
# EvolvingDomains.jl Test Suite
# =============================================================================

using Test
using EvolvingDomains
using Gridap
using GridapEmbedded
using StaticArrays

# Include individual test files
include("test_grid_info.jl")
include("test_external_solver_api.jl")
include("test_colliding_balls.jl")
