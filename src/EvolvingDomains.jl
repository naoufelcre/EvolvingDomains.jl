module EvolvingDomains

# Package-owned data model and Cut snapshots (included before operators).
include("DataModel.jl")
include("CutData.jl")

# Submodules. Imports are explicit: the public surface is declared once, below.
include("Geometric/Geometric.jl")
using .Geometric: Geometric, aggregateCutCells, reinitialize!, extend!,
                  filterUnresolvedFeatures!

include("Kinematic/Kinematic.jl")
using .Kinematic: Kinematic, eulerianTransport!, lagrangianTransport!, TransportMap,
                  CIPCache, WENO5Cache

# Multilevel transfer (topology lives with the data model in DataModel.jl).
include("Transfer.jl")

# Private native CutCells implementation. It returns the owned flat Cut data structure
# directly; there is no engine abstraction or selector.
include("Bindings/CutCells.jl")

"""
    cut(field::Field, level_id::Symbol)

Return a `Cut` from a field and a given level id: the nodal values attached to that level are
used to build a cut of their zero level set. They must be finite, on a uniformly spaced 2D or 3D
grid.
"""
cut(field::Field, level_id::Symbol) = CutCellsBinding._fromField(field, level_id)

# =============================================================================
# Exports: the single declaration site of the public surface.
#
# Every name here is also reachable qualified (EvolvingDomains.Kinematic.weno5Step!,
# EvolvingDomains.Kinematic.CIP, …). Anything not listed is internal and may change
# without a breaking release.

# Data model
export Grid, Discretization, Field

# Multilevel: topology, data wiring, navigation, algebra
export refine!, coarsen!, attach!, detach!
export immediateCoarser, immediateFiner, coarsest, finest
export transfer!, linearInterpolator

# Cuts: the owned snapshot, its interface surface and cell-state codes, and the producer
export Cut, cut, CutSurface, CUT_IN, CUT_OUT, CUT_CUT, subcells, cutCellIds

# Kinematics
export eulerianTransport!, lagrangianTransport!
export TransportMap, CIPCache, WENO5Cache

# Geometry
export aggregateCutCells, reinitialize!, extend!
export filterUnresolvedFeatures!

end # module EvolvingDomains
