module Geometric

using ..EvolvingDomains: Cut, Field, Grid, _UNIFORM_RTOL, _coordinateTolerance

include("GridInfo.jl")
include("GridInterpolation.jl")
include("Aggregation.jl")

using .GridInterpolation: _getInterpolator, ClampedBilinear, ClampedTrilinear

include("Reinitialization.jl")
using .Reinitialization: reinitialize!

include("ClosestPointExtension.jl")
using .ClosestPointExtension: extendClosestPoint!

include("ExtendCutValues.jl")
using .ExtendCutValues: extendCutValues, extendCutValues!

include("TopologyFilter.jl")
using .TopologyFilter: filterSmallPhaseIslands!, filterSubResolutionFeatures!

end # module
