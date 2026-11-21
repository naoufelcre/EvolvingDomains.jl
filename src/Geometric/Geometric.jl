module Geometric

using ..EvolvingDomains: Cut, Field, Grid, _UNIFORM_RTOL, _coordinateTolerance

include("GridInfo.jl")
include("GridInterpolation.jl")
include("Aggregation.jl")

using .GridInterpolation: _getInterpolator, ClampedBilinear, ClampedTrilinear

include("Reinitialization.jl")
using .Reinitialization: reinitialize!

include("Extension.jl")
using .Extension: extend!

include("UnresolvedFeatures.jl")
using .UnresolvedFeatures: filterUnresolvedFeatures!

end # module
