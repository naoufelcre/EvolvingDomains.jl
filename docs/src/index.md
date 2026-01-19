# EvolvingDomains.jl

**Hybrid Level-Set / CutFEM Simulation Framework.**

`EvolvingDomains.jl` couples grid-based level set evolution (via `LevelSetMethods.jl`) with unfitted Finite Element Methods (via `GridapEmbedded.jl`) to simulate problems on time-dependent domains.

## Installation

```julia
using Pkg
Pkg.add(url="https://github.com/naoufelcre/EvolvingDomains.jl")
```

## Quick Start

See the `examples/` directory in the repository.
