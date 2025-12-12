# Core Types and Functions

The core API for evolving domain simulations.

## Main Types

```@docs
EvolvingDiscreteGeometry
AbstractLevelSetEvolver
LevelSetMethodsEvolver
```

## Core Functions

### Geometry Evolution

```@docs
advance!
reinitialize!
invalidate_cache!
```

### Accessing Current State

```@docs
current_geometry
current_cut
current_time
current_levelset
```

### Evolver Interface

These are the methods that any level set evolver backend must implement:

```@docs
evolve!
current_values
supports_velocity_update
```
