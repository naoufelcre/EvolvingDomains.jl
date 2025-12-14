# Types and Functions

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

## Simulation Snapshotting

Types and functions for caching simulation state for visualization.

### Types

```@docs
SimulationFrame
SimulationResult
```

### Functions

```@docs
snapshot
```

## Visualization

These functions require `GLMakie` to be loaded.

```@docs
plot_levelset
plot_levelset!
viewer
view_live!
```
