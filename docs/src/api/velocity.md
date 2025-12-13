# Velocity Sources

## Velocity Source Types

```@docs
AbstractVelocitySource
StaticFunctionVelocity
TimeDependentVelocity
FEVelocitySource
```

## Velocity Extension Strategies

When using FE velocities that are only defined inside the domain, velocity extension strategies determine how to compute velocities at interface nodes.

```@docs
AbstractVelocityExtension
NoExtension
ZeroExtension
ConstantExtension
NarrowBandExtension
```

## Functions

```@docs
sample_velocity
is_time_dependent
update_velocity!
extend_velocity
extend_velocity_narrow_band!
update_levelset!
```
