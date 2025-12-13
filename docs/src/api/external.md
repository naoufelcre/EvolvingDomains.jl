# External Solver Interface

## Grid Information

```@docs
CartesianGridInfo
grid_info
```

## Domain Masks

These functions help identify which grid nodes are inside the domain or near the interface:

```@docs
domain_mask
narrow_band_mask
```

## Level Set Injection

For external solvers that evolve the level set independently:

```@docs
set_levelset!
set_values!
```

## Typical Usage Pattern

```julia
# Get grid metadata
info = grid_info(eg)
Δx, Δy = info.spacing
nx, ny = info.dims

# Get current level set and masks
ϕ = current_levelset(eg)
inside = domain_mask(ϕ)
near_interface = narrow_band_mask(ϕ, 6 * Δx)

# External solver computes new level set
ϕ_new = my_external_solver(info, ϕ, Δt)

# Inject back into EvolvingDomains
set_levelset!(eg, ϕ_new)
reinitialize!(eg)  # Restore signed distance property

# Continue with CutFEM
cut_geo = current_cut(eg)
```
