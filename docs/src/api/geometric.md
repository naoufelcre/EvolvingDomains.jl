# Geometric Quantities


## User-Facing Functions

```@docs
interface_curvature
interface_normal
curvature_at_band
```

## Cell Classification

```@docs
get_cut_cells
expand_to_band
cells_to_nodes
```

## Visualization

```@docs
plot_curvature
plot_curvature!
```

## Usage

```julia
using EvolvingDomains, Gridap, GridapEmbedded

# Setup
eg = EvolvingDiscreteGeometry(evolver, model)

# Get curvature as CellField (cached, narrow-band computation)
κ_Γ = interface_curvature(eg)

# Get normal from GridapEmbedded
n_Γ = interface_normal(eg)

# Surface tension in weak form
cut_geo = current_cut(eg)
Γ = EmbeddedBoundary(cut_geo)
dΓ = Measure(Γ, 2)
σ = 0.072

# Use in variational form
a(u, v) = ∫( σ * κ_Γ * (n_Γ ⋅ v) )dΓ
```

## Narrow Band Details

Curvature is computed only in a narrow band around the interface:

| Parameter | Default | Description |
|-----------|---------|-------------|
| `n_layers` | 2 | Layers of cells beyond CUT cells |


### Low-Level Access

```julia
# Get curvature values and band info
κ_values, band_nodes = curvature_at_band(eg; n_layers=2)

# κ_values: full-length Vector{Float64} (zeros outside band)
# band_nodes: Set{Int} of node indices where curvature was computed
```
