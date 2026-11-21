# EvolvingDomains

Julia library for composable and minimal orchestration of simulation components for writing evolving domain problems.

Its design aim for minimality and composability in the following sense:
- **Minimalism**: The package public API is centered around 3 structs that are basically just organized arrays. By doing so, we aim to have mathematical structures of package component clearly exposed.
- **Composable**: Because all that the package manipulate is bare-bones arrays, it's relatively easy to use the package as the orchestrator between different solvers. We provide a shared representation, but you can (and should) compose with methods of your choice.

## Overview

Let $d=2$ or $3$ and $\mathcal{O} \subset \mathbb{R}^d$. Our aim is to solve problems of the form

$$
\frac{d}{d\,t}S(t) = \mathcal{G}(t,S(t)),
$$
where for all $t$ in $[0,T]$: $S(t) = (\Omega(t),s(t))$ with $\Omega(t) \subset \mathcal{O}$ is the **evolving domain** and $s(t,\mathbf{x}) \in \mathcal{B}(\Omega(t))$ is a **state function**. We assume that $\Omega(t)$ can be described as the sub-zero level-set of a function $\phi:\mathcal{O}\to \mathbb{R}$ : $\Omega(t) = \{ \mathbf{x} \in \mathcal{O}, \phi(t,\mathbf{x}) < 0 \}$.

### API overview

The library is organized around the following three types:

```julia
Grid(axes::AbstractVector...)          # or Grid((v1, v2, ...)); Cartesian, strictly increasing axes
    .axes::NTuple{D,Vector{T}}   .mesh_id::UInt64    # fresh id per construction, survives copy

Discretization(axes::AbstractVector...)        # shorthand for Dict(:init => Grid(axes))
Discretization(levels::AbstractDict)           # :Symbol => Grid; all levels share one physical box
    .levels::Dict{Symbol,Grid{D,T}}

Field(d::Discretization, data::AbstractDict)   # :Symbol => dense Array; retained, not copied
    .discretization   .data::Dict{Symbol,Array{T,N}}
```
Schematically, those are nested structures.

### Typical structure of a program

A program written using `EvolvingDomains.jl` take the following form:

```julia
# Initialization
geometry = Field(grid, initial_levelset)
state    = Field(grid, initial_state)

config = (; geometry, state, parameters, time = 0.0)

function dynamics(config)
    # May involve external PDE solvers, geometric operators,
    # constitutive laws, constraints, or prescribed forcing.
    return rates
end

# Explicit simulation loop
while config.time < final_time

    rates = dynamics(config)
    velocity = velocityFromRates(rates, config)

    geometry = eulerianTransport(config.geometry, velocity, dt)
    
    state = lagrangianTransport(config.state, velocity, dt)  # optional
    state = evolve_state(state, geometry, rates, dt)        # optional
    
    config = update(config; geometry, state, time = config.time + dt)
end
```
In particular, we think verifiability of LLMs generated code is critical. To this end we propose a skill (see: https://agentskills.io/home) in `SKILL.md` for enforcing LLMs to comply to this target program structure.

## Composable 

Interfaces with the following libraries were tested: 

Trixi: for eulerianTransport and lagrangianTransport
Ferrite: for FEM solving
Gridap ecosystem: for FEM solving
FreeFem++: for FEM solving

## Examples



-------------------------------------------------------------- DO NOT EDIT WHAT IS ABOVE !!!!

Let  and $f(t,\mathbf{x}): [0,\infty) \times \mathcal{A} \subseteq \mathbb{R}^d $ be a scalar valued function defined over an ambiant space.  


## Core model

```julia
using EvolvingDomains

x = range(0.0, 1.0; length=65)
g = Grid((x, x))
d = Discretization(Dict(:init => g))
field = Field(d, Dict(:init => [hypot(a-0.5, b-0.5) - 0.25 for a in x, b in x]))
```

- `Grid` owns Cartesian axes and a unique mesh identity.
- `Discretization` maps stable `Symbol` levels to Grids.
- `Field` retains dense value arrays per level.

```julia
d.levels[:init]
field.data[:init]
field.discretization
```

## Cuts

`cut` produces an owned, backend-neutral `Cut` snapshot using the private native CutCells backend.

```julia
snapshot = cut(field, :init)
snapshot.bg_states
snapshot.subcell_states
snapshot.boundary
```

The classification is indexed, not just stored: `cut_cell_ids(snapshot)` lists the `CUT_CUT`
cells and `subcells(snapshot, cell)` returns the subcell ids that cell owns (a range lookup
over a CSR index built once at construction), so consumers need not re-filter
`cell_to_bgcell`. Subcells are grouped by ascending background cell in every `Cut`, which is
what makes those ids a plain range — no id array is stored.

```julia
for cell in cut_cell_ids(snapshot), subcell in subcells(snapshot, cell)
    snapshot.subcell_states[subcell] == CUT_IN || continue   # the in-phase piece of `cell`
    q = Int(snapshot.point_pointers[subcell])
    snapshot.physical_points[Int.(snapshot.point_ids[q:(q + 2)])]   # its triangle
end
```

There is no engine abstraction, selector, or public engine module. The native backend is
private; any compatible cutter can construct the same flat `Cut` fields directly.

A `Cut` owns its arrays — the checking constructor copies what it is given, the native
producer adopts the arrays it has just built — and retains neither the Grid nor the source
Field. Interface geometry is exposed as `CutSurface`, whose guarantees are stated in its
docstring. Two- and three-dimensional Cartesian cuts are supported on uniformly spaced axes;
the native backend is reached through a stable C ABI.

## Field evolution

WENO5 + SSP-RK3 for eulerian advection of a nodal field:

```julia
using StaticArrays: SVector
velocity = fill(SVector(0.1, -0.05), length(field.data[:init]))
eulerianTransport!(field, velocity, 1e-3; level_id=:init, cache=WENO5Cache())
```

Intensive CIP is available through `lagrangianTransport!`; conservative Cut-to-Cut transport
uses `TransportMap`. Its ordering is the load-bearing part and no type enforces it: the map
needs the cut *before* the level set moves and the cut *after*.

```julia
# `level`, `g`, `velocity` and `dt` from the examples above; two buffers, because source and
# target must not alias.
source, target = field, other
cache = WENO5Cache()
for step in 1:steps
    previous = cut(level, :init)                    # geometry at t⁻
    eulerianTransport!(level, velocity, dt; level_id=:init, cache=cache)
    current  = cut(level, :init)                    # geometry at t⁺
    map = TransportMap(previous, current, g, velocity, dt; cache_weights=true)
    lagrangianTransport!(target, source, map)
    source, target = target, source                  # source and target must not alias
    reinitialize!(level, current)
end
```

The map takes the `Grid`, not the `Field`: it derives node coordinates from the mesh and
never reads a level set.

## Multilevel transfer

Levels are ordinary `Symbol`s in the `Discretization`; `refine!`, `coarsen!`, `attach!`,
`detach!` and `transfer!` move data between them:

```julia
fine = refine!(d, :init; factor=2, name=:fine)
attach!(field, fine, zeros(129, 129))
transfer!(field, :init, fine)
```

`transfer!` interpolates with `linearInterpolator` by default, works on any grid, and
requires the two levels to be nested (one axis set a subsampling of the other). Fine→coarse
is point sampling, not conservative restriction.

Transfer is interpolation between function spaces, so it **cannot create information**.
Measured from a 9-node start refined twice to 33 nodes: volume 12% below the analytic value
and the interface 0.59 cells off, none of which decays during later advection because the
error was never in the transport. The two orders are not symmetric: coarse→fine→coarse
reproduces the coarse values bitwise, while fine→coarse→fine loses the prolongation error
(measured 9.2e-3 in 2D and 1.6e-2 in 3D at the resolutions tried). Zero drift in one round trip is therefore not evidence of a
lossless transfer.

## Conventions

- **Sign.** `φ < 0` is the material phase; `φ > 0` is outside. `CUT_IN` / `CUT_OUT` /
  `CUT_CUT` name the three background-cell states produced by `cut`.
- **Conservative transport.** `lagrangianTransport!` with a `TransportMap` conserves the
  nodal sum exactly — summing a source's contributions over all targets gives `1 ± 4e-16`,
  the invariant of Lentine, Grétarsson & Fedkiw 2011, §2 — and reproduces a constant field
  to rounding (measured 2e-16 at 17², 1.3e-15 at 65²) from six nodes inward. Within a boundary layer at the box edge and at the support
  edge the operator is no longer a partition of unity; totals stay quantitative, boundary
  values do not. See *Known limitations*.
- **Velocity.** Both operators want a nodal vector field in **node order with the first axis
  fastest**: `vec` of an `(nx, ny[, nz])` array of `SVector`s gives that order.
  `eulerianTransport!` requires exactly `Vector{SVector{N,Float64}}`, so a matrix, a view or
  `Float32` nodes is a `MethodError` there; `TransportMap` accepts any
  `AbstractVector{<:SVector}` and materializes it.
- **Boundary.** The WENO5 kernel clamps its stencil at the domain edge (zero-gradient
  outflow): nothing prevents mass leaving a shape that reaches the boundary, and no inflow
  condition is applied. CIP instead refuses a departure point outside the domain.
- **Resolution.** WENO5 is a five-point reconstruction and CIP a cubic Hermite profile, so
  they need at least five and three nodes per axis respectively. Shorter meshes are accepted
  (stencil indices clamp) but no longer represent the advertised method.

## Choosing the timestep

The operators do not choose a timestep and give no stability diagnostic: that is a property
of the pair (scheme, data) and belongs to the caller. They do refuse a step they cannot
represent — CIP raises when a departure point leaves the domain, `TransportMap` when a
characteristic misses the band — but neither is a CFL test. The estimate is the usual one,
and it is enough to measure it once per run:

```julia
using LinearAlgebra: norm
h  = minimum(minimum(diff(axis)) for axis in g.axes)       # smallest node spacing
vmax = maximum(norm, velocity)                            # frozen step velocity
dt = 0.4 * h / vmax                                       # CFL 0.4
steps = ceil(Int, t_end / dt)                              # t_end: the horizon you want
```

Measured on a rigid rotation of a disc (65², WENO5 + SSP-RK3), for the level-set update:

| `dt·|v|/h` | outcome |
|---|---|
| 1.0 | stable, bounded |
| 1.9 | still bounded |
| 2.8 | `max|φ|` reaches ~3e3 within 8 steps |
| 14 | reaches 3e20 within 8 steps |

Collapse is how it fails: the level set stops being a level set, the band no longer covers
the characteristics, and `TransportMap` then refuses to build a map for that step — its
error names the node and says so. If that message appears, check `dt` first, then whether
the band has thinned below the stencil.

## Geometry operators

```julia
roots = aggregateCutCells(snapshot, g)
reinitialize!(field, snapshot)
extend!(dest, values, snapshot, g)             # per-facet values -> nodes
extend!(dest, nodal, snapshot, g)              # nodal values -> closest-point interp.
filterUnresolvedFeatures!(levelset, level)
```

2D and 3D implementations exist for aggregation, reinitialization, and conservative
transport. The 3D cut uses a six-tetrahedron piecewise-linear model.

`extend!` is the single extension primitive: it projects grid nodes onto the reconstructed
interface and evaluates the interface data there, `u_ext(x) = u(P(x))` with
`P(x) = argmin_{p ∈ Γ_h} |x − p|` taken on the actual segments or triangles of
`snapshot.boundary`, never on their centroids. It needs no level set and no gradient, so it
works on any cut and is insensitive to the level-set scaling.

The source selects the representation by dispatch. A `Vector` with one entry per facet gives
the piecewise-constant extension, `u_ext(x) = values[f*(x)]`, where `values[f]` belongs to
facet `f` in the snapshot's own order (`snapshot.boundary.facet_to_bgcell[f]`): per-facet data
such as a curvature-weighted normal needs no mapping. A grid-shaped array of nodal values
(`Array{Float64,N}`, or `SVector`-valued) is interpolated at `P(x)` instead, which requires
the nodal data to carry a valid trace at the interface — interpolation across it mixes both
sides and does not reconstruct a one-sided trace. The `Field` method
`extend!(dest, source, snapshot)` works on the attached level and returns `dest`.

`band` (default three node spacings) is the distance around the interface where the extension
is written; nodes farther away, and every node when the cut has no facets, keep their
previous value, so pre-fill `dest` if you need a sentinel there. Equal distances go to the
lowest facet index, the cost is `O(nfacets × band nodes)`, and the projection point is
discontinuous across the medial axis by construction. The continuous target is the
normal-constant extension `∇(Eu)·∇d = 0`; on a piecewise-linear interface the max-norm error
is measured first order in `h` because the projection point is quantised by the facet width.

`filterUnresolvedFeatures!` removes interface structures the grid cannot resolve, following
the scale-separation method of Luo, Hu & Adams (2016). With `eps = 0.75h`, a cell cut by
`φ = 0` is *unresolved* when its 3-neighbourhood lacks one of the auxiliary level sets
`φ = ±eps`; its values are then re-distanced from the auxiliary level set that is still
visible within a 9-neighbourhood. A structure that no auxiliary level set reaches — an
isolated island whose phase never exceeds `eps` — is replaced by samples of the surrounding
phase; that fallback is this package's own addition and never deletes a filament attached to
resolved geometry, which is left to the published reconstruction. The method is local: a
neck only a few cells long may keep enough auxiliary cells in its neighbourhood to be judged
resolved.

It is **field hygiene, not an accuracy improvement**: it cannot tell numerical junk from a
genuine sub-resolution feature, it breaks the signed-distance property locally, and it
conserves neither area nor volume. It requires an approximately signed-distance input and
performs no reinitialization itself, so the call stays explicit and a fresh snapshot must be
built afterwards:

```julia
eulerianTransport!(field, velocity, dt; level_id=:init)
filterUnresolvedFeatures!(field, :init)
snapshot = cut(field, :init)          # any older snapshot describes the old geometry
```


`aggregateCutCells` returns one root per background cell, and **`0` means "no aggregate"**:
the cell is outside the phase. The aggregate set is `unique(r for r in roots if r > 0)`. A
cut cell is never left unassigned (the call errors instead), but a field with no zero
crossing at all has nothing to group, so every root is 0 and the call warns. The sentinel
follows the GridapEmbedded convention (`cell_to_cellin::Vector{Int32}` with 0 = unassigned),
but it is not a drop-in: that aggregator asserts every cut cell was assigned and aggregates
the opposite phase too, while here every `CUT_OUT` cell is 0.

Interface curvature and tangents are no longer part of this package: they moved to a
standalone package that consumes `Cut` and `Field`. The `Cut` snapshot still carries the
interface facets and normals they need.

## Scope

The transport operators are **reference implementations, not the state of the art**.
`eulerianTransport!` (WENO5 + SSP-RK3) is a standard level-set advection, and the
conservative Cut-to-Cut scheme redistributes mass exactly while smoothing the profile —
inherent to semi-Lagrangian methods. Measured on one full rotation of a disc (exact flow map
= identity): the smoothing is about one cell per step, so how much a feature loses is set by
how many cells resolve it. Measured at 65², CFL 0.4, peak after one revolution:
`σ/h = 3.2` → 0.04, `σ/h = 6.4` → 0.15, `σ/h = 12.8` → 0.57, with the total exact to 1e-15
in every case. Sub-cell features are destroyed within a turn; only features spanning tens of
cells keep their amplitude.
A user who needs a less diffusive or higher-order scheme is expected to write it against the
data model (`Grid`, `Field`, `Cut`, `CutSurface`); that stable contract, not the shipped
kernels, is the deliverable.

Implemented numerical methods currently require uniform 2D or 3D Cartesian nodal
`Float64` arrays, except multilevel `transfer!`, which only needs nested axes. Generic
Grid/Field storage is dimension-general, but that does not imply generic numerical
operators. Nonuniform evolution and welded/indexed 3D downstream topology remain future
work.

## Tests and case studies

```sh
julia --project=validation/env --startup-file=no test/runtests.jl
julia --startup-file=no test/bindings.jl
julia --project=validation/env --startup-file=no test/run_conformance.jl
julia --project=validation/env --startup-file=no test/three_dimensional.jl
```

The last three need the native CutCells library: build it with `bash scripts/build_cxx_bindings.sh`
and set `EVOLVINGDOMAINS_CUTCELLS_LIBRARY` to the installed dynamic library, or install a
published artifact. `test/runtests.jl` is pure Julia.

Scientific comparisons and long transients live in `validation/`; they are case studies,
not `Pkg.test` contracts, and each case directory (for example `validation/zalesak`) carries
its own `README.md` with its verdict.

The old `EvolvingDiscreteGeometry`, GridapEmbedded cut adapter, geometry cache, grid/mesh
transfer API, CSG helper API, and plotting extension have been removed. Kernel storage is a
flat value array plus a private grid-metrics value; `getInterpolator` remains an internal
helper, not public API. There is no package plotting extension.

## Known limitations

- **Conservative transport near an edge.** Totals are exact; the last few nodes at the
  computational box and at the support edge are not. At the box the deviation is 0.713 to
  1.403 at the edge and vanishes by the sixth node inward, falling 2 to 3 orders of magnitude
  per node. At a support edge (a density that ends at the interface) it is not erosion but
  overshoot: a constant 7 on a cut support spans (4.8, 10.1) after a single step, with the
  total exact. The source
  paper's remedy is a ghost-cell layer plus an inflow/outflow budget; it is not implemented,
  and the rescaled truncated stencil is an option the authors report but never tested.
  A density whose support ends at the interface therefore keeps its mass and not its
  interface values; a smooth field inside a support that stays away from the box is the case
  this operator is built for.
- **No conservative restriction.** `transfer!` fine→coarse is point sampling; `restrict!` is
  reserved for an averaging operator.
- **Transfer cannot create information.** A coarse level's error is reproduced faithfully on
  the fine level, and `coarse→fine→coarse` is exact while the reverse is not.
- **Levels are a mutable `Dict`.** The shared-domain invariant is enforced when a
  `Discretization` is constructed, not when `d.levels` is written by hand; a hand-grafted
  level is rejected at use.
- **`Grid.axes` is mutable.** Editing an axis invalidates cuts and maps derived from it;
  transport re-checks the geometry, the other operators trust `mesh_id`.
- **Native artifacts** are a release step: without them `cut` needs
  `EVOLVINGDOMAINS_CUTCELLS_LIBRARY` pointing at a library built by
  `scripts/build_cxx_bindings.sh`.
