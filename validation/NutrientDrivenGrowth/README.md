# NutrientDrivenGrowth — nutrient-driven free-boundary growth with a topology change

**Verdict: accepted** for the coupled two-dimensional free-boundary benchmark on the
`64²` baseline (and `32²`, `96²`, `128²` refinements).

The two initial discs **merged numerically**, with no artificial geometry: the level-set
transport produced the topological change on its own, at `t = 0.34` on `64²`
(`0.30 / 0.34 / 0.30 / 0.27` on `32² / 64² / 96² / 128²`), and the perimeter stops growing and
starts decreasing at exactly that moment.  The area-growth identity `d|Ω|/dt = ∫_Ω g(c)`
closes with a **normalized cumulative residual of 5.2 % at `128²`** (`10.4 %` at `64²`,
`17.5 %` at `32²`), and the two independent rate estimates — reconstructed-area difference and
interface flux `∫_Γ v·n` — agree with each other to **≤ 1 %**, so the residual sits in the
velocity recovery / discrete divergence theorem, not in the transport.

```sh
# native CutCells backend (needed by EvolvingDomains.cut); build once
PREFIX=$PWD/validation/lib/cutcells bash scripts/build_cxx_bindings.sh

export EVOLVINGDOMAINS_CUTCELLS_LIBRARY=$PWD/validation/lib/cutcells/lib/libevolvingdomains_cut_cells.dylib
julia --project=validation/env --startup-file=no \
      validation/NutrientDrivenGrowth/run.NutrientDrivenGrowth.jl 64 0.5 n64
julia --project=validation/env --startup-file=no \
      validation/NutrientDrivenGrowth/visualize.NutrientDrivenGrowth.jl n64
```

`run.NutrientDrivenGrowth.jl [ncells] [tfinal] [output-name] [dtmax]`.

## Package provenance is enforced, not assumed

Two EvolvingDomains are reachable from this directory, and only one of them is this project:

| environment | resolves `EvolvingDomains` to | version |
|---|---|---|
| default (`~/.julia/environments/v1.11`, stacked under the case env) | `…/Scientific stack/EvolvingDomains.jl/src/EvolvingDomains.jl` | 0.5.0 (pre-redesign) |
| `--project=validation/env` | `…/EvolvingDomains-Redesign/src/EvolvingDomains.jl` | **1.0.0** |

```sh
julia            -e 'println(Base.locate_package(Base.PkgId(Base.UUID("f9bd8acb-…"), "EvolvingDomains")))'
julia --project=validation/env -e 'println(Base.locate_package(Base.PkgId(Base.UUID("f9bd8acb-…"), "EvolvingDomains")))'
```

`checkEvolvingDomains()` runs before any physics and **refuses to start** unless
(a) `pkgdir(EvolvingDomains)` is this repository, (b) `pkgversion ≥ 1.0.0`, (c) `pathof` is under
that root, (d) the v1.0 surface (`Cut`, `CutSurface`, `cut`, `aggregateCutCells`,
`eulerianTransport!`, `lagrangianTransport!`, `CIPCache`, `CUT_IN/OUT/CUT`) is defined, and
(e) pre-1.0 names (`EvolvingDiscreteGeometry`, `CartesianGridInfo`,
`AbstractVelocitySource`, `aggregate_cut_cells`) are **absent**.  It records the root, version,
source path, `Base.load_path()` and SHA-256 of `src/EvolvingDomains.jl`, `src/DataModel.jl`,
`src/CutData.jl`, `src/Kinematic/Kinematic.jl`, `src/Geometric/Geometric.jl` at the top of every
`run_metadata.txt`.  `src/EvolvingDomains.jl` SHA-256 for all five runs below:
`77f7341b69ba19971f5d82621734661b69dfc55c52ddda56dbb96fc0425a50d4`.

The suite was rerun under the guard and every reported number came out **bit-identical** to the
earlier run (e.g. `64²`: `balance_mean = 0.10568766448280924`, `final_area = 0.2394063537658799`,
`merged_at = 0.34`), confirming that the case was already resolving the repository package and
that the run is deterministic.

## Model

| component | equation | boundary condition |
|---|---|---|
| nutrient | `-D_c Δc + λ c = 0`, `D_c = 0.02`, `λ = 0.05` | `c = 1` on `Γ` |
| growth | `g(c) = g₀ max(c − c*, 0)`, `g₀ = 1.5`, `c* = 0.2` | — |
| pressure | `−Δp = g(c)` (Darcy: `v = −∇p`) | `p = 0` on `Γ` |
| interface | `∂_t φ + v⋅∇φ = 0` (`V_n = v⋅n`) | `φ = min(‖x−c₁‖, ‖x−c₂‖) − R₀` at `t = 0` |
| marker | `∂_t q + v⋅∇q = 0` | `q₀` = 0.8/0.2 Gaussians at the two centres |

`c₁ = (0.33, 0.50)`, `c₂ = (0.67, 0.50)`, `R₀ = 0.15`, `D = [0,1]²`.

Parameters as required.  With them the nutrient is nearly uniform inside the tissue
(`k = √(λ/D_c) = 1.58`, diffusion length `0.63 ≫ R₀`), so `c ∈ [0.97, 1.01]` and `g` is almost
uniform: this is a *pressure-driven* growth, and the nutrient solve is genuine but smooth.

## Architecture, and how the representations are connected

Per step, `dynamics(levelset, backdrop, cfg)`:

1. `ED.cut(field, :init)` → owned `Cut` snapshot (public arrays only).  From its **triangular
   sub-cells** the physical area, inside-area centroid and inside fraction of every background
   cell are computed (`cutGeometryInfo`), the connected components are counted
   (`componentCounts`), and `aggregateCutCells` is probed on the live topology.
2. The same nodal array — `vec(φ)`, first axis fastest — becomes the free-dof vector of a
   Gridap `Q1` space on a `CartesianDiscreteModel` with the identical `ncells²` partition;
   `DiscreteGeometry(φh, model)` → `GE.cut` → `Triangulation(PHYSICAL_IN)`,
   `Triangulation(ACTIVE_IN)`, `EmbeddedBoundary`, `GhostSkeleton` and their measures.
   **Documented conversion**: the Q1 dof order of the model equals the node order of the
   EvolvingDomains array; this is asserted at start-up (`checkNodeOrdering`) and in
   `backgroundNodal`, and would raise instead of silently degrading.
3. `solveNutrient` → CutFEM (symmetric Nitsche + Burman–Hansbo facet ghost penalty, both
   scaled by `D_c`) on the active P1 space.
4. `computeGrowthRate` → `g(c)` as an FE field on the same space.
5. `solvePressure` → same CutFEM machinery with `−Δp = g` and `p = 0`.
6. `velocityFromPressure` → two objects: the smooth per-cell Darcy field `−∇p` sampled with the
   bilinear derivative at each cell's inside-area centroid (diagnostics, plots, divergence test)
   and the **nodal transport velocity**, a zero-order normal extension in which every node takes
   the velocity of the *nearest interface facet* (the owning cell's sample).  The nodal field is
   what `eulerianTransport!` and CIP consume.
7. `advanceGeometry!` → `ED.eulerianTransport!` (WENO5 + SSP-RK3) then `ED.reinitialize!`.
8. `advanceMarker!` → `ED.lagrangianTransport!` (intensive CIP) with **two** `Field` buffers on
   the same `Discretization`, swapped every step.

Nothing is reconstructed from the initial circles after `t = 0`; every step re-cuts the current
level set.

## Results

Baseline `64²`, `h = 1/64`, `dt = 5e-3` (fixed by `dtmax`; the velocity CFL allowed
`0.4·h/v_max ≈ 0.06`), 100 steps, `t ∈ [0, 0.5]`, wall time **54 s**, of which ~40 s is
first-call compilation.  Steady state: 0.46 s per step at `64²`, 0.54 s at `128²`, split roughly
as geometry re-cut + space setup 0.29 s, the two operator assemblies + level-set/marker transport
+ diagnostics 0.16 s, the two sparse direct solves 4 ms each (`52 / 61 ms` at `128²`).
`EvolvingDomains.cut` itself takes 1.3 ms per step at `64²`.

| quantity | 32² | 64² | 64², dt=2e-3 | 96² | 128² |
|---|---|---|---|---|---|
| `|Ω(0)|` (exact `2πR₀²` = 0.14137) | 0.14007 | 0.14106 | 0.14106 | 0.14123 | 0.14129 |
| `|Ω(0.5)|` | 0.22818 | 0.23941 | 0.23943 | 0.24454 | 0.24715 |
| `|Γ|` peak (time) | 2.179 (0.305) | 2.255 (0.345) | 2.255 (0.354) | 2.220 (0.300) | 2.190 (0.270) |
| first single-component time | 0.305 | 0.340 | 0.354 | 0.300 | 0.270 |
| `|Γ(0.5)|` after the neck rounds | 2.102 | 2.074 | 2.073 | 2.035 | 2.019 |
| mean step residual `|Δ|Ω|/Δt − ∫g| / |∫g|` | 0.177 | 0.106 | 0.107 | 0.071 | 0.053 |
| cumulative `Δ|Ω| / ∫∫_Ω g` over the recorded steps | 0.825 | 0.896 | 0.893 | 0.930 | 0.948 |
| mean `∫_Γ v·n / ∫_Ω g` | 0.833 | 0.900 | 0.898 | 0.931 | 0.948 |
| `c` range | 0.976–1.008 | 0.973–1.006 | 0.973–1.006 | 0.971–1.005 | 0.969–1.006 |
| `p` range | −3.9e−3–1.15e−2 | −2.8e−3–1.30e−2 | −2.8e−3–1.31e−2 | −2.5e−3–1.41e−2 | −2.9e−3–1.47e−2 |
| marker `q` range | −9.7e−7–0.796 | −7.6e−13–0.801 | −8.3e−13–0.801 | −5.0e−15–0.802 | −1.6e−27–0.801 |
| max algebraic residual (nutrient / pressure) | 3.7e−16 / 3.6e−15 | 3.8e−16 / 1.5e−14 | — | 4.2e−16 / 3.5e−14 | 4.5e−16 / 7.3e−14 |

Interpretation:

* **Transport is exact to ~1 %.**  `Δ|Ω|/Δt` and `∫_Γ v·n` agree at every resolution
  (0.825 vs 0.833, 0.896 vs 0.900, 0.930 vs 0.931, 0.948 vs 0.948).  WENO5+SSP-RK3 moves the
  reconstructed interface at the velocity it is given; `reinitialize!` changed the physical area
  by **0.0** in every step of every run.
* **The remaining 5–17 % is the recovered velocity.**  Both rate estimates are built from the
  same discrete `∇p`, so their common deficit against `∫_Ω g` measures how well the recovered
  field satisfies the discrete divergence theorem: a first-order-biased on the interface
  (the sample sits at the inside centroid, up to `h/2` inside; `v ∝ g r/2`), decreasing
  monotonically with `h` (17.5 → 10.4 → 7.1 → 5.2 %).  Halving `dt` at `64²` moves the residual
  by 0.1 % (0.106 → 0.107) and the final area by 8e−5 relative: the residual is spatial.
* **Elliptic solvers are accurate.**  On the initial geometry, against the closed-form solutions
  of the same operators on the two discs (`p = (R₀²−r²)/4`, `c = I₀(kr)/I₀(kR₀)`):
  `|Ω|` error 0.22 % at `64²` (0.06 % at `128²`); pressure nodal max error 2.0e−5
  (0.36 % of `p_max`), relative L2 error 1.3e−3; nutrient relative L2 error 1.1e−5.
* **Divergence mismatch** `‖∇·v − g‖`: central differences of the cell-sampled velocity on
  fully-inside cells (3×3 neighbourhood `CUT_IN`), against bilinear-interpolated `g`:
  mean 1.8–2.5 %, max 5.8–6.8 %.  It does not decrease with `h` — it is a property of the wide
  central stencil on Q1 gradients, not a residual of the solve.
* **Topology.**  `components_nodes` (4-connected nodes with `φ < 0`) flips `2 → 1` at the merge;
  `components_cells` (4-connected active cells, i.e. in + cut) flips about 2.5 grid cells
  earlier, because a cut band already connects the two regions while no node is negative yet.
  Both are resolution dependent; both are reported per step.
* **Small cut cells.**  `h²`-normalized minimum inside-area fraction reached 1e−9 in every run
  (an interface passing almost exactly through a node).  The Nitsche+ghost-penalty formulation
  stayed stable: `c` and `p` remained bounded and the algebraic residuals stayed at
  10⁻¹⁴–10⁻¹⁶.  `aggregateCutCells` ran on every snapshot (500–3900 distinct roots) but was not
  needed for stabilisation, so no aggregated space was built.
* **Post-merge velocity.**  After merging, the largest interface velocity samples (0.24 / 0.31
  at `64² / 128²) sit in the concave re-entrant corner at the neck, a crease inherited from
  `min(d₁,d₂)`; a corner singularity, resolved more sharply on finer grids, whose global effect
  on the area balance is below the reported residual.

## Figures

`output/n64/` and `output/n128/` (also `n32`, `n96`):

* `fig_domain_evolution.png` — 8-panel montage: growth, approach, merge, rounding of the neck.
* `fig_nutrient.png` — `c` on `Ω(t)`, masked outside the phase.
* `fig_pressure_velocity.png` — `|v|` (per-cell recovery), pressure contours, and velocity
  arrows of the nodal transport field in the interface band only.
* `fig_marker.png` — passive marker `q` (CIP).
* `fig_metrics.png` — area, perimeter (peak at the merge), the three area-growth rate
  estimates, and the component counts vs time.

Material boundary: solid black.  Computational box `∂D`: thin dashed grey.  Field values outside
`Ω` are `NaN` (transparent).

## Packages and how they were reached

* `EvolvingDomains 1.0.0` (**current repository**, `[sources] path = "../.."` in
  `validation/env`, enforced and hash-recorded by `checkEvolvingDomains`): `Grid`,
  `Discretization`, `Field`, `cut`, `Cut` (sub-cells, boundary facets + outward normals),
  `aggregateCutCells`, `eulerianTransport!`, `WENO5Cache`, `lagrangianTransport!`, `CIPCache`,
  `reinitialize!`.
* `Gridap 0.19.7` + `GridapEmbedded 0.9.8` (Algoim 0.2.2): read from the *stacked* default
  environment `~/.julia/environments/v1.11`, not from `validation/env` (which has no Gridap dep).
  `CartesianDiscreteModel`, `FESpace`/`FEFunction`, `DiscreteGeometry`, `cut`,
  `Triangulation(…, PHYSICAL_IN / ACTIVE_IN)`, `EmbeddedBoundary`, `GhostSkeleton`,
  `Measure`, `get_normal_vector`, `AffineFEOperator`, `solve`, `get_matrix`, `get_vector`.
* `CairoMakie 0.15.17`, `StaticArrays 1.9.24`, Julia 1.11.7.
* Native CutCells backend built locally into `validation/lib/cutcells/lib/` (no published
  artifact): `EVOLVINGDOMAINS_CUTCELLS_LIBRARY` must point at it.

No package source, project file or manifest was modified.  `StaticAssembler` and `Curvature`
were not used: the former was unnecessary once Gridap assembly worked, the latter is not present
in either environment, so the optional `p = γκ` experiment was not attempted.

## Public API limitations observed

* `Grid`/`Field` arrays carry no metadata; the Q1 dof ↔ node correspondence between
  EvolvingDomains and Gridap is an undocumented coincidence of ordering, verified at runtime
  here rather than guaranteed by either package.
* `Cut` is surface geometry only: it has no interface → cell adjacency beyond
  `facet_to_bgcell`, and its sub-cells are not a triangulation with dofs, so it cannot be fed
  to an FE assembly directly.  A doc conversion (level set array → `DiscreteGeometry`) is the
  only route, and it re-cuts at every step.
* `GridapEmbedded`'s cutter has no locality: the whole 4096-cell background model is re-cut every
  step (0.29 s per step at `64²` including space setup) although only 16–24 % of the cells are
  active, and no part of the `Cut` from the previous step can be reused.
* `aggregateCutCells` returns root ids for cells but no aggregated space, so a user who needs
  AgFEM-style stabilisation must build the aggregated basis themselves.
* The transport operators accept a single-valued nodal velocity only, which cannot represent two
  opposing interfaces closer than one cell; this is what limits the pre-merge gap resolution and
  is visible as the `components_cells` / `components_nodes` discrepancy.
