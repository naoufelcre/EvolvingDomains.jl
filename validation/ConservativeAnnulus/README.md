# ConservativeAnnulus — conservative cut-to-cut transport

**Verdict.** The conservative `TransportMap` path runs end to end and is exactly conservative
in the sense the package claims — the unweighted nodal sum is invariant to 1 ulp of itself.
That is *not* conservation of the physical mass: the same run loses 2.0% (64² cells) and
0.85% (128² cells) of `∫ ρ dΩ`, and the two facts are independent measurements here, not
assumptions. The dominant error source is the two-node edge layer of the transported band,
not the cut geometry: the same quadrature applied to the *exact* density on the *same*
reconstructed geometry drifts by 3.8e-4 / 6.6e-6.

## Files

| File | Contents |
|---|---|
| `run.ConservativeAnnulus.jl` | analytic model, multilevel study, `simulate` (the 7 §3 steps inline), diagnostics, entry point |
| `visualize.ConservativeAnnulus.jl` | CairoMakie figures, reads `output/results.jls` |
| `lib/libevolvingdomains_cut_cells.dylib` | native cut backend, built locally (the checkout publishes no artifact) |
| `output/` | `results.jls` and five PNGs |

## Running

```sh
julia validation/ConservativeAnnulus/run.ConservativeAnnulus.jl
julia validation/ConservativeAnnulus/visualize.ConservativeAnnulus.jl
```

Both activate `validation/env`. `run` sets `EVOLVINGDOMAINS_CUTCELLS_LIBRARY` to `lib/` when
the caller has not, so a machine without a published CutCells artifact still works. To
rebuild that library: `BUILD_DIR=/tmp/b PREFIX=validation/ConservativeAnnulus bash
scripts/build_cxx_bindings.sh` (cmake fetches CutCells), then keep `lib/`.

Wall clock: 0.33 s (quick) + 1.7 s (fine) + figure/transfer work, ~14 s total on an M-series
laptop; 545 MiB and 4.3 GiB allocated by the two loops.

## The distinction the case exists for

The operator moves one nodal vector (`vec(field.data[level])`) to another. Its invariant is
`Σ q` over the map's source support, and `q` is a nodal density with **no weights** — the
support is the set of nodes that are corners of an in-phase or straddling background cell.
Two independent measurements follow:

| quantity | 64² cells | 128² cells |
|---|---|---|
| Σ q drift per step | ≤ 3.3e-16 | ≤ 4.4e-16 |
| `h²·Σq` vs `∫ρ` at t=0 | +32.8% | +16.9% |
| physical `M_h` drift at t=1 | −2.00e-2 | −8.51e-3 |
| same, exact ρ on the same cut | −3.78e-4 | −6.59e-6 |
| same, void-side nodes taken from the exact solution | −2.16e-2 | −9.29e-3 |

`h²·Σq` is a uniform-weight Riemann sum over the *band*, so it counts the halo nodes that the
cut region does not contain; it is O(h) wrong and would be a 33% mass error at 64². The
physical integral instead uses the per-subcell triangles of the cut: `Material = Σ_K |K| ·
ρ_h(x_K)` over `CUT_IN` subcells, centroids of the three physical vertices, with `ρ_h` the
`linearInterpolator` of the nodal array. Centroid quadrature is exact for a linear density
and second order for this one. Replacing the void-side nodes by their exact values does not
help, so the loss is in the transported material values themselves, not in the shared
cut-cell corners.

## Recorded diagnostics

Geometry (`Cut` subcell triangles vs the analytic annulus): area error 3.80e-4 → 7.89e-6
relative, centre error 2.02e-5 → 1.89e-6, both over 55 / 109 steps.

Density (nodal values vs `ρ₀(x − t v)`), the cut band has a boundary layer about two nodes
thick and this geometry is only 8 cells wide at 64²:

| φ/h bin at t=1 (128²) | nodes | L∞ | L¹ |
|---|---|---|---|
| [−6, −4) | 649 | 7.90e-2 | 1.84e-2 |
| [−2, −1) | 290 | 2.35e-1 | 6.01e-2 |
| [−1, 0) | 299 | 3.01e-1 | 7.97e-2 |
| [0, +1) (void-side halo, not material) | 281 | 7.91e-1 | 1.50e-1 |

Positivity holds where it is meaningful: minimum over all steps and in-phase nodes 0.715
(128²) / 0.700 (64²); on the whole support, which includes void-side halo nodes, 0.480.
Extrema at t=1 on in-phase nodes [0.742, 1.301] against exact [1.000, 1.191], i.e. +0.110
overshoot. Constant-state preservation measured on the same map: `max|target − 1|` = 0.862 on
the band, 0.039 one node inside — the documented partition-of-unity defect at the band edge.

Multilevel: `transfer!` of `exp(sin2πx + cos2πy)` from 65² to 129² nodes gives L∞ 1.77e-2,
and from 129² to 257² gives 4.44e-3 — observed order 1.99. Fine → coarse is exact node
sampling (L∞ 0.0), confirming it is not a conservative restriction. Prolongating the quick
run's transported density to 129² differs from the fine run by L∞ 1.50e-1 on fine core nodes:
coarse truncation error is carried up, not removed.

Timestep: the nodal sum is invariant at CFL 0.2/0.4/0.8 (≤ 3.3e-16) while the mass drift moves
only from −1.90e-2 to −2.15e-2, so the reported conservation failure is not a dt artifact.

## Limitations of the public cut representation

- **No per-node material weights.** `Cut` exposes cell states, subcells and physical vertices,
  not the fraction of a node's dual cell inside the phase. A nodal quantity therefore cannot be
  turned into a physical integral; only cell-based quadrature is available (used here).
- **No ghost layer.** The map's support is one cell wide while its reconstruction stencil is
  four nodes wide, so the outermost two node layers are renormalized rather than conserved
  locally; that layer is where the density error and the mass drift come from. Values there
  reach 1.79–1.93 where the exact density is ≤ 1.19.
- **The band, not the interface, is the transport support.** Void-side nodes carry transported
  values that have no physical meaning; they are part of `Σ q` and of the pull stencils.
- The map's nodal sum would also silently drop any source value outside its support, so the
  source vector here is zeroed there explicitly (all dropped nodes lie outside the material).
- Copying a `Cut` snapshot per step doubles the cut cost; both are constructed per step here
  for the sake of an explicit loop.
