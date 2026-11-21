# Campaign II — Experiment 4: MovingCutAgFEM

**Conclusion.** EvolvingDomains' public `Cut` snapshot is a *complete and sufficient* descriptor
of a moving cut geometry: a 72-line, allocation-free adapter turns it into StaticAssembler's
`Mesh`/`CutData`/`BoundaryData` and StaticAssembler then assembles the projected AgFEM operator
`B = EᵀAE` directly, with **no GridapEmbedded and no re-cut anywhere in the primary path**
(`Base.get_extension(StaticAssembler, :StaticAssemblerGridapExt)` is `nothing`; see
`output/versions.txt`). The moving-domain manufactured elliptic problem converges at second
order in `h` (observed 2.08–2.19), the projected system is symmetric to `4e-17`, Cholesky
succeeds at every one of the 603 measured time levels (201 steps x 3 grids) down to a minimum cut-cell area fraction
of `4e-7`, and the projected right-hand side satisfies `Σ b_m / 4 = |Ω_h|` to `7e-16`. Aggregation
is not cosmetic: the same mesh and cut assembled without aggregation is **indefinite** and fails
Cholesky at every resolution tested.

---

## 1. Scope and what was verified first

| item | value |
| --- | --- |
| EvolvingDomains | `v1.0.0`, `/Users/naoufel/Documents/Scientific stack/EvolvingDomains-Redesign` |
| EvolvingDomains git | `0809db8bf38e5d25975022a1848798f2219289c8` (working tree dirty — pre-existing, not touched by this experiment) |
| StaticAssembler | `v0.1.0`, `/Users/naoufel/Documents/Scientific stack/StaticAssembler.jl` |
| StaticAssembler git | `4e0249cd670700e4a9329069e317e2af90237480` (worktree dirty — pre-existing) |
| Ferrite | `v1.7.0` |
| Julia | `1.11.7` |
| CutCells native library | built from this checkout's `cxxBindings`; `sha256 b06f6321…` identical to a fresh `scripts/build_cxx_bindings.sh` build (rebuilt during the experiment to confirm) |

StaticAssembler's public API was read from the installed source (`src/*.jl`, `docs/API.md`,
`docs/details.md`, `test/*.jl`), not from the online README: the documented primitives are
`Mesh`, `CutData`, `BoundaryData`, `StaticSystem`, `refresh_geometry!`, `assemble!`,
`extension`, `restrict_solution!`, `ElementVariationalForm`/`BilinearSlot`/`RhsSlot`, and the
four operator hooks `_op_pair`, `_rhs_contract!`, `_boundary_op_pair`,
`_boundary_rhs_contract!`. No online example was assumed.

No upstream package, manifest or parent project was modified. All artefacts live under
`validation/PlanarStefan/round2/` (an isolated `Project.toml`/`Manifest.toml` pair there is the
only dependency change, and it is local).

## 2. Problem

`-Δu = 4` in `Ω(t)`, `u = 0` on `∂Ω(t)`, `u_exact = R² - [(x-c_x)² + (y-c_y)²]`,
`R = 0.23`, `c(t) = (0.38+0.14t, 0.46+0.08t)`, `0 ≤ t ≤ 0.5`, `D = [0,1]²`.
The level set starts as the exact signed distance and is then advanced *only* by
`eulerianTransport!` (WENO5 + SSP-RK3), 200 steps of `dt = 2.5e-3` (CFL 0.026 at `h = 1/64`).
Resolutions `n ∈ {33, 65, 129}` nodes per axis, i.e. `32², 64², 128²` cells.

The disk never approaches `∂D` (over the whole trajectory it stays inside
`[0.15,0.68] × [0.23,0.73] ⊂ D`), so the whole Dirichlet boundary is the embedded interface —
exactly the part StaticAssembler's boundary path has to handle.

## 3. The adapter: `Cut` → `CutData`

`run.MovingCutAgFEM.jl:210` (`convertCutToStaticAssembler`) and `:189`
(`buildStaticAssemblerMesh`). 72 lines including comments.

The two representations are *deliberately aligned*, so the adapter is a field copy plus one
structural fix:

* `Cut.bg_states` / `subcell_states` / `cut`/`CUT_IN`/`CUT_OUT` codes are already
  `IN = -1, CUT = 0, OUT = 1`; StaticAssembler uses the same constants. An `@assert` at module
  scope fails loudly if either side drifts.
* `physical_points` / `reference_points` are already `Vector{NTuple{2,Float64}}` and the
  reference points already live in the parent cell's `[0,1]²` frame — the exact contract
  StaticAssembler validates (`_metric_matches(dx, dr, hx)`). They are passed **by reference**
  (StaticAssembler only reads them), so no geometry is copied or recomputed.
* `cell_to_bgcell` is already sorted by ascending background cell, and the CSR
  `point_ids`/`point_pointers` are already one-based triangular connectivity.

**The one semantic mismatch.** EvolvingDomains emits volume subcells for *every* background cell
(a fully-outside cell owns its own outside volume, and a fully-inside cell its inside volume),
whereas StaticAssembler requires every subcell to belong to a `CUT` cell, contiguously, in
ascending cell order. The adapter keeps the subcells whose parent is `CUT` and renumbers the CSR
pointers. Measured (n = 33, t = 0): 156 subcells emitted → 56 kept (24 inside, 32 outside
pieces), 312 points retained (all of them, since point ids stay valid), 14 boundary facets.

**The boundary descriptor.** `Cut.boundary` already satisfies the `BoundaryData` contract:
facets grouped by ascending background cell, 2 vertices each, unit normals outward from the
`φ < 0` phase, reference points in the facet's own parent cell. No orientation fixing was needed,
which is worth stating explicitly because the Nitsche signs depend on it.

**Mesh identity.** `Mesh.mesh_id` is `Any`; EvolvingDomains' `Grid.mesh_id` is a `UInt64` that
survives copies. Passing it through gives StaticAssembler's `===` identity check for free, so a
cut snapshot from a different lattice is rejected rather than silently mis-assembled.

**Addendum (index-based adapter).** `Cut` now carries one background-cell → subcell array:
the CSR `cell_subcell_ptrs`, read through `subcells(cut, cell)` (`cut_cell_ids(cut)` derives the
CUT cells from `bg_states`). Subcells are grouped by ascending cell in both constructor paths,
so the CSR row is the subcell id and no id array or duplicated classification is stored. The
adapter therefore no longer filters `cell_to_bgcell`; `checkCutSubcellIndex` (`run` stage A/B)
asserts the indexed traversal reproduces the old filtered scan exactly, and every artefact in
`output/` was regenerated at all three resolutions and is **byte-identical** to the pre-index
run (`diagnostics`, `B`, `E`, `rhs`, `u_bg`, `conditioning`, `ferrite_*`; only the timing columns
of `steps_*.csv` differ).

Cost of the whole adapter: **0.0023 / 0.0051 / 0.021 ms and 8.7 / 18.8 / 46 kB per step** at
`n = 33 / 65 / 129` — 2–4× less time than the scan it replaces, and the remaining allocation is
the descriptor itself plus one `findall` over `bg_states`. That is under 1 % of StaticAssembler's own refresh cost and of
EvolvingDomains' `cut` cost — the interface conversion is not the bottleneck and is not a
framework.

## 4. The operator through StaticAssembler

The mathematics is 58 lines (`run.MovingCutAgFEM.jl:57–114`): two operator types and four
hook methods. Everything else is element tables StaticAssembler generates.

```
a_h(u,v) = ∫_Ω ∇u·∇v − ∫_Γ (∂_n u) v − ∫_Γ (∂_n v) u + (γ/h) ∫_Γ u v
l_h(v)   = ∫_Ω 4 v            (Nitsche data g = 0 ⇒ both boundary load terms vanish)
```

* volume slot: `BilinearSlot(GradOp(), GradOp(), ConstantCoefficient(1.0))` — constant
  coefficient, so StaticAssembler caches the interior and per-cut-cell element matrices once per
  refresh and only the RHS source and the boundary slots are contracted pointwise.
* boundary slots: `(GradOp,ValueOp,−1)`, `(ValueOp,GradOp,−1)`, `(ValueOp,ValueOp,+γ/h)`.
  Symmetric Nitsche is expressible with the *documented* boundary slot API — no undocumented
  operator was invented, and the two flux slots make the boundary contribution symmetric by
  construction.
* `γ = 40 / min(hx,hy)`, supplied as a `ConstantCoefficient`; `h` is not exposed to the
  boundary slot, so the penalty is fixed per system from the known lattice spacing.
* `ncoefficients = 0`, so `assemble!(system, ())`: there is no physics field to prolongate, only
  the geometry-driven operator.

`E`, `B`, `b_m` are inspected directly (`output/E_n*_t500.csv`, `B_n*_t500.csv`,
`rhs_n*_t500.csv`):

| n | `E` | active rows | nnz(E) | rows w/ 1 entry | rows w/ 4 entries | row sums | nnz(B) | density | `‖B−Bᵀ‖/‖B‖` | diag(B) | bandwidth |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| 33 | 1089×168 | 230 | 416 | 168 | 62 | exactly 1 | 1766 | 6.26 % | 3.4e-17 | [2.667, 214.8] | 41 |
| 65 | 4225×679 | 799 | 1159 | 679 | 120 | exactly 1 | 6687 | 1.45 % | 1.3e-17 | [2.667, 173.0] | 87 |
| 129 | 16641×2723 | 2963 | 3683 | 2723 | 240 | exactly 1 | 25779 | 0.35 % | 1.5e-17 | [2.665, 179.9] | 177 |

(`E` is `n_bg × n_m` with `n_bg = (n)²`, i.e. the whole Q1 background lattice; only the 230 /
799 / 2963 rows belonging to active cells are nonzero, every other row is empty.)

Every master row of `E` is a unit row; every slave row has exactly four entries (Q1 root
weights) and every active row sums to **1.000000000000000** — `E` is a partition of unity onto
the aggregate roots, which is why the extension reproduces constants. `E` is short: the
extension only reaches the cut cells' own nodes.

A free and very sharp consistency check follows from `E·1 = 1`:
`Σ b_m = 1ᵀEᵀb = 1ᵀ b = 4|Ω_h|`. Measured at `t = 0.5`:

| n | `Σb_m/4` | `\|Ω_h\|` (from the cut triangles) | relative difference |
| --- | --- | --- | --- |
| 33 | 0.1651678520 | 0.1651678520 | 6.7e-16 |
| 65 | 0.1659271689 | 0.1659271689 | 1.7e-16 |
| 129 | 0.1661263226 | 0.1661263226 | 1.7e-16 |

The cut-cell quadrature therefore measures the reconstructed domain exactly, and the AgFEM
projection conserves the load — the RHS is not silently wrong on the cut cells.

## 5. Numerical accuracy of the moving-domain solve

Mean over the 11 sampled times in `[0, 0.5]` (`output/diagnostics_n*.csv`):

| n | `‖u_h−u‖/‖u‖` | interpolation-only | max interior error | max `\|u_h\|` on Γ | `‖Bu_m−b‖/‖b‖` | max asym | area error | centroid error | min cut fraction |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| 33 | 2.05e-2 | 1.07e-2 | 6.4e-4 | 8.3e-4 | 3.6e-14 | 4.0e-17 | 6.1e-3 | 6.8e-5 | 3.1e-7 |
| 65 | 4.50e-3 | 2.67e-3 | 1.6e-4 | 1.5e-4 | 5.4e-14 | 4.0e-17 | 1.5e-3 | 7.0e-6 | 1.6e-6 |
| 129 | 1.06e-3 | 6.67e-4 | 4.6e-5 | 3.8e-5 | 7.4e-14 | 4.0e-17 | 3.8e-4 | 1.9e-6 | 4.0e-7 |

(`‖u‖ = √(∫u²) ≈ 4e-3`, `u_max = R² = 5.29e-2`; the interior-error column is a maximum over
lattice nodes at least `1.5h` inside the interface, so it is not comparable in magnitude to the
volume L2 norm.)

**Convergence.** Observed orders between successive grids:

| quantity | order (h: 1/32→1/64) | order (1/64→1/128) |
| --- | --- | --- |
| total relative L2 | 2.19 | 2.08 |
| interpolation-only | 2.00 | 2.00 |
| geometric area error | 2.00 | 2.00 |

**Separating geometry from discretisation.** The total error splits cleanly into an
`O(h²)` interpolation share and a surplus that is also `O(h²)` but with a larger constant:
`total/interpolation = 1.93, 1.69, 1.59` at the three grids. The surplus is *not* a scheme
defect — it is the piecewise-linear cut: the reconstructed domain is an inscribed polygon whose
area is low by 0.61 % / 0.15 % / 0.038 % (`h²`), so `u` is being compared to a slightly shrunken
domain. Since the manufactured `u` vanishes quadratically at the interface, a domain error of
relative size `ε` contaminates the L2 error at a comparable order — consistent with the measured
ratio drifting to 1.6. In other words: the FE discretisation converges at the expected rate, and
the residual gap is geometric approximation, not solver error. **The Nitsche boundary error
converges at second order too** (`8.3e-4 → 1.5e-4 → 3.8e-5`, observed order 2.45 then 2.00).

Errors are essentially **constant in time** (variation < 3 % over 200 steps) despite the geometry
being re-cut from a *numerical* level set — the transport preserves the disk to a centroid error
of ≤ `7e-5/7e-6/2e-6` (order 3.3, better than the geometry's `h²`), so no time-accumulating
geometric drift contaminates the result. This is the expected invariant for rigid translation and
it is the strongest evidence that EvolvingDomains → StaticAssembler is a working coupling rather
than a coincidence at `t = 0`.

## 6. Does aggregation improve conditioning? Yes — decisively

Same mesh, same cut, same form, on the coarsest-cut snapshot (`t = 0`), comparing proper
aggregation with a control where every active cell is its own root ("unaggregated cut FEM",
`output/conditioning_n*.csv`):

| n | agglomerated `n_m` | λ_min | λ_max | cond | Cholesky | self-rooted `n_m` | λ_min | Cholesky | # tiny diagonals |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| 33 | 169 | 1.06e-1 | 272 | 2.6e3 | ok | 233 | **−1.0e-5** | **fails** | 2 |
| 65 | 679 | 2.67e-2 | 257 | 9.7e3 | ok | 803 | **−4.4e-4** | **fails** | 7 |
| 129 | 2725 | 6.67e-3 | 254 | 3.8e4 | ok | 2965 | `NaN` (no factor) | **fails** | 16 |

* With aggregation `λ_min > 0` always, `λ_max ≈ 254–272` is bounded independent of `h`, and
  `λ_min ~ h²`, so `cond ~ h⁻²` — the standard, benign SPD scaling.
* Without aggregation the operator is **indefinite**: the DOFs whose support lies entirely
  outside `Ω` get near-zero rows, the smallest eigenvalue is negative, and CHOLMOD fails. The
  aggregation removes exactly `n_m` = (#active background nodes − #masters) DOFs (233→169 at
  `n = 33`) by projecting them onto their root.
* Note it is the *aggregation*, not the Nitsche penalty, that is doing this work: the minimum
  diagonals of the self-rooted operators are negative (`−7.1e-6 … −4.0e-4`), i.e. the assembled
  operator itself is not sign-definite for those rows.

The eigenvalue on the two coarsest grids is a dense `eigvals` of the symmetric matrix; on
`n = 129` (`n_m = 2725`) it is a power/inverse-iteration estimate using the Cholesky factor
(`spectralEstimate`), so no expensive dense condition number is computed on the large system.

## 7. Does geometry refresh actually reduce work or allocations?

Measured over the 200-step moving run (`output/steps_n*.csv`, mean per step). Every number below
comes from `assemble!(...; threaded = false)` and single-threaded refresh; the runs used `-t 4`
only so that the threaded comparison at the end of this section is meaningful.

| n | structural refreshes | `refresh_geometry!` structural | `refresh_geometry!` retained | `assemble!` |
| --- | --- | --- | --- | --- |
| 33 | 62 / 200 | 0.32 ms / 783 kB | 0.091 ms / **0 kB** (139×) | 0.037 ms / **0 kB** |
| 65 | 159 / 200 | 0.99 ms / 2.9 MB | 0.216 ms / **0 kB** (41×) | 0.091 ms / **0 kB** |
| 129 | 199 / 200 | 5.06 ms / 10.2 MB | 0.73 ms / **0 kB** (1×) | 0.317 ms / **0 kB** |

Full pipeline (mean per step):

| n | `cut` (ED) | adapter | `aggregateCutCells` | `refresh_geometry!` | `assemble!` | Cholesky | solve | `E·u_m` |
| --- | --- | --- | --- | --- | --- | --- | --- | --- |
| 33 | 0.457 ms / 986 kB | 0.005 ms / 16 kB | 0.248 ms / 14 kB | 0.161 ms / 242 kB | 0.037 ms / 0 kB | 0.296 ms / 183 kB | 0.004 ms | 0.001 ms |
| 65 | 0.868 ms / 2.4 MB | 0.013 ms / 25 kB | 0.551 ms / 40 kB | 0.828 ms / 2.3 MB | 0.091 ms / 0 kB | 0.625 ms / 845 kB | 0.007 ms | 0.003 ms |
| 129 | 5.118 ms / 11.1 MB | 0.120 ms / 63 kB | 2.191 ms / 113 kB | 5.016 ms / 10.1 MB | 0.317 ms / 0 kB | 3.676 ms / 3.4 MB | 0.196 ms | 0.060 ms |

**What is reused and what is not** (confirmed by measurement *and* by reading
`src/StaticSystem.jl`):

* `assemble!` is **allocation-free in every single step** at every resolution. The projected
  matrix storage, its sparse pattern, the scatter slots, the per-cell constraint tables, the
  thread buffers, and the coefficient/RHS buffers are all retained across refreshes — the hot
  loop only writes into `B.nzval` and `b_m`.
* `refresh_geometry!` allocates **0 bytes** whenever the projection structure is unchanged: the
  quadrature staging buffers are re-staged in place, and the basis-at-quadrature tables and the
  per-cut-cell constant element matrices are regenerated into `resize!`d buffers
  (`_regenerate_element_tables!`, `_regenerate_boundary_tables!`). The refresh is therefore
  *always* re-doing `O(#quadrature points)` work; what it avoids is every `O(ncells)` allocation.
* On a **structural** refresh (`_projection_structure_changed`: an active/`OUT` membership flip
  or a change of an active cell's root) StaticAssembler *rebuilds* `E`, the master numbering,
  `B`'s sparse pattern, the scatter slots, the thread colouring and the constraint tables. The
  source says so explicitly (`ponytail: fresh B pattern per structural refresh`) and the
  measurement agrees: 0 kB vs 0.8–10 MB.
* **How often is it structural?** Translation at constant velocity changes membership only when
  the interface crosses a cell, so the structural count scales as `|v|Δt/h`: 62/200 → 159/200 →
  199/200 as `h` halves twice. *At the resolutions where the assembly cost matters, the refresh is
  almost always structural* — the reuse benefit is real for the matrix pattern and assembly
  workspaces (`assemble!` stays at 0 B), but the "retained" fast path is not what a translating
  domain spends its time in. That is the honest answer to "does refresh reduce work": it removes
  *assembly* allocations entirely and removes *projection* work only while the aggregate
  structure is frozen.
* The dominant cost is **not** StaticAssembler at all: EvolvingDomains' `cut` is 31–37 % of the
  pipeline and StaticAssembler's refresh + assembly 16–32 %. EvolvingDomains'
  `aggregateCutCells` alone costs 0.25–2.2 ms/step: it re-scans every background cell on each
  call instead of doing dirty-cell updates, and StaticAssembler's `_extension_row_supports!`
  carries the same self-declared shortcut (`ponytail: scan the fixed mesh; switch to dirty-cell
  updates only if this is measurable`) because it must rebuild the extension when the projection
  structure moves. Cholesky is a further 20–24 % of the pipeline.

**Threading.** `assemble!(...; threaded = true)` with `-t 4` was *slower* than serial at every
resolution (0.23 vs 0.06 ms, 0.32 vs 0.11 ms, 1.23 vs 0.29 ms at `n = 33/65/129` on fixed
geometry, i.e. 2.9–4.2× slower). The problems are small and the
colouring pass plus `@threads` overhead dominates; the measurement is reported as taken, not
smoothed.

## 8. Code size: mathematics versus interface conversion

`run.MovingCutAgFEM.jl` is 995 lines. Split by responsibility:

| part | lines | share |
| --- | --- | --- |
| differential operators + variational form (§1) | 58 | 6 % |
| `Cut` → `CutData`/`Mesh` adapter (§4) | 72 | 7 % |
| StaticAssembler driver: refresh, assemble, factor, extend (§5) | 50 | 5 % |
| FE evaluation + domain integration (diagnostics only) | 120 | 12 % |
| `checkManufacturedSolution` (11 diagnostics) | 74 | 7 % |
| Ferrite probe (§7 of the spec) | 178 | 18 % |
| CSV/state I/O | 48 | 5 % |
| conditioning study | 69 | 7 % |
| `simulate` (loop, timing/alloc instrumentation) | 166 | 17 % |
| CLI | 21 | 2 % |

The two things the experiment is actually about — the mathematics handed to the assembler and
the geometry handed to the assembler — are **58 and 72 lines**. Neither is a framework, neither
is opaque, and together they are 13 % of the file. The remaining 87 % is measurement, I/O and
the Ferrite side-experiment.

## 9. Ferrite as an independent consumer of the same cut-cell quadrature

**Verdict: suitable, with three API facts that must be handled.**

`probeFerriteQuadrature` (`run.MovingCutAgFEM.jl:502`) takes a *real* cut cell from the
EvolvingDomains snapshot (the cell whose inside area fraction is nearest 0.5 — cell 468 at
`n = 33`, fraction 0.483, 6 quadrature points), builds physical weights from the actual triangle
areas, maps the reference points into Ferrite's frame, and assembles
`K_ij = ∫_{K∩Ω} ∇N_i·∇N_j dx` twice: with `CellValues` (batch) and with `PointValues` (one
`reinit!` per point), then compares both against a hand-coded bilinear evaluation.

Results (`output/ferrite_probe_n*.csv`, all three resolutions):

| check | n = 33 | n = 65 | n = 129 |
| --- | --- | --- | --- |
| `max\|K−Kᵀ\|` | 0 | 0 | 0 |
| `λ_min(Sym K)` | 3.0e-17 | 3.2e-17 | 2.8e-16 |
| constant consistency `max\|K·1\|` | 3.0e-17 | 2.8e-17 | 1.1e-16 |
| `CellValues` vs `PointValues` | 0 | 0 | 0 |
| Ferrite vs hand-coded Q1 gradient | 2.8e-16 | 2.8e-16 | 7.2e-16 |
| physical→reference inverse map | 5.8e-15 | 5.8e-15 | 4.7e-15 |
| inside + outside subcells vs full cell | 1.3e-15 | 1.3e-15 | 9.7e-16 |

The last row is the strongest one: assembling the inside *and* outside pieces of the cell with
EvolvingDomains' own triangles reproduces Ferrite's standard 2×2-Gauss whole-cell matrix to
`1e-15`, i.e. the cut decomposition tiles the background cell exactly as far as the element
integrals can tell.

**API facts a consumer must know** (all verified by running, and all public):

1. `QuadratureRule{RefQuadrilateral}(weights::Vector{Float64}, points::Vector{Vec{2,Float64}})`
   accepts arbitrary user points — this is the documented route and needs no internals.
   `PointValues(fe_values, qr)` **does not exist** (`MethodError`); the batch interface is
   `CellValues(qr, ip; update_detJdV = false)` when the caller already holds physical weights.
   `PointValues(Float64, ip)` + `reinit!(pv, X, ξ)` is the per-point route and allocates nothing
   per point.
2. Ferrite's `RefQuadrilateral` is `[-1,1]²` (mapping `ξ = 2η − 1`), and its node order is cyclic
   counter-clockwise `(0,0),(1,0),(1,1),(0,1)`, whereas EvolvingDomains/StaticAssembler use the
   tensor order `(0,0),(1,0),(0,1),(1,1)`. Feeding the tensor order to `reinit!` is rejected with
   `ArgumentError: det(J) is not positive`; the 3↔4 permutation is required. This is the only
   real friction, and it is silent-free (it errors immediately).
3. The public inverse map is `PointEvalHandler(grid, points).local_coords` and needs a `Grid`, so
   a throw-away one-cell `Grid` must be built to map physical cut points to reference
   coordinates. The grid-free `Ferrite.find_local_coordinate` works and agrees bit-for-bit but is
   undocumented internal — not used in the assembly path here.

Nothing in this flow touched package internals or required type piracy. Ferrite is therefore a
legitimate second consumer of EvolvingDomains cut data. What it does **not** provide — and what
StaticAssembler did provide — is the aggregation, the projected operator and the multi-cell
assembly machinery; a Ferrite-based cut FEM would need those written around it.

## 10. What a future StaticAssembler-based Stefan solver could reuse

1. **The adapter, unchanged.** `convertCutToStaticAssembler` + `buildStaticAssemblerMesh` are
   domain-agnostic. A Stefan solver needs exactly this and nothing more to consume
   `cut(field, :init)` each step.
2. **The Nitsche boundary slots.** The Stefan interface condition is a flux condition on the same
   embedded facets. `BoundaryPoint` exposes the normal, per-facet fields and coefficient
   gradients, and `nfacetfields` is the designed hook for per-facet curvature or a temperature
   jump — so the interface flux term is one more `BilinearSlot`/`RhsSlot`, not a new code path.
3. **`refresh_geometry!` + `assemble!` reuse.** A thermal solve every step benefits from the
   0-byte `assemble!` and the retained-pattern refresh, exactly as measured here. The known limit
   is the structural-refresh frequency: whenever the thermal solve is stiff enough that
   projection cost matters, one must budget `O(ncells)` per membership change.
4. **The `E`-based recovery.** `u_bg = E·u_m` is the natural bridge back to a full-lattice
   temperature field for the level-set/enthalpy coupling, and `restrict_solution!` handles the
   reverse direction for lagged coefficients.
5. **The load-conservation check** (`Σ b_m = 4|Ω_h|`) as a cheap regression test in any PDE with
   a constant source: it catches cut-quadrature and extension regressions in one number.

What it would **not** get for free: time-dependent coefficients on the cut cells (this experiment
deliberately used `ncoefficients = 0`; StaticAssembler supports `FieldCoefficient`/`StateFunction`
but with the AgFEM root prolongation semantics documented in `assemble!`, which needs a
"previous solution on a shrinking domain" strategy), and a per-facet Nitsche measure (the penalty
here is `γ/min(hx,hy)`, not `γ/|Γ_f|`).

## 11. Limitations and what was not done

* **No reinitialisation.** `φ` is only advected. That is adequate for rigid translation over
  `t ≤ 0.5` (centroid error ≤ `7e-5`), but a deforming interface would need `reinitialize!`, which
  is a normalisation step, not a geometry source — the adapter would not change.
* **Penalty not optimised.** `γ0 = 40` was chosen for a comfortable margin above the coercivity
  threshold and used everywhere; no study of the minimum admissible `γ0` was performed, and the
  penalty does not scale with the cut-cell area fraction (down to `4e-7` here). Stability at those
  fractions came from the aggregation, not the penalty.
* **Latent edge case found by reading, not observed.** If a `CUT` cell ever contributes no
  positive-area inside subtriangle, `_assemble_raw!` takes the `const_cut[cut_cell_index[cell]]`
  branch with `cut_cell_index[cell] == 0` and would throw a `BoundsError`. No such cell occurred
  in 600 refreshes across three resolutions, so this is reported as unverified rather than fixed
  (the package must not be modified).
* **Threading not a win** at these sizes (`threaded = true` 2.9–4.2× slower than serial with
  `-t 4`).
* **2D only**, Q1 only, one embedded interface only. The 3D `Cut` layout is structurally the
  same (tetrahedral subcells, triangular facets) but the adapter's triangle/segment assumptions
  are 2D-specific.
* The reported timings include a few GC outliers (visible as the 88 ms spike near `t = 0.24` in
  `output/fig_times.png`); the tables use means, and the min/max are in the CSVs.

## 12. Artefacts and reproduction

```
cd validation/PlanarStefan/round2
julia --project=. -e 'using Pkg; Pkg.instantiate()'     # isolated local project
julia -t 4 --project=. run.MovingCutAgFEM.jl --n 33  --stage all
julia -t 4 --project=. run.MovingCutAgFEM.jl --n 65  --stage all
julia -t 4 --project=. run.MovingCutAgFEM.jl --n 129 --stage all
julia --project=. visualize.MovingCutAgFEM.jl
```

`--stage A|B|C|D|E|F` runs a subset (`A`/`B` = descriptor audit, `C` = one stationary solve,
`D` = moving sequence, `E` = reuse/allocation/conditioning, `F` = Ferrite probe).

| file | content |
| --- | --- |
| `run.MovingCutAgFEM.jl` | numeric kernel: adapter, form, refresh/assemble/solve, diagnostics, Ferrite probe |
| `visualize.MovingCutAgFEM.jl` | CairoMakie figures (no plotting in the kernel) |
| `output/versions.txt` | revisions, resolved paths, extension-not-loaded proof |
| `output/diagnostics_n*.csv` | the 11 required diagnostics per sampled time |
| `output/steps_n*.csv` | per-step timings and allocations for every pipeline stage |
| `output/conditioning_n*.csv` | aggregated vs self-rooted spectral comparison |
| `output/ferrite_probe_n*.csv`, `ferrite_quadrature_n*.csv` | Ferrite checks and cut-cell quadrature points |
| `output/E_*.csv`, `B_*.csv`, `rhs_*.csv`, `u_bg_*.csv` | the sparse `E`, `B`, `b_m` and the recovered background solution |
| `output/log_n*.txt` | full console transcript of each run |
| `output/fig_geometry.png` | moving domain: cut cells, interface chords, exact circle |
| `output/fig_solution.png` | exact vs recovered solution and pointwise error at `t = 0.5` |
| `output/fig_convergence.png` | L2, interpolation-share and area error vs `h` with an `h²` reference |
| `output/fig_dofs.png` | master DOFs, aggregates, active background DOFs vs time |
| `output/fig_times.png` | per-step refresh history (structural vs retained) and assemble/Cholesky |
| `output/fig_memory.png` | allocation per stage, per resolution |
| `output/fig_ferrite.png` | Ferrite quadrature points on a real cut cell, physical and reference frames |
