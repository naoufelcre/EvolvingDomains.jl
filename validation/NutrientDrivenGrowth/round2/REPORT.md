# NutrientGrowthRefined — nutrient-driven free-boundary growth, corrected baseline and velocity-recovery study

Campaign II, Experiment 5.  Round-two reimplementation of `validation/NutrientDrivenGrowth`
(read-only here) on the EvolvingDomains redesign checkout, with the numerical diagnostics of that
round corrected and the interface-velocity question investigated directly.

---

## 0. Verdict

1. **The correct package was loaded; that is now proved and recorded, not assumed.** The run
   refuses to start unless `pkgdir(EvolvingDomains)` is this repository, `pkgversion ≥ 1.0`,
   `pathof` is inside that root, and the v1.0 surface is present while the pre-1.0 names are
   absent.  The provenance block (root, git revision, dirty-path count, version, native CutCells
   library path *and* SHA-256, Julia/Gridap/GridapEmbedded versions, source hashes) is printed and
   written into every `run_metadata.txt` and `output/verification_provenance.txt`.

2. **The round-one baseline is reproduced bit for bit.** `n64_A` gives
   `final_area = 0.2394063537658799` and `mean balance residual = 0.10568766448280925`, against
   `0.2394063537658799` and `0.10568766448280924` recorded in `../README.md`.  The first round
   therefore *did* execute this redesign (a 0.5.0 checkout cannot even call the API this case
   uses), so the "wrong package version" premise of this campaign is not supported by the
   evidence; what was wrong in round one was a *diagnostic*, not the environment.  The working
   tree has nevertheless changed since then (three different SHA-256 of `src/EvolvingDomains.jl`:
   `77f7…` recorded by round one, `9df7…` at `HEAD`, `fd7a…` loaded now), so identity rests on the
   runtime guard plus this numerical reproduction, not on the version string.

3. **Interface-located recovery (Method B) removes the recovery bias completely in the
   pre-merger phase.**  Averaged over `t < 0.20`, the interface flux `∫_Γ v·n ds` equals the
   growth integral `∫_Ω g(c)` to **0.1 %** with Method B on `96²…128²` (0.2 % at `64²`) against a
   4–17 % deficit with the original centroid recovery (A).  On a manufactured test with a known
   exact interface velocity B is 3–31× more accurate than A in `|V_n|` and 5–47× in the flux.

4. **The remaining area-growth error is *not* the recovery.**  On the *merged* geometry of the
   baseline run, with a smooth manufactured pressure, B still recovers the interface flux to
   0.19–0.42 %.  On that same geometry, with the actual nutrient-driven pressure, B's flux falls
   4–7 % below `∫_Ω g` after `t = 0.35`.  The deficit therefore lives in the *coupled pressure
   solution* in the post-cusp void wedge, where the two walls sit one to two cells apart with
   opposite normals and `V_n` differs by a factor of eleven within a two-cell neighbourhood, and the
   raw Nitsche boundary flux is low by the same amount as the reconstructed one.

5. **Post-merger, the more accurate method is the one with the larger error.**  Method A's
   centroid-spike artifact partially cancels that solution deficit, so for `t > 0.35` A is closer
   to `∫g` at every resolution (e.g. 0.952 vs 0.930 at `128²`).  Over the full interval B still
   wins on the time-averaged residual (`64²`: 0.075 vs 0.106; `128²`: 0.043 vs 0.053) and on the
   cumulative area ratio (`128²`: 0.953 vs 0.946).  Reporting only the time average would hide the
   trade-off; both phases are tabulated below.

6. **Everything is built on the public surface of both packages.**  No upstream file was
   modified; the composition *state = level set + independent Gridap solves + explicit adapter*
   survived, with three adapter facts documented as limits (level-set array ↔ Q1 dof ordering,
   single-valued nodal transport velocity, no interface→cell adjacency in `Cut`).

---

## 1. Verified package provenance

Printed at every start and stored in each `run_metadata.txt`:

| item | value |
|---|---|
| EvolvingDomains root | `/Users/naoufel/Documents/Scientific stack/EvolvingDomains-Redesign` |
| EvolvingDomains source | `…/src/EvolvingDomains.jl` (inside the root, checked) |
| version | `1.0.0` |
| git revision | `0809db8bf38e5d25975022a1848798f2219289c8` (`Fixing remaining things to clean for release`) |
| working tree | **60 modified/deleted paths** — dirty; identity is enforced by the runtime guard, not by the version |
| `src/EvolvingDomains.jl` SHA-256 | `fd7a176d…e7431de2` (round one recorded `77f7…`; `HEAD` holds `9df7…`) |
| `src/DataModel.jl` | `24c6a1ed…228065f` |
| `src/CutData.jl` | `57c496e3…65c3b1` |
| `src/Kinematic/Kinematic.jl` | `c5c846a5…bc53d8` |
| `src/Geometric/Geometric.jl` | `6ab8961b…a91a87cfc` |
| native CutCells library | `validation/lib/cutcells/lib/libevolvingdomains_cut_cells.dylib`, **built locally** (`scripts/build_cxx_bindings.sh`), provenance `env-override`, SHA-256 `b06f6321…23d25c` — no published artifact is used |
| Julia | `1.11.7` |
| Gridap | `0.19.7` (`~/.julia/packages/Gridap/jIgaN`, from the *stacked* default environment) |
| GridapEmbedded | `0.9.8` (`~/.julia/packages/GridapEmbedded/PIAqm`) |
| load path | `validation/env` → `~/.julia/environments/v1.11` → stdlib |

The guard also rejects the pre-1.0 checkout that is reachable from the stacked default
environment: `EvolvingDiscreteGeometry`, `CartesianGridInfo`, `AbstractVelocitySource`,
`aggregate_cut_cells` must be **absent**, and `Cut`, `CutSurface`, `cut`, `aggregateCutCells`,
`eulerianTransport!`, `lagrangianTransport!`, `WENO5Cache`, `CIPCache`, `TransportMap`,
`extendClosestPoint!` must be **present**.  Note the trap this creates for the operator: a
relative `--project=../env` from *this* directory resolves to a non-existent path
(`…/NutrientDrivenGrowth/env`), Julia then silently falls back to the default environment, and
the 0.5.0 checkout is loaded.  Every command in §12 uses the absolute project path.

---

## 2. The problem

Two discs of tissue grow, consume nutrient and merge, in `D = [0,1]²`:

| | equation | boundary |
|---|---|---|
| nutrient | `-D_c Δc + λ c = 0` in `Ω(t)` | `c = 1` on `Γ(t)` |
| growth | `g(c) = g₀ max(c − c\*, 0)` | `D_c = 0.02`, `λ = 0.05`, `g₀ = 1.5`, `c\* = 0.2` |
| pressure | `−Δp = g(c)` in `Ω(t)` | `p = 0` on `Γ(t)` |
| velocity | `v = −∇p` (Darcy) | `V_n = v·n` |
| interface | `∂_t φ + v·∇φ = 0`, `φ < 0` inside | `φ₀ = min(‖x−c₁‖, ‖x−c₂‖) − R₀` |
| marker | `∂_t q + v·∇q = 0` (intensive) | `q₀` = 0.8/0.2 Gaussians, width 0.06 |

`c₁ = (0.33, 0.50)`, `c₂ = (0.67, 0.50)`, `R₀ = 0.15`, `t ∈ [0, 0.5]`.  Nothing analytic is
imposed after `t = 0`: every step re-cuts the current level set.  The marker is a passive
intensive scalar; **no mass conservation is claimed or measured for it** (its domain grows, so its
integral is not conserved even in the continuum).

Central identity, used as the primary physical diagnostic:

```
d|Ω(t)|/dt = ∫_Ω(t) g(c) dx        (from ∇·v = g and a material boundary)
```

## 3. Software decomposition and the sequence actually executed

EvolvingDomains owns the geometry: `cartesian Grid`/`Discretization`/`Field`, the level-set field,
the native `cut` snapshot (sub-cell triangles, boundary facets, outward normals), `reinitialize!`,
`eulerianTransport!` (WENO5 + SSP-RK3) and `lagrangianTransport!` (CIP, separate source/target
buffers).  Gridap + GridapEmbedded own the two elliptic solves: `CartesianDiscreteModel` with the
same `ncells²` partition, `DiscreteGeometry` ← the same nodal level-set array, `GE.cut`,
`PHYSICAL_IN`/`ACTIVE_IN`/`EmbeddedBoundary`/`GhostSkeleton` triangulations, CutFEM with symmetric
Nitsche (`γ_D = 10/h`) and Burman–Hansbo ghost penalty (`γ_g = 0.1 h`), both scaled by the
operator's diffusion coefficient.  The two elliptic solves are separate `AffineFEOperator`
assemblies and solves; nothing is shared between them except the growth field passed explicitly.

Per step, `dynamics(levelset, backdrop, cfg)` executes exactly the prescribed sequence, with each
stage timed and its allocations counted:

1. `ED.cut` → owned `Cut` snapshot; physical area/centroid per cell, cut-cell counts,
   `aggregateCutCells` probe.
2. Interface facets: midpoints, outward normals, lengths, owning background cell (public arrays
   only).
3. Gridap embedded geometry + active space (re-cut from the same nodal array).
4. Nutrient CutFEM solve → `c`.
5. `g(c)` on the same space.
6. Pressure CutFEM solve → `p`.
7. Velocity recovery, Method A **or** B (the only difference between the two runs).
8. Extension of the interface velocity to every background node.
9. Topology (three measures).
10. Diagnostics: divergence mismatch, geometric growth integral, FE area/perimeter.
11. Loop: CIP marker transport, level-set transport, `reinitialize!` (with distinct
    before/after snapshots), then `dynamics` **again on the updated geometry** so that every saved
    frame and every metric row refers to one physical time.

`dt = min(0.4 h/‖v_h‖∞, dt_max, t_final − t)`; `0.4` is the documented conservative CFL fraction.
With `dt_max = 5·10⁻³` the CFL branch never bound (`min dt = 5·10⁻³` in every run), so the runs are
`dt_max`-limited, and the time-step sensitivity study below is what establishes that temporal error
is not the issue.

## 4. The two velocity-recovery methods

Both produce exactly one velocity per interface facet, use the **same** extension operator, the
same geometry, the same pressures and the same diagnostics.  Only the sample point and the
reconstruction differ.

**Method A — original centroid recovery (`method = A`).**  For every active background cell,
`−∇p_h` is evaluated with the bilinear derivative of the Q1 pressure *at the area centroid of that
cell's physical part* — a point up to `h/2` inside Ω, never on Γ.  Each interface facet then
inherits the velocity of its owning background cell.  This is the round-one method, kept as the
baseline, and it is described here as what it is: a *cell-centroid* gradient, not a gradient at
the facet midpoint.

**Method B — interface-aware recovery (`method = B`).**  For every interface facet:

* samples: for each fully-inside (`CUT_IN`) cell, the average of its four nodal pressure dofs at
  the cell centre (second-order accurate for a Q1 field; no interface cell, no weakly imposed
  boundary value, no ghost term);
* **visibility filter**: a sample is used only if the straight segment from the facet midpoint to
  it stays inside the phase, tested on the same bilinear level-set interpolant the cut is built
  from.  Without this, a patch on the facing side of one disc reaches *into the other disc* (the
  gap is 2.6 cells at `64²`), where `p` belongs to a different analytic branch — a measured
  failure mode, not a hypothetical one;
* fit: weighted least squares of `p ≈ a₀ + a₁X + a₂Y + a₃X² + a₄XY + a₅Y²` in coordinates scaled
  by `h`, Wendland weight `(1 − (d/r)²)²`, radius growing from `3.5h` until ≥ 8 visible samples
  (ceiling `6h`), degree dropped to linear, then to Method A, only if that fails (counted);
* solve: the weighted design matrix is column-equilibrated and solved by **QR**, *not* by normal
  equations, and the fit is rejected unless `σ_min/σ_max ≥ 10⁻³`.  This matters: with normal
  equations the same reconstruction produced condition numbers up to `10¹⁷` and interface speeds
  an order of magnitude too large (e.g. `V_n = 0.199` where the neighbours give `0.094`), which
  then poisoned the single-valued nodal extension and made the CIP marker transport abort with a
  departure point outside the box.  The QR variant has `cond ≤ 10³` on every facet of every run;
* evaluation at the facet midpoint, i.e. on Γ, giving `v_f = −∇p(x_Γ)` and `V_n = v_f·n_h`.

**Extension (identical for both).**  Zero-order closest-facet extension: each background node
takes the velocity of the nearest interface-facet midpoint, ties broken by the smallest facet
index.  Documented limitation: a single-valued nodal field cannot represent two opposing
interfaces closer than one cell; the nearest facet wins.  §11 measures the cost of this adapter
against the first round's `O(nodes × facets)` search.

## 5. Independent verification on the initial geometry

Both elliptic solvers, and both recovery methods, are scored on the initial separated discs where
closed forms exist.  Two *separate* problems are used; the second is never presented as the
solution of the coupled system.

* nutrient: `c(r) = I₀(kr)/I₀(kR₀)`, `k = √(λ/D_c) = 1.581`;
* manufactured pressure: `−Δp = 1`, `p = 0` on Γ ⟹ `p = (R₀² − r²)/4`, `V_n = R₀/2 = 0.075`,
  `∫_Γ V_n = πR₀²` per disc.  This fixed reference is what makes the recovery question answerable.

`output/verification.csv`, `output/figures/fig_convergence.png` (right panel):

| n | area err (ED / FEM) | `p` nodal max / rel-L2 | `c` nodal max / rel-L2 | `−∫_Γ ∂_n p` | Nitsche penalty | `∫_Ω 1` | A: mean `|ΔV_n|` | A: flux err | B: mean `|ΔV_n|` | B: flux err |
|---|---|---|---|---|---|---|---|---|---|---|
| 16 | 3.6e−2 / 2.9e−2 | 4.9e−4 / 6.6e−2 | 1.2e−3 / 4.9e−4 | −0.0911 | 0.0461 | 0.1372 | 2.38e−2 | 0.313 | 8.66e−3 | 6.2e−2 |
| 32 | 9.2e−3 / 7.4e−3 | 1.8e−4 / 1.6e−2 | 4.5e−4 / 1.3e−4 | −0.1144 | 0.0260 | 0.1403 | 1.28e−2 | 0.171 | **4.16e−4** | **3.6e−3** |
| 64 | 2.2e−3 / 1.8e−3 | 2.0e−5 / 1.3e−3 | 5.1e−5 / 1.1e−5 | −0.1317 | 0.0094 | 0.1411 | 4.49e−3 | 5.8e−2 | **4.29e−4** | **4.8e−3** |
| 96 | 1.0e−3 / 8.1e−4 | 9.4e−6 / 5.3e−4 | 2.3e−5 / 4.3e−6 | −0.1350 | 0.0063 | 0.1413 | 3.01e−3 | 3.9e−2 | **3.05e−4** | **3.0e−3** |
| 128 | 5.6e−4 / 4.5e−4 | 5.3e−6 / 2.6e−4 | 1.3e−5 / 2.1e−6 | −0.1366 | 0.0047 | 0.1413 | 2.24e−3 | 2.9e−2 | **3.15e−4** | **9.7e−4** |

Exact values: area `2πR₀² = 0.141372`, perimeter `4πR₀ = 1.88496`, `V_n = 0.075`,
`∫_Γ(−∇p)·n = 0.141372`.

What this shows:

* Both solves converge at a rate consistent with second order (per-doubling error ratios of
  2.5–12 across this range) and the algebraic residuals stay at `10⁻¹⁶`/`10⁻¹⁴`.
  Perimeter relative error 8.4e−3 (16²) → 1.3e−4 (128²); nutrient raw flux `0.226 / 0.284 / 0.327 /
  0.335 / 0.339` against the analytic `0.3510`, i.e. the weak boundary condition biases the *raw*
  boundary gradient low by 36 % at `16²` and 3.4 % at `128²`.
* The discrete **Nitsche identity** holds to machine precision in both problems:
  `−∫_Γ n·∇p + (γ/h)∫_Γ p = ∫_Ω 1` (e.g. `0.1372423 = 0.1372423` at `16²`, residual `1.5·10⁻¹⁵`),
  and `D_c(−∫_Γ n·∇c + (γ/h)∫_Γ(c−1)) = −λ∫_Ω c` with recorded residual `0.0` (i.e. ≤ 1e−16).
  Those are algebraic identities of the formulation, **not** evidence that the physical boundary
  flux is accurate — the raw flux column above is 3–36 % away from the exact one while the
  identity is exact.  Both are reported separately for that reason.
* Method A's error is one-signed **bias**: its mean absolute error equals its mean signed error at
  every resolution (−2.4e−2 / −1.3e−2 / −4.5e−3 / −3.0e−3 / −2.2e−3 for `16²…128²`), a 32 %
  undershoot at `16²` and 3 % at `128²`, and it is `O(h)` — consistent with sampling `∇p` at a point
  up to `h/2` inside a domain where `|∇p|` grows linearly from the boundary.
* Method B's error is 3–31× smaller in `|ΔV_n|` (10× at `64²`, 7× at `128²`) and 5–47× smaller in
  the flux (12× at `64²`, 30× at `128²`), and its mean signed error is two orders of magnitude
  smaller than A's (−5.1e−6 at `32²`, +4.2e−4 at `64²`, +1.2e−4 at `128²`).  Its residual sensitivity to the patch
  radius is small: mean `|ΔV_n|` at `64²` over `r/h ∈ {1.5, 2, 2.5, 3, 3.5, 4.5, 6}` is
  `(4.2, 4.1, 4.2, 4.1, 4.3, 3.7, 2.4)·10⁻⁴`, with no fallbacks and at most 2 linear facets; the
  default `3.5h` is not tuned.

## 6. Manufactured test on the *merged* geometry

`run.NutrientGrowthRefined.jl merged` re-reads the final level set of an existing run (no
re-evolution, no analytic geometry) and solves on it the *separate* problem

```
−Δp = 1  in Ω,   p = 0.25 − |x − (½,½)|²/4  on Γ      ⟹  −∇p = (x − (½,½))/2 exactly.
```

Two references are reported independently.  `Σ_f (−∇p_exact(x_f)·n_f) ℓ_f` uses the *discrete*
facet normals and lengths and equals half the polygon area by the divergence theorem whatever the
normal errors are (ratio to `polygon_area/2`: `1.000000` on every geometry, i.e. the facet set is a
closed polygon); `∫_Ω 1 dΩ` on the Gridap domain is the geometric representation used by the
solvers.  `output/merged_manufactured.csv`:

| geometry from | `|Ω|` ED / FEM | A: mean / max `|ΔV_n|` | A: flux err | B: mean / max `|ΔV_n|` | B: flux err | raw Nitsche flux err |
|---|---|---|---|---|---|---|
| `n64_A` (merged, `t=0.5`) | 0.23941 / 0.23943 | 4.08e−3 / 5.92e−3 | 3.4e−2 | **5.04e−4 / 1.64e−3** | **4.2e−3** | 3.9e−2 |
| `n64_B` (merged, `t=0.5`) | 0.24410 / 0.24413 | 3.97e−3 / 5.92e−3 | 3.6e−2 | **4.17e−4 / 1.13e−3** | **4.2e−3** | 4.1e−2 |
| `n96_B` (merged) | 0.24814 / 0.24816 | 2.75e−3 / 4.58e−3 | 2.2e−2 | **3.87e−4 / 1.25e−3** | **3.0e−3** | 2.6e−2 |
| `n128_B` (merged) | 0.24889 / 0.24890 | 2.05e−3 / 2.91e−3 | 1.7e−2 | **2.41e−4 / 8.63e−4** | **1.9e−3** | 1.9e−2 |

Conclusions available from this table alone:

* Method B is still 8–10× more accurate on a **merged, non-convex, post-cusp** geometry with a
  smooth exact solution: 0.19–0.42 % on the flux.  The recovery operator is therefore not what
  fails after merging.
* The raw Nitsche boundary flux is 1.9–4.1 % low even here, matching Method A's error rather than
  B's: the *boundary gradient of the CutFEM solution* is the biased quantity, and the interior
  extrapolation is what recovers the physical flux.
* The EvolvingDomains sub-cell area and the Gridap measure agree to `1·10⁻⁵` relative, so the two
  cut representations are consistent and the "ED vs Gridap geometry" hypothesis for the balance
  error is dead.

## 7. Area-growth balance: results

`ε_A^n = |Δ|Ω|/Δt − ∫_Ω g| / max(|∫_Ω g|, ε_scale)` with
`ε_scale = 10⁻³ g₀ π R₀² = 1.06·10⁻⁴`, a fixed positive scale built from the parameters so the
denominator never vanishes; it never bound (`eps_scale_binding_steps = 0` in every run, `|∫g| ≈
0.17–0.29`).  `Δ|Ω|` is the difference of reconstructed cut areas at successive synchronised
times; `∫_Ω g` is the Gridap volume quadrature, cross-checked against an independent EvolvingDomains
cut-cell quadrature that agrees to `5·10⁻⁴` relative at `64²` and `1·10⁻⁴` at `128²` (`1.8·10⁻³` at
`32²`), and `∫_Γ v_h·n ds` is
evaluated on the actual recovered facet velocities.  Three different discrete quantities, compared
as such.

Full runs, `100` steps, `dt = 5·10⁻³`, `t ∈ [0, 0.5]` (`output/*/metrics.csv`,
`output/figures/fig_balance.png`, `fig_convergence.png` left panel):

| run | mean `ε_A` | mean `|flux−∫g|/∫g` | mean `|flux−Δ|Ω|/Δt|/|Δ|Ω|/Δt|` | cumulative `Δ|Ω|/∫∫g` | `|Ω(0.5)|` | `|Γ(0.5)|` |
|---|---|---|---|---|---|---|
| `n32_A` | 0.1766 | 0.1666 | 0.0128 | 0.8234 | 0.22818 | 2.1016 |
| `n32_B` | **0.1047** | **0.0991** | 0.0208 | **0.8873** | 0.23853 | 2.3376 |
| `n64_A` | 0.1057 | 0.1003 | 0.0066 | 0.8936 | 0.23941 | 2.0739 |
| `n64_B` | **0.0748** | **0.0764** | 0.0061 | **0.9165** | 0.24410 | 2.2957 |
| `n96_A` | 0.0710 | 0.0679 | 0.0039 | 0.9279 | 0.24454 | 2.0348 |
| `n96_B` | **0.0478** | **0.0484** | 0.0044 | **0.9481** | 0.24814 | 2.1094 |
| `n128_A` | 0.0534 | 0.0508 | 0.0029 | 0.9457 | 0.24715 | 2.0190 |
| `n128_B` | **0.0433** | **0.0435** | 0.0031 | **0.9534** | 0.24889 | 2.0799 |

Split by phase (mean of the per-step ratio `∫_Γ v·n ds / ∫_Ω g`, and mean `ε_A`):

| run | pre-merge (`t<0.20`): flux/∫g | `ε_A` | post-merge (`t>0.35`): flux/∫g | `ε_A` | post: `Δ|Ω|/Δt / ∫g` |
|---|---|---|---|---|---|
| `n32_A` | 0.8306 | 0.1783 | 0.8310 | 0.1797 | 0.8203 |
| `n32_B` | **0.9727** | **0.0337** | 0.8146 | 0.2064 | 0.7936 |
| `n64_A` | 0.9053 | 0.0995 | 0.8966 | 0.1089 | 0.8911 |
| `n64_B` | **0.9977** | **0.0040** | 0.8444 | 0.1570 | 0.8430 |
| `n96_A` | 0.9437 | 0.0583 | 0.9321 | 0.0711 | 0.9289 |
| `n96_B` | **0.9996** | **0.0031** | 0.9202 | 0.0814 | 0.9186 |
| `n128_A` | 0.9616 | 0.0393 | **0.9518** | **0.0505** | 0.9495 |
| `n128_B` | **1.0006** | **0.0025** | 0.9299 | 0.0722 | 0.9278 |

Fitted convergence of the time-averaged residual over `32²→128²`: A `0.177→0.053`, slope `h^0.86`;
B `0.105→0.043`, slope `h^0.64`.  Method B's rate is set by the post-merger phase, which is
first-order-ish and much slower than the pre-merger behaviour (`ε_A ≈ 0.003` flat in `h`).

Time-step sensitivity at `64²` (same grid, same physics):

| `dt_max` | steps | `n64_A` mean `ε_A` | `n64_A` `|Ω(0.5)|` | `n64_B` mean `ε_A` | `n64_B` `|Ω(0.5)|` |
|---|---|---|---|---|---|
| `5·10⁻³` | 100 | 0.1057 | 0.23941 | 0.0748 | 0.24410 |
| `2.5·10⁻³` | 200 | 0.1069 | 0.23940 | — | — |
| `1·10⁻³` | 500 | 0.1078 | 0.23937 | 0.0766 | 0.24403 |

Refining `dt` by 5 changes `ε_A` by 2 % relative and the final area by `4·10⁻⁵` relative.  The
residual is spatial, as in round one.

Three independent statements about the same step, `n64_B`, `t = 0.5`:

| quantity | value | ratio to `∫_Ω g = 0.291` |
|---|---|---|
| `∫_Γ v_h·n ds` (recovered facet velocities) | 0.278 | 0.956 |
| `Δ|Ω|/Δt` (reconstructed areas) | 0.277 | 0.952 |
| `−∫_Γ n·∇p_h` (raw CutFEM boundary flux) | 0.246 | 0.846 |
| `(γ/h)∫_Γ p_h` (Nitsche penalty) | 0.045 | 0.154 |
| `−∫_Γ n·∇p_h + (γ/h)∫_Γ p_h` (discrete identity) | 0.291 | 1.000 (by construction) |

The recovered velocity and the reconstructed area agree to 0.4 % at this instant (`0.3–2 %` averaged
over all runs): geometry transport and time integration faithfully realise the velocity they are
given.  The raw Nitsche flux is 15 % low, and the weak identity shows the difference is exactly the
boundary-penalty term — which is why the identity's near-zero residual says nothing about
`V_n`.

## 8. Attribution of the remaining error

Ordered by what the evidence excludes:

| candidate | verdict | evidence |
|---|---|---|
| recovery *location* (the A bias) | **resolved**, was 4–17 % | A's error is `O(h)` bias; B has no systematic sign and is 0.1–0.2 % accurate pre-merge (§5, §7) |
| recovery *operator* | excluded post-merge | on the identical merged geometry B is 0.19–0.42 % accurate with a smooth field (§6) |
| geometry representation (ED cut ↔ Gridap cut) | excluded | sub-cell area vs FE measure agree to `1·10⁻⁵`; facet polygon closed to 1e−16 relative (§6) |
| geometry transport / time integration | excluded | `|flux − Δ|Ω|/Δt| / |Δ|Ω|/Δt|` = 0.3–2 % everywhere; `dt` sensitivity 2 % (§7) |
| diagnostic quadrature (`∫g`) | excluded | Gridap volume quadrature vs independent cut-cell quadrature agree to `5·10⁻⁴` (`64²`) and `1·10⁻⁴` (`128²`) |
| **coupled pressure solution in the post-cusp void wedge** | **remaining** | see below |

The post-merger geometry is a peanut whose two *cusps* — the points where the surfaces of the two
lobes meet — leave a void wedge one to two cells wide.  Within two cells of `x = ½` on `n128_B` at
`t = 0.5` there are 22 interface facets in near-coincident pairs at `x = 0.500 ± 0.001` with
*opposite* normals, and `V_n` on them ranges from `0.020` to `0.225` — a factor of eleven inside a
two-cell neighbourhood, on walls that are one to two cells apart.  That is a thin-wedge boundary
layer of the pressure problem, and its discrete normal derivative is where the deficit appears: at
`t = 0.5` the reconstructed flux and the *raw* Nitsche flux agree to 2.4 % (`0.826` vs `0.846` of
`∫g` at `64²`), i.e. both interior extrapolation and boundary gradient see the same deficit, while
on the same geometry with a smooth field both are accurate.  Consistent with that reading, the
deficit decreases with `h` (`0.844 → 0.920 → 0.930` for B post-merge) and the maximum `|v_h|` grows
(`0.12 → 0.24 → 0.29`), i.e. the peak flux is being resolved more sharply.

Two mechanisms cannot be separated by this experiment and both act in the same place: (i) the
CutFEM solution of `−Δp = g` in the thin, nearly closed void wedge (a Dirichlet-only channel a cell
wide), and (ii) the single-valued transport velocity, which must give both walls of the wedge one
value.  The transport-consistency diagnostic (`|flux − Δ|Ω|/Δt| = 0.3–2 %`) shows the extension
faithfully realises whatever single value it chose, so the ambiguity costs little *in the recorded
quantities*; but it does mean the wedge's two walls cannot both move at their own `V_n`, which is
exactly where the missing flux is.

Two honest qualifications.  First, this is inference from a combination of measurements, not a
proof: the decisive extra experiment would be a manufactured solution with a cusp-like singularity
on the same geometry, which was not built.  Second, the *shapes* after merging are not converged at
`64²`: B's final perimeter is `2.34 / 2.30 / 2.11 / 2.08` at `32/64/96/128²` while A's is
`2.10 / 2.07 / 2.03 / 2.02`, so at coarse resolution the two methods deliver genuinely different
merged geometry, and the post-merge residual is a functional of that geometry as much as of the
solve.

## 9. Topology

Three measures are tracked every step, and they are *not* the same statement
(`output/figures/fig_topology.png`):

| measure | definition | first single-component time, A | B |
|---|---|---|---|
| `components_nodes` | 4-connected set of grid nodes with `φ < 0` | 0.305 / 0.340 / 0.300 / 0.270 | 0.215 / 0.215 / 0.215 / 0.215 |
| `components_geometric` | edge-connected components of the `PHYSICAL_IN` cut triangulation (vertices welded at `h·10⁻⁹`; triangles *and* quads, so sub-cells and whole cells are handled together) | 0.305 / 0.340 / 0.300 / 0.270 | 0.215 / 0.215 / 0.215 / 0.215 |
| `components_cells` | 4-connected set of active background cells (`state ≠ CUT_OUT`) | `1` from `t = 0` at `32²`; 0.055 / 0.110 / 0.140 | `1` from `t = 0` at `32²`; 0.055 / 0.110 / 0.135 |

(resolutions `32² / 64² / 96² / 128²`.  Times are *post*-step: the first step whose updated geometry
is already connected, from `run_metadata.txt`.  The `components_*` columns of `metrics.csv` are the
counts on the geometry at `t_start`, i.e. the geometry that row's fields belong to, so they lag by
one step; that is the convention of the whole row, not a second measurement.)

* The `cells` measure is **degenerate**: at `32²` the two cut bands are already 4-connected at
  `t = 0` (the gap is 1.3 cells), and at `64²…128²` it flips at `0.055…0.140`, long before the
  material connects, because a cut band one cell thick on each side already joins the two regions.
  It must not be reported as a merger time — round one quoted it as one, and it is the reason its
  merger times looked resolution-consistent.
* Method B's merger time is `0.215` at every resolution.  This is *not* a convergence result: for a
  symmetric configuration with an even grid there is a grid node exactly at `(½, ½)`, where
  `min(d₁,d₂) − R = gap/2 − R`, so the node-sign criterion flips precisely at the continuum contact
  time — which for a nearly uniform `V_n ≈ 0.09` and an initial gap of `0.04` is `≈ 0.22`.  The
  agreement is a property of the measure, not of the discretisation, and no claim of a converged
  merger time is made.
* Method A's merger time is later (`0.270–0.340`) and **non-monotone in `h`**.  Its velocity is
  4–6 % low and its near-cusp velocities are noisy, so the last sub-cell of the gap closes at a
  resolution-dependent time.  Topology timing with the biased velocity is not converged, and no
  single-resolution merger time should be quoted from it.

## 10. Corrected diagnostics

**Reinitialisation, measured with two distinct snapshots.**  The first round returned the
*pre*-`reinitialize!` snapshot as the post-reinitialisation measurement, so its reported change was
identically zero by construction.  Here `pre = ED.cut(...)` is taken after transport, φ is copied,
`reinitialize!(levelset, pre)` runs, and `post = ED.cut(...)` is a fresh snapshot.  Over all 11 runs:

| observable | max over all steps and runs |
|---|---|
| `|Δ|Ω||` across reinitialisation | `2.8·10⁻¹⁷` (round-off) |
| one-sided interface displacement, pre → post facet midpoints: mean / max | `5.8·10⁻¹⁸` / `1.1·10⁻¹⁶` |
| sign flips of `φ` | `0` |
| RMS change of `φ` itself | `4.9·10⁻³` (32²) → `1.4·10⁻³` (128²) |

So `reinitialize!` does not move this interface at all — it preserves the lattice crossing fractions
exactly — while it does restore the distance property off the interface.  Round one's conclusion
was therefore right by accident, and the corrected measurement confirms it quantitatively and adds
the field-level change that makes it a meaningful statement.

**Snapshot synchronisation.**  The loop evaluates `dynamics` again on the updated geometry, so a
saved frame's `c`, `p`, `g`, `v` and marker are all solved on that frame's own `φ` and time; the
metrics row of step `n` carries `t_start`, `t_end` and both areas, and the residual compares the
area difference with the growth integral measured at `t_start`.  No frame mixes fields from one
time with geometry from another, so no staggered label is needed.  Frames are stored every 10th
step plus the last (`frame_count = 11`).

**Velocity-recovery documentation.**  Method A's gradient is evaluated at the physical-part
centroid of each background cell, up to `h/2` inside Ω; the facet inherits its owning cell's value.
Method B's gradient is a local polynomial reconstruction evaluated *at the facet midpoint on Γ*.
The code and the CSV (`plot_*` fields, `order` strings, §4) say which is which; no cell-centroid
gradient is described as an interface gradient anywhere in this case.

**Other corrections relative to round one.**  The nutrient diagnostic now reports the raw boundary
flux, the Nitsche penalty term and the discrete identity separately instead of treating
`|D_c∫_Γ ∂_n c − λ∫_Ω c|` as a solver error (it is a Nitsche effect); the growth integral is
cross-checked with an independent quadrature; the topology count gained the geometric measure; the
extension search is validated against brute force; every stage is timed and its allocations
counted.

## 11. Cost: package kernels versus application adapters

Medians over steps ≥ 3 (compilation excluded), `output/profile.csv`, `fig_performance.png`:

| stage | `n64_A` | `n64_B` | `n128_A` | `n128_B` | allocation (64²/128²) |
|---|---|---|---|---|---|
| Gridap embedded re-cut (`GE.cut` + geometry) | 18.3 ms | 18.4 ms | 20.9 ms | 20.8 ms | 21 / 26 MiB |
| pressure CutFEM (assembly + solve) | 16.0 | 15.6 | 29.4 | 29.7 | 36 / 79 MiB |
| nutrient CutFEM (assembly + solve) | 14.6 | 15.2 | 28.4 | 28.0 | 34 / 75 MiB |
| **velocity recovery** (A: centroid+extension, B: patch+extension) | **0.50** | **2.49** | **3.59** | **7.60** | 0.4 / 1.1 (A), 2.5 / 5.9 (B) MiB |
| topology diagnostics (3 measures) | 1.10 | 1.09 | 3.86 | 3.85 | 2.2 / 6.5 MiB |
| EvolvingDomains `cut` | 0.70 | 0.73 | 1.67 | 1.67 | 2.6 / 11 MiB |
| Gridap space construction | 0.63 | 0.63 | 1.46 | 1.48 | 1.1 / 3.3 MiB |
| growth rate field | 0.44 | 0.45 | 0.93 | 0.94 | 0.8 / 2.4 MiB |
| cut-cell geometry (area/centroids from `Cut`) | 0.22 | 0.22 | 0.54 | 0.54 | 0.2 / 0.8 MiB |
| everything else (area, perimeter, facets, divergence, …) | < 0.2 | < 0.2 | < 0.4 | < 0.4 | — |
| **sum of stage medians** | **52.7 ms** | **55.1 ms** | **91.1 ms** | **94.9 ms** | 99 / 206 MiB |

Outside `dynamics` (from `metrics.csv`): level-set transport + `reinitialize!` 2.9 ms (64²) /
9.5 ms (128²), CIP marker transport 0.15 ms / 0.55 ms.  Full step 60 ms (64²) / 110 ms (128²);
whole `t ∈ [0, 0.5]` run 12 s (`64²`) and 18 s (`128²`) including compilation.

Reading:

* The package kernels are *not* the cost.  `ED.cut` is 0.7 ms of a 60 ms step (1 %);
  `eulerianTransport!` + `reinitialize!` are 2.9 ms (5 %).  Gridap's embedder re-cuts the whole
  background model every step and costs 18 ms at `64²` — and only 21 ms at `128²`, i.e. it is
  dominated by the *setup* of `DiscreteGeometry`/`EmbeddedDiscretization`, not by the cell count.
  That is the single largest fixed cost of the composition, and it is an upstream (GridapEmbedded)
  property, not something this case can fix.
* Method B costs 2.0 ms more than A at `64²` (2.49 vs 0.50 ms) and 4.0 ms more at `128²` — 4–5 %
  of a step, for a 3–31× accuracy gain in `|V_n|` in the pre-merger phase.  The extra cost is the visibility
  filter (≈ 6 segment samples × candidates per facet) plus the 6×6 fits.
* **The interface→grid extension is not a bottleneck, and indexing it does not pay.**  Measured
  (`run.NutrientGrowthRefined.jl extension 64 128 192`, `output/extension_benchmark.csv`):

| n | facets | nodes | brute force (used) | CSR + expanding-ring index | evaluation ratio | `max|Δv|` between them |
|---|---|---|---|---|---|---|
| 64 | 152 | 4 225 | **0.41 ms** | 4.8 ms | 9.4× fewer | 0 |
| 128 | 308 | 16 641 | **3.16 ms** | 63 ms | 11.9× fewer | 0 |
| 192 | 460 | 37 249 | **10.5 ms** | 286 ms | 13.3× fewer | 0 |

  The index visits 9–13× fewer distance evaluations and is 3–20× *slower* in wall time: the dense
  sweep is a contiguous SIMD-friendly loop (0.6 ns per evaluation) while the ring scan spends
  ~60 ns per visit on bounds checks and bookkeeping.  Both are kept: the pipeline uses brute force
  (`shortcut:` comment in `nearestFacetExtension`), the indexed variant is retained only as the
  measured comparison, and the two are asserted identical (`max|Δv| = 0`, ties broken by the
  smallest facet index so that "identical" means identical *velocities*, not merely identical
  distances).  Round one's `O(nodes × facets)` search was therefore not the inefficiency it was
  suspected to be, at least up to `192²`.
* Adapter versus kernel is separated in the table: the application-level adapter
  (facet extraction, recovery, extension, topology) is 1.8 ms (A) / 3.8 ms (B) at `64²` and 8.0 /
  12.0 ms at `128²`, i.e. 3–13 % of the step; the rest is packages.

## 12. The composition survived

The architecture of `intended.md` — *geometry is a level set field, state is a field, a plain
`dynamics` function returns the rate, the loop transports geometry and state* — holds with two
independent discretisation libraries and needed **no modification to either package**:

* EvolvingDomains is used only through `Grid`, `Discretization`, `Field`, `cut`, `Cut`,
  `CutSurface`, `CUT_IN/CUT_OUT/CUT_CUT`, `aggregateCutCells`, `reinitialize!`,
  `eulerianTransport!`, `WENO5Cache`, `lagrangianTransport!`, `CIPCache`.
* Gridap/GridapEmbedded only through `CartesianDiscreteModel`, `FESpace`/`FEFunction`,
  `DiscreteGeometry`, `cut`, `Triangulation(..., PHYSICAL_IN/ACTIVE_IN)`, `EmbeddedBoundary`,
  `GhostSkeleton`, `Measure`, `get_normal_vector`, `AffineFEOperator`, `get_matrix`, `get_vector`,
  `interpolate`, `get_node_coordinates`, `get_cell_coordinates`.
* All recovery, extension, topology and diagnostics are application-level functions in one file;
  no new solver framework wraps the sequence.

Documented limits of the seam (unchanged from round one, now with measured consequences):

* `Grid`/`Field` arrays carry no metadata: the Q1 dof order of Gridap's model equals `vec` of the
  `(nx, nx)` node array (first axis fastest).  This is asserted at start-up
  (`checkNodeOrdering`, max deviation ≤ 100 eps) and again in `backgroundNodal` (bijection checked);
  it would raise rather than degrade silently, but it is still an ordering coincidence of two
  packages.
* `Cut` is surface geometry: sub-cells are not a triangulation with dofs, there is no interface →
  cell adjacency beyond `facet_to_bgcell`, and Gridap must re-cut from the nodal array every step.
* The transport operators accept a single-valued nodal velocity, which cannot represent two
  opposing interfaces closer than a cell.  This is the one place where the composition, not the
  numerics, limits the post-merger phase: at `t = 0.5` the two walls of the post-cusp void wedge sit
  one to two cells apart with opposite normals, and the nearest-facet rule necessarily gives both
  walls one value.

## 13. Remaining uncertainties

1. **The post-merger deficit is attributed, not proved.**  A manufactured solution with a
   cusp-like singularity on the merged geometry would separate the CutFEM thin-wedge boundary layer
   from the single-valued-extension ambiguity.  Until then, "the coupled pressure solution in the
   post-cusp void wedge" is the best-supported of several surviving explanations.
2. **Post-merger shapes are not converged at `64²`.**  The final perimeter of method B still changes
   by 9 % between `64²` and `128²` (`2.296 → 2.080`; method A: 2.6 %), and the two methods deliver
   visibly different merged shapes at coarse resolution; any post-merger metric inherits that.
3. **The merger time is not a converged quantity.**  `0.215` for B at four resolutions is a property
   of the symmetric node-sign criterion (see §9), not a discretisation result.
4. **The CIP marker is not conserved and is not checked against anything.**  It is transported with
   the same single-valued velocity, so it inherits the same limitation; its range (`−1.6·10⁻²⁷ …
   0.801`) is reported only as a boundedness indicator.
5. **The divergence mismatch is not method-specific.**  `divergenceMismatch` differentiates the
   *cell-centroid* field, which both methods share (B reuses A's cell samples), so the columns
   differ only because the two trajectories differ.  Its value (1–8 %) does not converge cleanly and
   it is reported as a diagnostic of the centroid field, not as a method comparison.
6. **The working tree is dirty (60 paths) and `src/EvolvingDomains.jl` has changed since round one.**
   The numerics of the A path are bit-identical to round one, which bounds what those changes can
   have touched, but this is evidence, not a guarantee.
7. **Not attempted:** `aggregateCutCells`-based AgFEM stabilisation (the Nitsche + ghost-penalty
   formulation stayed stable down to inside-fractions of `10⁻⁹`); curvature-driven `p = γκ`
   (no curvature API in either environment); grids finer than `128²`; 3D.

## 14. Reproducing

```sh
cd "…/EvolvingDomains-Redesign/validation/NutrientDrivenGrowth/round2"
PROJ="/Users/naoufel/Documents/Scientific stack/EvolvingDomains-Redesign/validation/env"
export EVOLVINGDOMAINS_CUTCELLS_LIBRARY="/Users/naoufel/Documents/Scientific stack/EvolvingDomains-Redesign/validation/lib/cutcells/lib/libevolvingdomains_cut_cells.dylib"

# elliptic verification + recovery scoring on the initial geometry (16²…128²)
env "$EV" julia --project="$PROJ" --startup-file=no run.NutrientGrowthRefined.jl verify

# one coupled run: ncells tfinal method(A|B) dtmax tag
env "$EV" julia --project="$PROJ" --startup-file=no run.NutrientGrowthRefined.jl run 64 0.5 B 5.0e-3 n64_B

# manufactured test on the merged geometry of a finished run
env "$EV" julia --project="$PROJ" --startup-file=no run.NutrientGrowthRefined.jl merged n64_A n64_B n96_B n128_B

# extension adapter: brute force vs CSR ring index
env "$EV" julia --project="$PROJ" --startup-file=no run.NutrientGrowthRefined.jl extension 64 128 192

# whole matrix (11 runs, ~10 min) and the figures
bash run_suite.sh
env "$EV" julia --project="$PROJ" --startup-file=no visualize.NutrientGrowthRefined.jl
```

Outputs (`output/`): `verification.csv`, `verification_provenance.txt`,
`merged_manufactured.csv`, `extension_benchmark.csv`, one directory per run with `metrics.csv`
(one row per step), `profile.csv` (one row per stage per step), `frames.jls` (11 synchronised
frames), `verification.txt` and `run_metadata.txt` (provenance + all aggregates), and
`figures/*.png`:

| figure | content |
|---|---|
| `fig_domain_evolution.png` | 8-frame montage of `Ω(t)` with cut-facet interface, area and the three component counts |
| `fig_fields.png` | nutrient `c` and pressure `p` on synchronised geometry, 3 frames |
| `fig_velocity.png` | cell speed, nodal velocity arrows, and `V_n` on the interface |
| `fig_marker.png` | passive marker `q` with interface, 3 frames |
| `fig_balance.png` | `d|Ω|/dt`, `∫g`, `∫_Γ v·n` for A and B at `64²`, plus both residuals |
| `fig_area.png` | physical (and FE) area versus time, all runs |
| `fig_topology.png` | the three component counts versus time, A versus B |
| `fig_convergence.png` | time-averaged `ε_A` versus `h` for both methods; verification recovery errors and fluxes versus `h` |
| `fig_performance.png` | per-stage time and allocations per step |

The interface in the montage is drawn from cut-facet midpoints rather than from a `φ = 0` contour:
the `min()` of two distances is non-smooth on the medial line, so the *bilinear* interpolant of `φ`
has a spurious zero crossing in the void there.  That crossing produces no cut facet and never
enters the computation, but it does appear in a contour plot.
