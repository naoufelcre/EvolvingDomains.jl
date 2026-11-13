# CurvatureExtension3D — Curvature.jl-driven three-dimensional mean-curvature flow

**Case directory:** `validation/MeanCurvatureSphere3D/round2/`
**Deliverables:** `run.CurvatureExtension3D.jl`, `visualize.CurvatureExtension3D.jl`, `REPORT.md`,
`output/` (`n32/`, `n40/`, `stageA/`, `stageC/`)
**Model:** DeepSeek V4.1 Flash. Delegation was attempted once and impossible in this environment
(`subagent` returned `spawn pi ENOENT`); everything below was produced by the single main agent.

Every number in this report is read from the artifacts this run wrote (`timeseries.csv`,
`snapshots.csv`, `cost.csv`, `reinit.csv`, `stageC/extension_test.csv`) or was measured by the
probe in `output/stageA/probe.jl`. Claims are labelled **measured**, **analytical**, **code fact**,
**hypothesis**, or **open**.

---

## 1. Executive summary

Curvature.jl now **drives** the interface. At every one of the 14 (32³) and 22 (40³) actual time
steps, the pipeline is

```text
cut(field) -> CutSurface -> facet vertices/normals -> Curvature.jl interfaceCurvature
           -> V_f = -mu*kappa_f -> closest-point extension of V_f * n_f to band nodes
           -> eulerianTransport! (WENO5 + SSP-RK3)
```

with `kappa` never taken from a finite-difference operator in the `cv` run, and no analytic
trajectory injected anywhere (the exact radius `R(t) = sqrt(R0^2 - 4 mu t)` appears only in
diagnostics).

**Sphere control (A), measured.** At 40³ the Curvature.jl-driven interface shrinks by
0.0106846 against the exact 0.0106829 — a **+0.01 % displacement error** — and ends at
`R = 0.218766` against the exact `0.219317` (−0.251 %), where −0.239 % of that end error is already
present in the initial piecewise-linear cut. The independent nodal-finite-difference run of the same
setup gives +0.51 % displacement error at the same resolution. At 32³ the two are −0.50 % (Curvature.jl)
and +0.59 % (FD).

**Nonspherical case (B), measured.** The P2-perturbed sphere evolves only through its computed
curvature: mean radius 0.231264 → 0.220414 (32³), area 0.67328 → 0.61171 (monotone decrease, −9.14 %),
polar radius 0.263120 → 0.250204, equatorial 0.211822 → 0.202916, so the anisotropy falls from
0.051298 to 0.047288 (−7.8 %) while the interface stays smooth (sphericity 6.85 × 10⁻² → 6.58 × 10⁻²).
The two curvature operators agree to 8.6 × 10⁻⁵ in mean radius (0.04 %) after 14 steps, and their
trajectories do not separate.

**Cost, measured.** Curvature.jl is 95 % of the wall clock: 1623 ms/step at 32³ and 2069–2444 ms/step
at 40³, against **17–28 ms/step** for the whole facet-to-grid extension. The bridge is 89–95× cheaper
than the curvature it transfers.

### The five required answers

1. **Does an appropriate public extension operator exist?** No. Neither package, today or in the
   reachable history, exposes a facet→background-node transfer. `extendClosestPoint!`
   (`src/Geometric/ClosestPointExtension.jl:31`, exported) is nodal→nodal along `∇φ`; the deleted
   `src/Transfer/ExtensionOperator.jl` was the same operation in 2D with Gridap types. Curvature.jl
   offers only `curvatureAt` (a pointwise *re-fit*, not a value transfer). Details and citations in
   `output/stageA/AUDIT.md`.
2. **Can Curvature.jl drive EvolvingDomains without upstream modifications?** Yes, and it did. The
   only glue is `CutSurface` → plain `NTuple{3,SVector{3,Float64}}` facets plus `interfaceCurvature`
   (one call per step), and an application-local extension written in ordinary arrays. This experiment
   wrote only under `round2/`: the 59 pre-existing working-tree entries outside `validation/` are
   unchanged, no manifest or `Project.toml` was touched, and a later rerun of the 32³ case reproduced
   its radii bit-for-bit.
3. **How much code is genuinely required to bridge facet and nodal representations?** 122
   non-comment lines: facet CSR index over cut cells (22), exact closest-point-on-triangle (31),
   narrow-band nearest-facet search (42), velocity assembly with taper (27). The rest of the 970-line
   script is diagnostics, timing, CSV output and self-checks.
4. **Curvature versus extension runtime?** Curvature.jl 95.1–95.7 % of the step, extension 1.0–1.1 %;
   ratio 89–95×. Cut reconstruction is 2.7–3.4 %, WENO5 transport 0.3–0.5 %, reinitialisation 0.2 %.
5. **Does the perturbed geometry evolve consistently?** Yes, within the accuracy of the operators:
   smooth radius/area/volume histories, monotone area decrease, decaying anisotropy, no non-finite
   curvature, no unresolved band node, and CV-vs-FD trajectories that differ by 0.04 % in mean radius
   at both resolutions.

---

## 2. Verified identities (Stage A prerequisite)

| Item | Value |
|---|---|
| EvolvingDomains loaded from | `/Users/naoufel/Documents/Scientific stack/EvolvingDomains-Redesign/src/EvolvingDomains.jl` |
| EvolvingDomains revision | git `0809db8bf38e5d25975022a1848798f2219289c8`, **plus** the uncommitted redesign working tree (that tree is the subject under test; it was not modified) |
| Curvature.jl loaded from | `/Users/naoufel/Documents/Scientific stack/Curvature.jl` (source checkout on `LOAD_PATH`, no `.git`, version 0.1.0) |
| Environment | `validation/env` (`EvolvingDomains = {path = "../.."}`), Julia 1.11.7 |
| Native cut library | `$TMPDIR/evolvingdomains-cutcells-install/lib/libevolvingdomains_cut_cells.dylib` |
| Public Curvature.jl call used | `interfaceCurvature(facets, normals; spacing, origin, radius) -> Vector{Float64}` (`src/Curvature.jl:65`); exports at `src/Curvature.jl:52` |

Both facts were checked before execution: `pathof(EvolvingDomains)` resolves through `validation/env`'s
path entry to the directory above, and the run scripts discover the Curvature.jl checkout by walking up
from the case directory (`locateCurvatureSource`), so a different checkout cannot be loaded silently.

One observation about the checkout, disclosed because it bears on reproducibility rather than on the
results: during this session an outside process rewrote `src/EvolvingDomains.jl`,
`src/Geometric/Geometric.jl`, `src/Geometric/ClosestPointExtension.jl`, `README.md` and
`test/contracts/*` (file mtimes 02:21–02:22, i.e. after the session began, while this experiment was
writing only under `round2/`). Line counts and the exported surface are the same as read at the start of
the session (`EvolvingDomains.jl` 51 lines, `Geometric.jl` 17, `ClosestPointExtension.jl` 61), and the
32³ case was rerun afterwards with bit-identical radii, so the package behaviour under test did not
change across the reported runs. The revisions in the table are therefore the ones that produced the
numbers below.

---

## 3. Model, level-set convention, and the sign chain

```text
D = [-0.5, 0.5]^3,  R0 = 0.23,  mu = 0.01
phi < 0 inside,  n = grad(phi)/|grad(phi)| outward,  kappa = div n  (>0 convex)
V_n = -mu kappa        (analytical definition)
d_t phi - mu kappa |grad phi| = 0
R(t) = sqrt(R0^2 - 4 mu t)   ->  R(0.12) = 0.2193172,  kappa_exact(0.12) = 2/R = 9.11922
```

Geometry A: `phi0 = |x| - R0` (exact signed distance).
Geometry B: `phi0 = r - R0 (1 + eps P2(z/r))`, `P2(s) = (3s² - 1)/2`, `eps = 0.15`.

At `r = 0` the direction `z/r` is undefined; every directional limit of `phi0` there is negative
(between `-R0(1 + eps)` on the axis and `-R0(1 - eps/2)` in the equatorial limit), so the sign — and
therefore the interface — is unaffected. The stored value is the angular average
`-R0(1 - eps/2)`, which keeps the interior field finite with no division by zero (code fact,
`perturbedField`).

**Geometry B is not a signed-distance function.** Its gradient norm is 1.0091 at 32³ and 1.0107 at
40³ (measured, band mean; geometry A reads 0.9962 / 0.9975). Two consequences, both handled
explicitly:

- the metric of the field and of the represented surface differ, so the *cut* (which only depends on
  the sign pattern and the vertex interpolation) is correct while `|∇φ| ≠ 1`; the velocity is applied
  through `v·∇φ`, so the interface speed is `V_f (n_f·n_φ)` and does not need `|∇φ| = 1`;
- `reinitialize!` is applied every 10 steps and measured before/after with independent cuts (§8).

The sign chain that makes the sphere contract is a code fact in one line:
`rate = -cfg.mu * kappa[q]`, then `v = chi * rate * n_f` with `n_f` the outward facet normal, so a
convex facet (`kappa > 0`) moves opposite to `n`, i.e. inward.

---

## 4. What runs at every time step

```text
per step, method = :cv
  geo      = interfaceData(field, x)          # cut -> CutSurface -> facets, normals, areas, cellof
  kappa    = facetCurvature(geo, x, cfg)      # Curvature.jl, one value per facet, radius = 4h
  index    = FacetIndex(geo.cellof, ncells)   # CSR of facets grouped by cut cell
  band     = assignNearestFacets!(...)        # exact point-to-triangle nearest facet per band node
  velocity = facetVelocity!(...)              # v = chi(|phi|) * (-mu*kappa_f) * n_f
  eulerianTransport!(field, velocity, dt)     # WENO5 + SSP-RK3, first axis fastest
```

### 4.1 The extension, stated exactly

```text
v_h(x) = chi(|phi(x)|) * (-mu * kappa_{f(x)}) * n_{f(x)}
f(x)   = argmin_f  dist(x, triangle_f)        (exact closed-form point-to-triangle distance)
chi(d) = 1 for d <= 1.5h,  cos²(pi/2 (d-1.5h)/(1.5h)) for 1.5h < d < 3h,  0 for d >= 3h
```

- **Search is indexed, never all-pairs.** Candidates come from the CSR block of background cells
  within 2 cells of the node's own cell (widened once to 6 cells if that block is empty). Cost is
  proportional to band nodes × facets of a few cells.
- **Narrow band.** Nodes with `|phi| > 3h` get exactly zero velocity; nothing is extrapolated to the
  far field. The zero level set lies in the full-speed zone, so the interface speed is the facet speed
  up to the angle between `n_f` and `∇φ`.
- **Piecewise-constant by construction.** The extended field is constant in the normal direction *of
  its own facet*, so it is O(h)-smooth across facet boundaries. This is the price of using facet data;
  it is visible in §6 as the difference between the scalar and vector tests.
- **Non-finite curvature is never hidden.** A facet with `NaN` curvature contributes zero velocity at
  its nodes and increments a counter that is written to `timeseries.csv`
  (`nonfinite_velocity`). Measured: **0 at every step of every run** (all four configurations), and 0
  `NaN` facet curvatures at every snapshot.
- **Truncation of the extension is measured, not assumed**: 0 band nodes without a facet in the
  searched block in any run (`band_nodes_unreached` column).

### 4.2 Time stepping

`dt = C h²/μ` with `C = 1/12` (documented factor-two margin over the `μ dt/h² ≤ 1/6` diffusion scale;
**not** a rigorous bound for the assembled WENO5 + extension + reinit algorithm). This gives
`dt = 8.571 × 10⁻³` (32³, 14 steps) and `5.455 × 10⁻³` (40³, 22 steps), each with the last step
shortened to land exactly on `tmax = 0.12`. Per-step interface displacement is ≈ 0.023 h, far inside
the 3h band, so the taper never lags the front (`shortcut:` note in `facetVelocity!`).

### 4.3 The comparison method (`fd`)

An independent copy of the level set, advanced with the same grid, same `dt` sequence, same transport
kernel, same band, same taper and the same reinitialisation schedule; only the curvature operator
differs: nodal normals by central differences (`eps_gradient = 0.5`) and `kappa_h = div_h n_h`. This
is round 1's operator, kept unchanged so the comparison isolates the curvature source.

---

## 5. Extension-operator investigation (summary; full text in `output/stageA/AUDIT.md`)

Four operations that are routinely confused, and what exists publicly:

| Operation | Public API | Verdict |
|---|---|---|
| AgFEM master/slave DOF extension | none | **absent**; different semantics anyway (copies FE DOFs across a skeleton during assembly) |
| Closest-point extension of interface **values** | `extendClosestPoint!` (`src/Geometric/ClosestPointExtension.jl:31`, exported) | **exists but nodal→nodal**, `u_ext(x) = u(x − φ∇φ/|∇φ|)`; no facet argument |
| Normal extension with `n·∇V = 0` | none | **absent**; the projection above is one explicit characteristic step, first order |
| Narrow-band taper of already-nodal values | none as an operator | **absent by design**; a taper invents nothing (it is the application's choice here) |

History (deleted in the redesign working tree, read via `git show HEAD:`): `src/Transfer/ExtensionOperator.jl`
was `ClosestPointExtension` with `extend_field(op, field::CartesianMeshField)` — 2D, Gridap-typed,
nodal→nodal, i.e. the same operation, not a facet→grid transfer. `src/Geometric/Curvature.jl` was the
2D predecessor of Curvature.jl and its role has *moved out* of EvolvingDomains, which is why the
cross-package call is natural today. `src/Kinematic/VelocitySource.jl` sampled velocities from
functions; it never extended interface values.

Curvature.jl's closest thing to a grid query is `curvatureAt(s, x...)` (`src/Query.jl:82`): exact
closest point on the facet set, then a *re-fit* there. Measured at 32³: 0.215 ms/node → 0.86 s for a
4000-node band, against 1.71 s for `interfaceCurvature` over 5808 facets. So a per-node query route is
about 2× cheaper than fitting every facet, but it can only return curvature (it refits internally); it
cannot transfer an arbitrary per-facet quantity. It is a viable *alternative* bridge for curvature
fields alone and was not used here.

**Conclusion:** the facet→node bridge had to be written in the application. It is 122 non-comment
lines and needs nothing from either package beyond what is already exported.

---

## 6. Stage C — the extension tested on a known scalar surface function

The bridge is exercised before any flow is run, against closed-form ground truth on the analytic
sphere, at three resolutions (`stageC/extension_test.csv`, figure `stageC/fig_extension_test.png`):

| n | facets | band nodes | unresolved | max \\|g−1\\| | max \\|g_x − x\*\\| | rms \\|g_x − x\*\\| | max \\|\\|v\\| − 1\\| | max \\|v − n_analytic\\| |
|---|---|---|---|---|---|---|---|---|
| 28 | 4200 | 3192 | 0 | 0.00e+00 | 2.72e−02 | 5.66e−03 | 2.2e−16 | 1.50e−01 |
| 32 | 5808 | 3968 | 0 | 0.00e+00 | 1.92e−02 | 4.65e−03 | 2.2e−16 | 1.26e−01 |
| 40 | 9000 | 6296 | 0 | 0.00e+00 | 1.84e−02 | 3.67e−03 | 2.2e−16 | 8.83e−02 |

- `g_f = 1` on every facet is reproduced **exactly** at every band node, and `v = n_f` is therefore a
  unit vector to round-off: the search, the index and the assembly are correct.
- `g_f =` facet-centroid x-coordinate must equal `R0 · p_x/|p|` at the node, the x-coordinate of the
  closest point of the *smooth* sphere. Measured error is first-order and decreasing
  (rms 5.66 → 4.65 → 3.67 × 10⁻³); it is the O(h) tangential offset between a piecewise-flat facet and
  the sphere, not an index error. Its max value saturates because the largest offsets sit at the band
  edge.
- `|v − n_analytic|` is an order of magnitude larger (0.15 → 0.088, first order) and is **not** an
  extension error: it is the O(h) error of the six-tetrahedron facet normals themselves (documented in
  `CutSurface`: differencing them gives zeroth-order curvature; using them for orientation only is the
  supported use). Its consequence for the flow is second order, because only the normal component of
  the velocity moves the interface.
- One transport-order self-check also runs at every start: a linear level set advected by a constant
  velocity reproduces the exact shift with error `0.000e+00`, which pins the "first axis fastest"
  velocity ordering that the assembly assumes.

---

## 7. Geometry A — sphere control

`R0 = 0.23`, `tmax = 0.12`, exact `R(0.12) = 0.2193172`. Radius is the area-weighted mean distance of
the cut-surface centroids; volume is the divergence-theorem volume of the same triangulation. Both are
zero-crossing estimates of the discrete interface, so they carry the piecewise-linear cut bias already
at `t = 0`.

| n | method | R(0) measured | R(0.12) | rel. error vs exact | shrink | exact shrink | shrink error | volume rel. error |
|---|---|---|---|---|---|---|---|---|
| 32 | Curvature.jl | 0.2291147 | 0.2184849 | −0.379 % | 0.010630 | 0.010683 | **−0.50 %** | −0.975 % |
| 32 | nodal FD | 0.2291147 | 0.2183690 | −0.432 % | 0.010746 | 0.010683 | +0.59 % | −1.133 % |
| 40 | Curvature.jl | 0.2294504 | 0.2187658 | −0.251 % | 0.0106846 | 0.0106829 | **+0.01 %** | −0.642 % |
| 40 | nodal FD | 0.2294504 | 0.2187123 | −0.276 % | 0.0107381 | 0.0106829 | +0.51 % | −0.715 % |

Reading these numbers honestly:

- The initial cut under-reports the radius by 0.386 % (32³) and 0.239 % (40³) **before any evolution**
  (measured). Roughly that much of every end error is geometry reconstruction, so the
  bias-cancelling quantity is the *displacement*: the Curvature.jl-driven front moves the exact
  distance to 0.5 % at 32³ and to 0.01 % at 40³, and it converges, while the FD operator stays at
  ≈ 0.5 % at both resolutions.
- Curvature against the exact value `2/R(t)`: 32³ `8.6610 ± 0.0172` at `t = 0` (exact 8.6957) and
  `9.0799 ± 0.3191` at the end (exact 9.1192, −0.43 %); 40³ `8.6824 ± 0.0079` → `9.1009 ± 0.0130`
  (−0.20 %). The facet-to-facet scatter is grid dependent (0.319 vs 0.013 for the same smooth
  geometry), which marks it as a lattice/fit-window effect, not physics: the *geometry* stays
  spherical meanwhile (sphericity 1.31 × 10⁻³ and 8.52 × 10⁻⁴ respectively). The FD band mean reads
  9.0107 (32³) and 9.0835 (40³).
- Surface area decreases monotonically in every run: 0.66140 → 0.60148 (−9.06 %, 32³) and
  0.66264 → 0.60248 (−9.08 %, 40³), as it must for a regular mean-curvature flow.
- Cross-sections stay circular to the printed precision of the crossings:
  `R_x, R_y, R_z = 0.218247, 0.218248, 0.218248` at 32³ and identical to 6 digits at 40³.
- `|∇φ|` in the band (measured): Curvature.jl run 0.9962 → 1.0023 (32³), 0.9975 → 1.0009 (40³); FD run
  0.9962 → 0.9857 and 0.9975 → 0.9938. Both are restored at each reinitialisation (§8). The FD velocity
  compresses the field more, which is a property of that operator, not of the transport.

---

## 8. Reinitialisation, measured with independent cuts

Every 10 steps the field is reinitialised as `reinitialize!(field, cut(field, :init))`, and the
represented geometry is measured from **three separately constructed** cut snapshots: cut #1 before,
cut #2 the one passed into `reinitialize!`, cut #3 after. Never the same snapshot twice.

| n | case | step | ΔR | ΔV/V | Δarea | facets |
|---|---|---|---|---|---|---|
| 32 | A cv | 10 | 0.000e+00 | 0.000e+00 | 0.000e+00 | 5352 → 5352 |
| 32 | A fd | 10 | 0.000e+00 | 0.000e+00 | 0.000e+00 | 5352 → 5352 |
| 32 | B cv | 10 | 0.000e+00 | 0.000e+00 | 0.000e+00 | 5384 → 5384 |
| 40 | A cv | 10, 20 | 0.000e+00 | 0.000e+00 | 0.000e+00 | 8808 → 8808, 8280 → 8280 |
| 40 | B cv | 10, 20 | 0.000e+00 | 0.000e+00 | 0.000e+00 | 8880 → 8880, 8480 → 8480 |

The operation clearly *does* change the field — the band-mean `|∇φ|` jumps from 0.9955 to 1.0025 at
32³ (A) and from 1.0080 to 1.0028 at 32³ (B) — while leaving the represented zero crossing bit-identical
at the printed precision. That is the intended contract ("seeding nodes of explicit CUT cells so the
represented zero crossings remain fixed"). The before/after pair is read from cuts #1 and #3, and the
cut handed to `reinitialize!` is #2, so no snapshot is reused as its own reference.

---

## 9. Geometry B — nonspherical perturbation, driven only by computed curvature

No analytic trajectory exists or is imposed: `radius_exact` is `NaN` in every B row. `Radius` below is
the area-weighted mean centroid radius, which for this shape is ≈ `R0` by construction of `P2`
(measured 0.231264 at 32³, 0.231595 at 40³ — the +0.6 % offset is the area weighting, not an error).

| quantity (32³, Curvature.jl) | t = 0 | t = 0.12 | change |
|---|---|---|---|
| mean radius | 0.231264 | 0.220414 | −4.69 % |
| volume | 5.11611e−02 | 4.43410e−02 | −13.33 % |
| surface area | 0.67328 | 0.61171 | −9.14 % |
| polar radius `R_z` | 0.263120 | 0.250204 | −4.91 % |
| equatorial radius | 0.211822 | 0.202916 | −4.20 % |
| anisotropy `R_z − R_eq` | 0.051298 | 0.047288 | **−7.82 %** |
| sphericity (std R / R) | 6.850e−02 | 6.581e−02 | −3.92 % |
| facets | 5784 | 5320 | −8.0 % |
| band-mean curvature ± std | 8.6525 ± 1.0750 | 9.0664 ± 1.0720 | +4.8 % |

At 40³ the same quantities move the same way (area −9.13 %, anisotropy 0.051468 → 0.047441, i.e.
−7.82 %, sphericity 6.849e−02 → 6.596e−02), so the qualitative behaviour is resolution independent.

Interpretation, with the required caution:

- For a smooth closed surface under mean-curvature flow the area must decrease while the solution is
  regular, and it does, monotonically, in both geometries and both resolutions. This is a qualitative
  check, not a reference solution.
- The P2 mode is smoothed faster than the mean radius shrinks: the anisotropy decays 7.8 % while the
  mean radius falls 4.7 %. That is the expected qualitative behaviour (higher harmonics of a sphere
  decay faster), and it is **not** a quantitative agreement with theory — no exact solution for this
  shape is claimed.
- Here the curvature std (≈ 1.08) is *physical*: the perturbed surface genuinely carries a spread of
  mean curvature, unlike the sphere, where the same estimator scatters by 0.017 (32³). It is not a
  discretisation artefact, and no closed-form curvature range is quoted because the P2 shape's
  principal curvatures have no simple `2/R` proxy.

---

## 10. Curvature.jl versus nodal finite differences

Same grid, same `dt` sequence, same WENO5 transport, same band, same reinitialisation schedule.

| n | geometry | R_end (Curvature.jl) | R_end (FD) | difference | κ CV (area mean) | κ FD (band mean) | rms CV−FD at facet centroids |
|---|---|---|---|---|---|---|---|
| 32 | A | 0.2184849 | 0.2183690 | 1.16e−04 | 9.0799 | 9.0107 | 0.370 |
| 32 | B | 0.2204140 | 0.2203280 | 8.60e−05 | 9.0664 | 8.9240 | 0.177 |
| 40 | A | 0.2187658 | 0.2187123 | 5.35e−05 | 9.1009 | 9.0835 | 0.078 |
| 40 | B | 0.2207513 | 0.2206858 | 6.55e−05 | 9.0596 | 9.0078 | 0.102 |

- The two driven trajectories never separate: the mean-radius difference is 0.03–0.05 % of `R` and
  *decreases* with refinement in A (1.16e−4 → 5.35e−5). Neither method is declared superior: on the
  sphere they bracket the exact shrink (−0.50 % / +0.59 % at 32³), which is what a resolution study of
  two consistent discretisations looks like.
- The two curvature *fields* differ more than the trajectories do (rms 0.37 at 32³, 0.078 at 40³,
  against a curvature scale of 9.1). This is the expected pattern: the interface position integrates
  the velocity, so a zero-mean curvature error largely cancels, while the pointwise difference does not.
- One inconsistency is left open, as measured: at 32³ the FD band mean curvature is 1.2 % *below* exact
  while its radius shrinks 0.6 % *faster* than exact. A band mean over `|φ| ≤ 0.8h` and the curvature
  actually seen by the zero crossing are not the same estimator, so these two numbers need not agree to
  the percent level. **open** — resolving it needs the curvature sampled *on* the interface for both
  operators, which the FD field only provides by interpolation.
- The histograms in `fig_curvature_distribution.png` show the character of the two operators. At 40³ the
  Curvature.jl facet curvatures of the sphere form a single narrow peak (std 0.0079) while the nodal FD
  values are visibly multi-modal over ≈ 8.55–8.80: differences of a piecewise-linear normal field
  inherit the lattice, a local polynomial fit of the facet geometry does not. At 32³ that ordering
  *reverses* — the Curvature.jl scatter is 0.319 against a similar FD spread — which is the same
  unexplained lattice/fit-window sensitivity recorded under §7 and §13, not a property of either
  operator that can be stated once and for all.

---

## 11. Performance (measured, warm, single-threaded, `JULIA_NUM_THREADS=1`)

Per step, milliseconds. Timings are accumulated around each operator after a warm-up cut, curvature call
and zero-length transport step, so first-call compilation is not charged to the loop. `cut` includes the
two extra cuts taken on reinitialisation steps and the per-snapshot cut.

| n | case | cut | Curvature.jl | extension | FD curvature | WENO5 | reinit | total |
|---|---|---|---|---|---|---|---|---|
| 32 | A cv | 54.7 | **1622.8** (95.2 %) | **17.2** (1.01 %) | 0.10 | 5.0 | 3.9 | 1703.7 |
| 32 | B cv | 53.1 | **1784.0** (95.7 %) | **18.9** (1.01 %) | 0.12 | 7.1 | 1.7 | 1864.9 |
| 40 | A cv | 71.1 | **2069.1** (95.1 %) | **22.8** (1.05 %) | 0.12 | 8.4 | 4.8 | 2176.3 |
| 40 | B cv | 83.7 | **2443.5** (95.1 %) | **27.6** (1.07 %) | 0.12 | 10.9 | 4.7 | 2570.4 |
| 32 | A fd | 57.0 | (522.7, snapshots only) | — | 0.43 | 7.5 | 2.2 | 589.9 |
| 40 | B fd | 78.4 | (493.6, snapshots only) | — | 0.61 | 11.8 | 6.5 | 591.0 |

- **Curvature.jl / extension ratio: 89–95×.** The bridge is not the bottleneck and is not close to it.
- The extension costs 4.3 µs per band node at 32³ and 3.6 µs at 40³ (3968 and 6296 band nodes).
- In the `fd` rows the `curvature` column is *not* the driver — it is the five snapshot evaluations of
  Curvature.jl kept for comparison (each ≈ 1.5–2.5 s). The genuine FD loop cost is
  cut + fd + transport + reinit ≈ 67 ms/step at 32³, so the Curvature.jl-driven step is ≈ 25× more
  expensive end to end.
- Allocation (probe, 32³): one `interfaceCurvature` call allocates ≈ 82 MB, which is why the in-loop
  per-step cost (1.6 s) is 60 % above the best-of-three probe value (1.0 s) — GC pressure, not a
  different algorithm.
- Fit-window sensitivity (probe, sphere at `t = 0`): `radius = 3h` gives 8.688, `4h` gives 8.661, `6h`
  gives 8.495 against an exact 8.696, at 0.59 s / 1.03 s / 2.97 s per call. `4h` was kept for
  comparability with round 1; the probe shows `3h` is both more accurate and ~1.7× cheaper, so a
  follow-up should use it.

---

## 12. Figures (CairoMakie only; produced by `visualize.CurvatureExtension3D.jl` from saved outputs)

`output/n32/` and `output/n40/`, same seven figures each:

| file | content |
|---|---|
| `fig_interfaces_3d.png` | 3D interface at t = 0, 0.25, 0.5 `tmax` for A and B, facets coloured by the Curvature.jl mean curvature (polar high / equatorial low on B, as it must be) |
| `fig_cross_sections.png` | orthogonal cross-sections (`z = 0`, `y = 0`, `x = 0`) of the final level set with the zero contour: circular for A, prolate for B |
| `fig_curvature_distribution.png` | per-snapshot histograms of the Curvature.jl and nodal-FD curvature on the interface, plus the pointwise CV-vs-FD scatter at the same facet centroids |
| `fig_radius_vs_exact.png` | A: measured radius against `sqrt(R0² − 4μt)` and the absolute radius error for both methods |
| `fig_anisotropy.png` | B: polar and equatorial crossings, and the anisotropy `R_z − R_eq` with sphericity, per method |
| `fig_cv_vs_fd.png` | radius difference between the two drivers, and the surface-area history of both |
| `fig_cost.png` | stacked per-step cost breakdown for all four runs |
| `output/stageC/fig_extension_test.png` | Stage C extension error versus `h` (log-log) |

---

## 13. Limits — what this does not establish

- **The facet velocity is piecewise constant.** The extension is a closest-point transfer of facet data,
  so the velocity field is discontinuous across facet boundaries and only first-order accurate as a
  field. It is adequate here because the interface integrates the velocity; it would not be adequate
  for a problem that differentiates the velocity field.
- **The taper is not a normal extension.** Beyond `3h` the velocity is zero, not extended. This is
  documented in code and is safe while the band exceeds the per-step displacement (0.023 h).
- **`CutSurface` normals are first order**, so the velocity direction carries an O(h) error
  (measured 0.15 → 0.088 in Stage C). Its effect on the interface is second order, but no claim of
  second-order interface accuracy is made; the measured interface errors are 0.5 % (32³) and 0.01 %
  (40³) in displacement, which is consistent with first-order-plus-cancellation, not with a proven rate.
- **Two resolutions only.** A convergence study needs at least three; 24³/28³/32³/40³ are cheap
  enough (curvature cost permitting) and are the obvious next step.
- **Time step chosen for stability, not accuracy.** `C = 1/12` over-damps the temporal error relative to
  the spatial one (0.023 h per step); the `h²` scaling is a diffusion-scale heuristic, and no rigorous
  bound is claimed for the assembled algorithm.
- **The 32³ curvature scatter (0.32 vs 0.013 at 40³) is unexplained beyond "lattice/fit-window effect"**
  — **open**. It does not corrupt the mean (0.43 % / 0.20 %) but it deserves the 3h window test.
- Curvature.jl is used as a black box with the documented public call only; no upstream change, patch
  or private function was involved, and none was needed.

---

## 14. Reproduction

```bash
cd "/Users/naoufel/Documents/Scientific stack/EvolvingDomains-Redesign"
JULIA_NUM_THREADS=1 julia --project=validation/env \
    validation/MeanCurvatureSphere3D/round2/run.CurvatureExtension3D.jl 32 0.12   # primary, ~1 min
JULIA_NUM_THREADS=1 julia --project=validation/env \
    validation/MeanCurvatureSphere3D/round2/run.CurvatureExtension3D.jl 40 0.12   # refinement, ~2 min
JULIA_NUM_THREADS=1 julia --project=validation/env \
    validation/MeanCurvatureSphere3D/round2/visualize.CurvatureExtension3D.jl \
    validation/MeanCurvatureSphere3D/round2/output
```

The run script locates the native cut library and the Curvature.jl checkout itself
(`EVOLVINGDOMAINS_CUTCELLS_LIBRARY` / `CURVATURE_JL_DIR` override them), runs the transport-ordering
check and the Stage C extension test first, then the four `(geometry, method)` runs, then writes
`comparison.csv`. Set `CURVATUREEXT_SKIP_FIGURES=1` to run the simulation without figures.
Reruns are bit-reproducible (verified: identical radii to 10 digits before and after a column-only
change to the CSV writer).
