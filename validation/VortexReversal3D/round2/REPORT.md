# StrongDeformation3D — Campaign II, Experiment 2

Three-dimensional reversible transport on `[0,1]³` with an off-centre ellipsoid and a
Gaussian marker, transported by WENO5 (level set and scalar) and CIP (scalar), plus an
isolated steady-state allocation study of every operation of the time loop.

Artifacts: `run.StrongDeformation3D.jl`, `visualize.StrongDeformation3D.jl`, `output/`
(and `output_cfl02/`, the CFL-halving refinement). EvolvingDomains was used unmodified.

---

## 1. Conclusions

**The geometry now deforms strongly, and the deformation is measured, not asserted.**
At `t = 0.5` the symmetric difference between the initial and the midpoint domains is
`D_sym = 1.40` (smeared indicator) / `1.41` (voxel indicator) on a scale where two disjoint
sets of equal volume give `2.0`; the Jaccard index is `0.13`, i.e. only 13 % of the union is
shared. The interface's reconstructed vertex count grows by 44 % at mid-interval
(16 872 → 24 216 at 64³), so the surface is genuinely stretched, not merely translated. At
`t = 0.5` all three orthogonal slices show a non-elliptic, lobe-and-neck configuration
(`output/slices_n64.png`), and the five-time montage
(`output/montage_interface_n64.png`) shows the ellipsoid returning to a compact shape at
`t = 1`. `A = 2` was retained: no run produced a non-finite value or a CIP
characteristic-crossing error.

**The return error is dominated by numerical diffusion of an over-stretched interface, not
by temporal error.** Halving `Δt` at fixed 48³ (CFL 0.4 → 0.2, 472 → 940 steps) changes the
level-set return error in the fourth significant digit (`5.0479e-3 → 5.0478e-3`) and leaves
the WENO5 scalar error and the volume drift unchanged; only the CIP scalar moves, by 7.1 %
in L1. Spatial resolution is what the error responds to: L1 volume drift at `t = 1` falls
from −30.3 % (48³) to −21.9 % (64³).

**The entire steady-state per-step allocation of the loop is CIP.** After compilation and
cache initialisation, at 64³:

| operation | bytes | bytes/node/step | s/step |
|---|---:|---:|---:|
| velocity update (`dynamics`) | 0 | 0 | 0.0161 |
| one WENO5 geometry step | 2 832 | 0.011 | 0.0474 |
| one WENO5 scalar step | 2 832 | 0.011 | 0.0537 |
| one CIP scalar step (loop path) | 270 552 112 | 1 032.1 | 0.1831 |
| CIP source/target copy-back | 0 | 0 | 3.1e-5 |
| snapshot `copy` | 0 | 0 | 4.2e-8 |
| one `cut` construction (diagnostic) | 431 386 960 | 1 645.6 | 0.1139 |
| whole loop, measured | 171 032 362 240 | **1 032.3** | 0.2672 |

`1 032 B/node/step` is *exactly* the CIP row, and it is resolution independent
(48³: 1 032.5 B/node/step). Velocity construction, source/destination buffer management and
snapshots allocate **nothing**: the application-level share of the allocation is
`2 × 2 832 B` per step, i.e. 0.005 %.

**Inside CIP the bytes are fully accounted for**: `1 032 B/node = 129` full-grid `Float64`
arrays per step, of which 108 are `_derivative` outputs (36 per directional sweep) and 21
are the scaffolding of the velocity-partials table. A single `CIPCache` in-loop call
performs 644 heap allocations in total (64³), i.e. ~128 grid-sized arrays plus rounding.
`CIP` is *not* "allocation-heavy because the simulation allocates a lot": it is
allocation-heavy because it materialises 36 derivative fields per sweep instead of the 19 it
uses, and because `reduce(+, vals)` produces five temporaries.

**Adapter costs are negligible in both bytes and time.** The non-aliasing CIP interface costs
one `copyto!` per step: 0 bytes, 3.1e-5 s (0.017 % of the CIP step). The one adapter cost a
caller can *create* is the interpolation-profile rebuild: if the buffer passed as `source` is
not the previous `target`, CIP re-derives the profile from nodal values, which costs 7 extra
grid-sized arrays = `56 B/node/step` (+5.4 % over the loop path), measured as
`cip_scalar` vs `cip_scalar_loop`.

**`cut` is not cheap**: 1 646 B/node per construction at 64³ (431 MB), 1.6 × a CIP step's
bytes and 0.6 × its wall time, called 5 times here. It is diagnostic-only and outside the
timed loop, but a per-step interface extraction would dominate everything else.

---

## 2. Environment and source verification (Stage A)

Printed by `verifySource()`; the same text is saved as `output/verification.txt`.

```text
EvolvingDomains source : /Users/naoufel/Documents/Scientific stack/EvolvingDomains-Redesign/src/EvolvingDomains.jl
  package dir          : /Users/naoufel/Documents/Scientific stack/EvolvingDomains-Redesign
  version (Project.toml): 1.0.0
  git revision         : 0809db8bf38e5d25975022a1848798f2219289c8
                         (60 modified/untracked paths, pre-existing in the checkout)
  CutCells native lib  : .../validation/lib/cutcells/lib/libevolvingdomains_cut_cells.dylib
  julia 1.11.7, 1 thread(s), StaticArrays 1.9.24
```

`pathof(EvolvingDomains)` and `pkgdir(EvolvingDomains)` were compared with
`realpath(@__DIR__/../../..)` and the run aborts if they differ. `CairoMakie 0.15.17` is used
by the visualization only. Public API used, nothing else:

`Grid`, `Discretization`, `Field`, `eulerianTransport!`, `lagrangianTransport!`,
`WENO5Cache`, `CIPCache`, `cut`, `Cut`, `CUT_OUT`, plus `field.data[:init]` for the plain
arrays. No package file was modified; the protocol requires no recompilation of the package.

Two environment facts that shaped the measurements: `Threads.nthreads() == 1` (so the
`@threads` in `weno5Rhs!` takes the single-thread branch, and its 1 872 B of allocation per
step is `@threads` bookkeeping, not numerical work), and the CutCells library is not
resolvable from an artifact in this checkout, so the documented
`EVOLVINGDOMAINS_CUTCELLS_LIBRARY` override points at the prebuilt library in
`validation/lib/`. The fallback path (report `missing` volume) exists but was never taken.

---

## 3. Problem, initial data, and why the midpoint now deforms

### 3.1 Flow

```text
v(x,t) = A cos(πt) V(x),   A = 2,   0 ≤ t ≤ 1

        ⎛  2 sin²(πx) sin(2πy) sin(2πz) ⎞
V(x) =  ⎜ −sin(2πx) sin²(πy) sin(2πz)  ⎟
        ⎝ −sin(2πx) sin(2πy) sin²(πz)  ⎠
```

Divergence-free and tangent to `∂D` by construction (each divergence term carries a
`sin(2πx_i)` factor and cancels against the other two). Because
`cos(π(1−t)) = −cos(πt)`, the flow map over `[0,1]` is the identity, hence
`u(·,1) = u(·,0)` exactly for both the level set and the scalar. No manufactured solution is
needed: every error reported below is a *return* error.

**Discrete check** (`verifyFlow`, sampled on the 64³ nodes):

```text
max |∇·V| = 3.55e-15      max |V·n| on ∂D = 0.00e+00
max |V|   = 1.997841      vmax = A·max|V| = 3.995682
```

### 3.2 Initial data

```text
φ₀(x) = a_min ( √( ((x−c_x)/a_x)² + ((y−c_y)/a_y)² + ((z−c_z)/a_z)² ) − 1 )
c  = (0.41, 0.48, 0.56),  (a_x,a_y,a_z) = (0.14, 0.11, 0.09),  a_min = 0.09
c₀(x) = exp( −|x − (0.42,0.48,0.56)|² / 0.10² )
V₀    = 4π a_x a_y a_z / 3 = 5.8056632e−3
```

The zero level set is *exactly* the ellipsoid, negative inside; `φ₀` is **not** a signed
distance function. Everywhere in the box `|∇φ₀| ∈ [a_min/a_max, a_min/a_min] = [0.643, 1.0]`
(along the `i`-th semiaxis it equals `a_min/a_i`), so `|φ₀(p)|` *underestimates* the Euclidean
distance to the interface by up to a factor `1/0.643 = 1.56`; the function also has a conical
non-differentiable point at the ellipsoid centre, harmless here because it is deep inside the
negative phase. Every distance below is
reported normalised by `a_min = 0.09` and inherits this one-sided bias; this is stated again
where the numbers are used.

The ellipsoid contains the box centre `(0.5,0.5,0.5)`
(`Σ((0.5−c_i)/a_i)² = 0.89 < 1`) and lies strictly inside the box
(`x ∈ [0.27,0.55]`, `y ∈ [0.37,0.59]`, `z ∈ [0.47,0.65]`), so `wall_contacts = 0` for every
snapshot and the cut-cell volume is a legitimate enclosed volume.

### 3.3 Why the previous geometry was easy and this one is not

The previous campaign observed that its centred sphere sat almost on an invariant surface of
the old velocity. That structure is **still present** in the new velocity. With

```text
I(x) = sin²(πx) sin²(πy) sin²(πz),   sⱼ = sin(πxⱼ),  cⱼ = cos(πxⱼ)
V·∇I = 2π ( 16 − 8 − 8 ) s_x³c_x s_y³c_y s_z³c_z = 0
```

using `sin(2πx) = 2 s_x c_x` in each of the three terms. `I` is therefore a first integral of
the flow: every trajectory stays on its level set of `I`, and those level sets are
approximately spheres about the box centre near `I ≈ 1`.

What changed is the *geometry*, not the flow. The ellipsoid is displaced by
`|c − (½,½,½)| = 0.112` from the centre of the invariant shells and its radial extent in that
frame spans roughly `−0.03 … 0.20`, so its surface cuts across the whole shell family instead
of hugging one shell. Transport along the shells then shears and stretches the surface: the
measured `D_sym(0.5) = 1.40` and the 44 % area growth are the consequence. This is a
*stronger* test than the centred sphere (large tangential displacement, large shear,
thin filaments), but it is honest to state its nature: the exact flow still preserves `I`, so
the deformation is tangential shear of a non-invariant surface, not the imaging of a surface
into a topologically different object. The numerical difficulty it creates — inverting that
shear — is real: it is what the return errors and the volume drift measure.

---

## 4. Numerical methods and the executed loop

Time step from a CFL bound on the *nodal* velocity maximum at `t = 0` (the worst instant,
`cos = 1`), rounded down so that `steps` is a multiple of 4 and `t = 0.25, 0.5, 0.75` are
step boundaries:

```text
Δt = CFL · h / vmax,   CFL = 0.4,   N = 4·⌈1/(4Δt)⌉,   Δt ← 1/N
```

| nodes | h | Δt | steps | loop wall time | GC | allocated |
|---|---:|---:|---:|---:|---:|---:|
| 48³ | 2.128e-2 | 2.11864e-3 | 472 | 54.2 s | 1.6 s | 53.9 GB |
| 64³ | 1.587e-2 | 1.58228e-3 | 632 | 168.9 s | 16.7 s | 171 GB |

Both intervals ran to `t = 1`; all saved values are finite.

The loop is exactly the prescribed one:

```julia
for step in 1:steps
    v = dynamics(config, (step - 0.5) * dt)     # midpoint velocity, one array, reused
    eulerianTransport!(geometry, v, dt; level_id=:init, cache=weno)    # φ, WENO5+SSP-RK3
    eulerianTransport!(scalarEuler, v, dt; level_id=:init, cache=weno) # c, WENO5+SSP-RK3
    lagrangianTransport!(scalarCipTarget, scalarCip, v, dt; level_id=:init, cache=cip)
    copyto!(scalarCip.data[:init], scalarCipTarget.data[:init])        # CIP is non-aliasing
end
```

The same physical velocity array is passed to all three transports. Midpoint freezing is
second order in `Δt` and *exactly* time-symmetric (`cos(π(1−t)) = −cos(πt)`), so it does not
by itself break reversibility. **No reinitialization is used**: `reinitialize!` reconstructs a
distance function, which would move the zero level set and destroy the quantity under test.
Its absence is not free — the level set loses amplitude (see §7) — but a reinitialized run
would not be a reversibility measurement.

### 4.1 Velocity node ordering (Stage B)

The ordering convention ("first axis fastest", one `SVector` per node in `vec` order) is
untyped and load-bearing. The test uses an affine field and a velocity that depends on `y`
only, `v = (0.01, 0, 0.006)(0.5 + y)`, for which the exact solution is known in closed form,
`φ(x,t) = φ₀(x) − t(0.5+y)(2·0.01 + 3·0.006)`, because `v_y = 0` makes `y` a trajectory
invariant. Both kernels transport it exactly (WENO5 reconstructs a linear function exactly,
SSP-RK3 integrates the per-node linear ODE exactly, every CIP Hermite profile is a straight
line). The same test is then repeated with the velocity array attached to a transposed node
index, which must fail:

| resolution | WENO5 (correct) | CIP (correct) | WENO5 (transposed) | CIP (transposed) |
|---|---:|---:|---:|---:|
| 48³ | 4.72e-9 | 1.78e-15 | 5.99e-3 | 5.98e-3 |
| 64³ | 3.19e-9 | 1.78e-15 | 6.40e-3 | 6.39e-3 |

A transposed ordering is six orders of magnitude worse in both kernels, so the test is
sensitive to what it claims to detect. The CIP residual is at round-off; the WENO5 residual of
a few `1e-9` is not a transport error but a boundary effect of the clamped stencil next to the
zero-velocity shell (widening the comparison window to 7 nodes on every axis drives it to
`1.8e-15`, measured), which is why the assertion threshold is `1e-7` rather than machine
epsilon.

---

## 5. Geometric validation

All numbers in `output/n{48,64}/metrics.csv`, `interface.csv`, and printed in `metrics.txt`.
Three independent volume estimators are carried: the nodal voxel count (`count(φ<0)h³`), a
one-cell smeared Heaviside (`Σ clamp(½ − φ/h)`), and a divergence-theorem integral
`V = ⅓∮ x·n dA` over the facets of the package's own `cut` snapshot.

### 5.1 Volume and drift

Exact initial volume `4πabc/3 = 5.8056632e-3`.

| t | V_voxel | V_smeared | V_cutcell | drift voxel | drift smeared | drift cut |
|---:|---:|---:|---:|---:|---:|---:|
| 0.00 | 5.77891e-3 | 5.85207e-3 | 5.74322e-3 | — | — | — |
| 0.25 | 5.23502e-3 | 5.35913e-3 | 5.17491e-3 | −9.41 % | −8.42 % | −9.90 % |
| 0.50 | 4.85509e-3 | 4.88706e-3 | 4.62279e-3 | −15.99 % | −16.49 % | −19.51 % |
| 0.75 | 4.62713e-3 | 4.70048e-3 | 4.54914e-3 | −19.93 % | −19.68 % | −20.79 % |
| 1.00 | 4.51515e-3 | 4.56660e-3 | 4.42675e-3 | −21.87 % | −21.97 % | −22.92 % |

(64³; 48³ gives −30.3 % / −30.5 % / −32.4 % at `t = 1`.) The three estimators agree to
0.2–1.0 percentage points on the *drift*, which is what matters; their common bias
(the initial `cut` volume is 1.1 % low, the smeared one 0.8 % high) largely cancels in the
ratio. `wall_contacts = 0` at every snapshot, so no boundary facet contaminates the
cut-cell integral.

The exact flow preserves volume, so the whole of this drift is numerical: the interface is
thinned below the grid scale by the shear and its zero contour retreats. The drift roughly
halves between 48³ and 64³ (−30.3 % → −21.9 %, ratio 1.39 for an `h` ratio of 1.34), i.e.
approximately first order in `h` over this range; with two resolutions no order is claimed.

### 5.2 Midpoint deformation

`D_sym = |Ω(0) △ Ω(t)| / |Ω(0)|`, both sets discretised on the same nodal grid (smeared
indicator; the voxel (sign-count) version is the cross-check). `Jaccard = |A∩B|/|A∪B|` and
`overlap = |A∩B|/|Ω(0)|` are added because a *shrinking* `Ω` inflates `D_sym` on its own: the
Jaccard index is the volume-loss-independent part of the statement.

| t | D_sym smeared | D_sym voxel | Jaccard | overlap | interface vertices |
|---:|---:|---:|---:|---:|---:|
| 0.00 | 0 | 0 | 0.838 | 0.912 | 16 872 |
| 0.25 | 1.291 | 1.282 | 0.192 | 0.308 | 21 216 |
| 0.50 | **1.401** | **1.407** | **0.132** | **0.214** | 24 216 |
| 0.75 | 1.266 | 1.260 | 0.172 | 0.265 | 18 168 |
| 1.00 | 0.223 | 0.219 | 0.685 | 0.724 | 16 908 |

(64³.) `D_sym ≤ 2` with equality for disjoint sets of equal volume; at `t = 0.5` the two
configurations share 13 % of their union, and even the overlap of the *initial* set alone
(0.214) leaves no room for the shape change to be a volume-loss artefact. The midpoint
configuration is therefore demonstrably different from the initial ellipsoid.

The vertex count of the reconstructed interface is an area proxy (facets per unit area
scales with area): 16 872 → 24 216 (+44 %) at `t = 0.5`, back to 16 908 (+0.2 % of the
initial count) at `t = 1`. So the surface is stretched by about 44 % at mid-interval — a
deformation measure that is independent of the set's displacement.

`t = 1` is not a deformation statement: the exact `D_sym(1)` is 0, so the residual 0.223 is
pure numerical error of the return.

### 5.3 Interface return distance

Interface points are the reconstructed vertices of `cut` (`boundary.physical_points`,
16 908 at 64³), never an analytic surface. Two directions are evaluated, so the measure is
two-sided:

* **forward** — every reconstructed vertex `p` of `φ(t)` evaluated in the *exact* algebraic
  `φ₀`: `|φ₀(p)|`. No interpolation error in the reference.
* **backward** — 4 000 deterministic Fibonacci points of the exact initial ellipsoid, where
  `φ(t)` is sampled by trilinear interpolation: `|φ(t)(p)|`.

Values are `distance / a_min`; multiply by 0.09 for the absolute error.

| t | mean_fwd | max_fwd | mean_bwd | max_bwd |
|---:|---:|---:|---:|---:|
| 0.00 | 1.69e-3 | 8.00e-3 | 3.60e-3 | 6.40e-3 |
| 0.25 | 0.5673 | 1.5874 | 0.6012 | 1.4473 |
| 0.50 | 0.5758 | 1.3401 | 0.7073 | 1.5061 |
| 0.75 | 0.5739 | 1.4740 | 0.6244 | 1.2102 |
| 1.00 | **0.0938** | **0.7271** | 0.1157 | 1.2019 |

(64³; 48³ `t = 1`: 0.1175 / 0.7791 / 0.1619 / 1.2560.) At `t = 1` the forward error is a mean
of `8.4e-3` (`0.53 h`) and a maximum of `6.5e-2` (`4.1 h`).

**What these numbers are and are not.** They are a *sampled, two-sided, level-set-based
surrogate*: sampled because only the extracted vertices are tested (no maximisation over the
surface), level-set-based because `φ` is neither a signed-distance function nor, after
300–600 steps of transport, a good approximation of one. The `t = 0` row is the diagnostic
floor (`1.7e-3` mean, `8.0e-3` max) and shows the extraction and the interpolation agree with
the exact geometry. Two systematic biases are quantified rather than hidden:

* the forward direction underestimates distance by up to 1.56× because `|∇φ₀| ≥ 0.643`;
* the backward direction is inflated relative to the forward one (0.707 vs 0.576 at
  `t = 0.5`) because `|∇φ|` decays under diffusion — `min φ` falls from `−0.0843` to
  `−0.0618`, a 27 % amplitude loss — so the same geometric displacement reads as a larger
  `|φ|` when the reference is the numerical field rather than the exact one.

A true Hausdorff distance would need both surfaces and a maximisation over a dense covering of
each; that is not what is reported, and the numbers above should not be quoted as one.

---

## 6. Return errors and resolution study

Relative return errors at `t = 1` over the nodal values (ε_L1 = Σ|u_N−u₀|/Σ|u₀|,
ε_L2 = ‖u_N−u₀‖/‖u₀‖):

| quantity | 48³ L1 | 48³ L2 | 64³ L1 | 64³ L2 | ratio 48→64 (L1) |
|---|---:|---:|---:|---:|---:|
| φ (WENO5) | 5.0479e-3 | 1.9269e-2 | 3.1406e-3 | 1.5284e-2 | 1.61 |
| c (WENO5) | 2.1562e-1 | 2.9022e-1 | 1.4399e-1 | 2.2991e-1 | 1.50 |
| c (CIP) | 1.1255e-1 | 1.4129e-1 | 7.9745e-2 | 1.0768e-1 | 1.41 |

CIP is 1.5–1.8 × more accurate than WENO5 for the same scalar and the same velocity, at both
resolutions, and both scalar errors are far larger than the level-set error — the smooth
Gaussian is a harder test than the interface because it is stretched into sub-cell filaments
(visible as the detached lobe in `output/scalar_schemes_n64.png` at `t = 0.5`).

With two resolutions and `h` ratio 1.34 the L1 ratios 1.41–1.61 correspond to apparent orders
1.2–1.6, and the L2 ratios 1.26–1.31 to apparent orders 0.8–0.9. **No asymptotic order is
claimed.** The errors are dominated by under-resolution of thin structures, which is exactly
where a formal order does not apply.

### 6.1 Temporal versus spatial error (CFL halving at 48³)

`output_cfl02/` is the same 48³ experiment at `CFL = 0.2` (940 steps, `Δt` halved, everything
else identical; produced with `SD3D_CFL=0.2 SD3D_OUT=_cfl02`):

| quantity | CFL 0.4 | CFL 0.2 | change |
|---|---:|---:|---:|
| φ return L1 | 5.0479e-3 | 5.0478e-3 | < 0.01 % |
| φ return L2 | 1.9269e-2 | 1.9269e-2 | < 0.01 % |
| c WENO5 L1 | 2.1562e-1 | 2.1562e-1 | 0 % |
| c CIP L1 | 1.1255e-1 | 1.0459e-1 | −7.1 % |
| c CIP L2 | 1.4129e-1 | 1.4045e-1 | −0.6 % |
| volume drift at t=1 | −30.3 % | −30.3 % | 0 % |
| D_sym(0.5) | 1.343 | 1.343 | 0 % |

The WENO5/SSP-RK3 paths are time-step independent at this CFL: their error is spatial. Only
CIP responds, by 7 % in L1, consistent with the first-order directional splitting in time.
This is the clean separation the first campaign could not make, and it justifies spending the
resolution budget on `h` rather than on `Δt`.

---

## 7. Allocation investigation

Everything below is measured after compilation and cache initialisation, on dedicated work
fields, as the **minimum `@allocated` over 6 calls** (the allocation counter is deterministic
here: repeated runs give identical bytes) and the **median `@timed`** of 5 calls.
`BenchmarkTools` is not in the pinned `validation/env` manifest and was deliberately not
added (no environment changes); repeated `@allocated`/`@timed` gives exact byte counts, and
the accepted noise is only in the seconds column (±10–20 % run to run).
Raw data: `output/n{48,64}/allocs.csv`, `allocs_profile.csv`.

### 7.1 Per-operation cost

| operation | 48³ bytes | 48³ B/node | 64³ bytes | 64³ B/node | 64³ s/step |
|---|---:|---:|---:|---:|---:|
| velocity update | 0 | 0 | 0 | 0 | 0.0161 |
| WENO5 geometry step | 2 832 | 0.026 | 2 832 | 0.011 | 0.0474 |
| WENO5 scalar step | 2 832 | 0.026 | 2 832 | 0.011 | 0.0537 |
| CIP step, loop path | 114 150 448 | 1 032.2 | 270 552 112 | 1 032.1 | 0.1831 |
| CIP step, profile rebuilt | 120 344 272 | 1 088.2 | 285 232 848 | 1 088.1 | 0.1905 |
| source→target copy-back | 0 | 0 | 0 | 0 | 3.1e-5 |
| snapshot `copy` (2 MB array) | 0 | 0 | 0 | 0 | 4.2e-8 |
| `cut` construction | 208 409 056 | 1 884.5 | 431 386 960 | 1 645.6 | 0.1139 |
| whole loop body | 114 156 112 | 1 032.2 | 270 557 776 | 1 032.1 | 0.2917 |
| loop as measured (`@timed`, 632 steps) | 53 898 496 768 | 1 032.5 | 171 032 362 240 | 1 032.3 | 0.2672 |

The isolated loop-body measurement (1 032.1 B/node) and the `@timed` total over 632 steps
(1 032.3 B/node) agree to 0.02 %, at both resolutions. **Bytes per step** are `114.2 MB`
(48³) and `270.6 MB` (64³); **bytes per node per step** are `1 032` in both cases; **time per
step** is `0.115 s` and `0.267 s`.

Shares of the per-step allocation at 64³: CIP `270.55 MB` (**99.995 %**), the two WENO5 calls
`5 664 B` (0.005 %), velocity/copy-back/snapshot `0 B`. Shares of the per-step wall time
(from the medians): CIP 61 %, the two WENO5 steps 34 %, velocity update 5 % (the sum of the
medians exceeds the loop mean by 12 %, within the measurement noise).

**Memory**: peak RSS (`Sys.maxrss`, read immediately after the time loop of the 64³ run, in a
process that had done nothing heavier before) is **1.51 GB**, against a live set of ≲100 MB
(4 fields 8.4 MB, `CIPCache` 16 arrays 33.6 MB, velocity 6.3 MB, 5 snapshots × 3 fields
31.5 MB). The 15× ratio is the footprint of the 270 MB/step churn, i.e. of garbage that GC
has not yet returned. GC itself costs 4.6 s (2.7 % of the loop) in the most favourable run
and 16.7 s (9.9 %) in the least; large arrays are cheap to free, so the churn is visible in
bytes and RSS but not catastrophic in time.

### 7.2 Where the CIP bytes come from (`Profile.Allocs`, sample_rate = 1)

One in-loop CIP step at 64³: **644 heap allocations, 270 550 312 B**. (This is the honest way
to answer "is CIP itself allocation-heavy": the profile is attached to a single public call,
and the package was not modified or instrumented.)

| site (innermost package frame) | bytes | allocations | B/node | grid-sized arrays |
|---|---:|---:|---:|---:|
| `_derivative` (CIP.jl) | 226 500 480 | 405 | 864.0 | **108** |
| `_velocityPartials3` | 44 049 000 | 180 | 168.0 | **21** |
| `_sweep3!` | 408 | 12 | 0.0 | 0 |
| `_stride` | 360 | 45 | 0.0 | 0 |
| `fieldCipLagrangianTransport!` closure | 64 | 2 | 0.0 | 0 |

`226 500 480 / (8·64³) = 108.0` and `44 049 000 / (8·64³) = 21.0`: the bytes are exactly
whole grid-sized `Float64` arrays, so the accounting closes.

Per step: 3 directional sweeps × 43 arrays = 129 arrays = **1 032 B/node**. Per sweep:

| contribution | arrays | note |
|---|---:|---|
| ∂_i v_a (first order) | 3 | needed |
| ∂²v_a (second order) | 6 | needed |
| ∂³v_a, non-mixed third order | 9 | needed |
| ∂³v_a/∂x∂y∂z via 6 permutation orders | 18 | 6× redundant, then averaged |
| `values = [v[a] for v in velocity]` | 1 | full copy of one velocity component |
| `reduce(+, vals)` for the 6 permutations | 5 | pairwise addition, no accumulator |
| `./ length(vals)` | 1 | |

So 55.8 % of the per-sweep temporaries (24 of 43) belong to the symmetrised mixed third
derivative, and the derivative *reconstruction* of the transported profile (`φ`, `∂φ`, and
the 15 mixed derivatives, 7 arrays) is not the problem — it is skipped entirely on the loop
path because `CIPCache` correctly recognises that its snapshot still describes the source.
That is the 5.4 % difference between the `cip_scalar_loop` and `cip_scalar` rows: the cache
contract works, and it is worth `56 B/node/step` only when the caller's loop leaves the source
equal to the previous target, which the natural copy-back loop does.

No conclusion in this report is drawn from the *total* allocation of the simulation alone:
each operation was measured in isolation, and inside CIP the attribution is by allocation
site.

### 7.3 `cut` and the adapter surface

`cut` allocates 431 MB (1 646 B/node) and 0.114 s per construction at 64³ —
`Profile.Allocs` attributes it to 1 639 338 allocations, dominated by `_append3DPoints!`
(225 MB in 80 allocations) and `_append3DFullCell!` (166 MB in 1.49 M small allocations, the
per-cell interface polygons). It is *not* in the timed loop: five calls for the snapshots, plus the
handful of calls made by the allocation study and the profile. A workflow that extracted the interface
every step would pay more for the diagnostic than for the transport.

Application-level (adapter) costs, measured separately: velocity array fill 0 B; CIP
copy-back 0 B and 3.1e-5 s per step (0.017 % of the CIP step); snapshot copies 0 B (15 array
copies = 31.5 MB over the whole 64³ run = 2e-5 % of the loop's 171 GB); `Field` construction
~2 kB. Nothing in the application layer is a meaningful share of bytes or time. The only
adapter cost with a number attached is the *optional* profile rebuild described above.

---

## 8. Separation of the error sources

| source | evidence | magnitude at 64³ |
|---|---|---|
| numerical diffusion of the interface (spatial) | volume drift grows monotonically to −21.9 %; `min φ` decays 27 %; WENO5 scalar peak decays 5.3 %; unchanged when `Δt` is halved | dominant: L1 return error 3.1e-3, max interface displacement 4.1 h |
| geometric error (interface position, two-sided) | mean forward distance 0.53 h, max 4.1 h; scales ~1st order in the mean, weakly in the max | 8.4e-3 mean / 6.5e-2 max |
| temporal error | CFL halving at 48³ changes φ/WENO5 by < 0.01 %, CIP by 7 % | ≤ 7 % of the CIP scalar error, negligible for φ |
| allocation cost | 1 032 B/node/step, 99.995 % from CIP; 171 GB per 64³ run; peak RSS 1.51 GB | 5–10 % of wall time in GC |
| adapter cost | copy-back 0 B / 0.017 % of a step; velocity 0 B; snapshots 2e-5 % | negligible |
| CIP/CIP-cache | profile rebuild skipped on the loop path (5.4 % of the step saved) | 56 B/node/step if the caller's buffer discipline differs |

Dispersion is not separated from diffusion: the return errors combine both, the peak decay is
a diffusion indicator, and the small negative undershoot of CIP (`min c = −0.023` at `t = 1`
vs `−2e-6` for WENO5) is its dispersion indicator. A dedicated dispersion measurement was not
performed.

---

## 9. Figures

All produced by `visualize.StrongDeformation3D.jl` from the saved arrays; the module contains
no transport and no physics, and every interface drawing is the saved `cut` point set — no
analytic surface is ever drawn.

| figure | content |
|---|---|
| `output/montage_interface_n{48,64}.png` | five-time montage of the true interface point clouds, identical axes and shared colour scale |
| `output/slices_n{48,64}.png` | three orthogonal φ slices at `t = 0.5` (`z=0.56`, `y=0.48`, `x=0.41`) with the zero contour |
| `output/scalar_schemes_n{48,64}.png` | WENO5 (top) vs CIP (bottom) on the `z = 0.56` slice at `t = 0, 0.5, 1`, one colour scale |
| `output/metrics_n{48,64}.png` | `D_sym` (both estimators) and relative volume drift (three estimators) versus time |
| `output/interface_distance_n{48,64}.png` | mean and max, forward and backward interface distances |
| `output/allocation_n{48,64}.png` | per-kernel bytes/node/step and seconds/step |
| `output/convergence.png` | return error versus resolution, L1 and L2, all three quantities |

Readings that carry information the tables do not:

* the montage shows the ellipsoid at `t = 0`, a sheared lobe at 0.25–0.75, and a compact
  blob again at `t = 1` — the return is visually convincing and visibly imperfect at the rim;
* the `yz` slice at `t = 0.5` shows the interface compressed into a thin inclined ellipse,
  i.e. the filament that the grid cannot resolve, which is the origin of the volume drift;
* the scalar slices show a *detached second lobe* at `t = 0.5` (the Gaussian is torn into two
  pieces) for both schemes, with CIP's core visibly brighter — the peak-decay difference of
  §8 seen directly;
* the allocation figure makes the three-order-of-magnitude gap between the transport kernels
  and everything else explicit, including the two zero bars.

Limitation of the plotting: the montage draws points, so its visual density reflects the
number of reconstructed vertices (~17–24 k per snapshot), not a rendered surface. CairoMakie's
`scatter!` inside `Axis3` was verified to work; no volumetric rendering was attempted.

---

## 10. Limitations and what is not shown

* Two resolutions only (48³, 64³): no asymptotic order is claimed, and the observed apparent
  orders (L1 0.9–1.6, L2 0.7–0.9) are reported as data, not as convergence rates.
* One temporal study (48³, CFL 0.4 vs 0.2). CIP's 7 % response is consistent with first-order
  splitting but is not enough to fit an order.
* The exact flow still has the first integral `I` of §3.3. The deformation is a strong
  tangential shear of a surface that cuts across the `I`-shells, not a flow without invariant
  surfaces. Both the `D_sym` values and the interface distances are properties of *this*
  benchmark.
* The interface distance is a sampled, two-sided, level-set-based surrogate with two
  quantified biases (up to 1.56× underestimate from `|∇φ₀|`, and an inflation of the backward
  direction from the 27 % amplitude decay). It is not a Hausdorff distance.
* `cut`'s facet normals are only used for orientation in the volume integral, and the volume
  estimator is the divergence-theorem workaround of the previous campaign, not a subcell
  volume (the `Cut` surface exposes no subcell weights).
* A=2 was not reduced. Its effect is visible as interface diffusion at 48³; an A-sweep would
  be the natural follow-up and was not performed (`A` is a single constant in the script).
* Timings are single-machine, single-thread, unpinned, and vary by ±10–20 % between runs;
  the byte counts do not vary at all.
* No reinitialization was used, so the level-set results include the amplitude decay that a
  reinitialized run would suppress — deliberately, since reinitialization would move the zero
  level set and confound the reversibility measurement.

---

## Appendix — artifacts and reproduction

```text
round2/
  run.StrongDeformation3D.jl          driver: dynamics, loop, diagnostics, allocation study
  visualize.StrongDeformation3D.jl    CairoMakie figures only
  REPORT.md                           this file
  output/verification.txt             Stage A/B evidence (source, flow, ordering)
  output/convergence.csv              return error vs resolution
  output/n48/, output/n64/            params.txt, fields.bin, metrics.csv, metrics.txt,
                                      errors.csv, interface.csv, allocs.csv,
                                      allocs_profile.csv, interface_points_t{0,0.25,0.5,0.75,1}.csv
  output/*.png                        13 figures (six per resolution + convergence)
  output_cfl02/n48/                   the CFL = 0.2 refinement run (same file layout)
```

`fields.bin` holds 17 `Float64` records of `nodes³` in column-major order (first axis
fastest, so `reshape(v, n, n, n)` is indexed `(x,y,z)`): `phi0`, `c0`, then
`phi, c_weno5, c_cip` at `t = 0, 0.25, 0.5, 0.75, 1`. `interface_points_t<S>.csv` are the
reconstructed interface vertices per snapshot, `x,y,z` per line.

```bash
cd ".../EvolvingDomains-Redesign"
julia --project=validation/env validation/VortexReversal3D/round2/run.StrongDeformation3D.jl 48 64
SD3D_CFL=0.2 SD3D_OUT=_cfl02 \
  julia --project=validation/env validation/VortexReversal3D/round2/run.StrongDeformation3D.jl 48
julia --project=validation/env validation/VortexReversal3D/round2/visualize.StrongDeformation3D.jl
```

Resolutions run finest-first, because `Sys.maxrss()` is a process high-water mark: the fine
case's memory reading is only free of the profiler's own footprint if it runs before anything
heavier. The run is deterministic; repeated runs reproduce every CSV byte for byte.
