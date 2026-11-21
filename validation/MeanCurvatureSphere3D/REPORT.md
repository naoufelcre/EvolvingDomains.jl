# MeanCurvatureSphere3D — Final Report

**Benchmark:** three-dimensional mean-curvature flow of a sphere in EvolvingDomains.jl
**Case directory:** `validation/MeanCurvatureSphere3D/`
**Deliverables examined:** `run.MeanCurvatureSphere3D.jl`, `visualize.MeanCurvatureSphere3D.jl`,
`output/n48/`, `output/n64/`

## Scope, evidence, and reading conventions

This report documents what the benchmark actually did and measured. Every numerical value below
is quoted from the artifacts written by the two runs:

- `output/{n48,n64}/timeseries.csv` — recorded at t = 0 and at the three in-loop snapshots;
- `output/{n48,n64}/summary.md` — end-of-run summary of the same run;
- `output/{n48,n64}/interface_t*.csv`, `phi_n*_t*.bin`, `snapshots.csv` — saved geometry;
- `output/n64/fig_*.png|pdf` — figures (see §12).

Claims are labelled where the distinction matters: **measured** (read from those files or
observed in the session record), **analytical** (exact solution or definition), **code fact**
(visible in the implementation), **hypothesis** (possible explanation, not established), and
**open** (unresolved).

No simulation was rerun and no code, manifest, or output was modified for this report. No new
figures were produced.

---

## 1. Executive summary

Mean-curvature flow of a sphere was solved in three dimensions on 48³ and 64³ uniform grids over
0 ≤ t ≤ 0.6, using EvolvingDomains `Grid`, `Discretization`, `Field`, `eulerianTransport!`, `cut`
and `reinitialize!`. The interface at every step was advanced by a velocity built from the
*discrete curvature of the evolving level set* — nodal normals from central differences, mean
curvature from their divergence — never from the analytic radius, which appears only in
diagnostics.

The measured final radii are 0.22060418 (48³) and 0.22080003 (64³) against the exact
0.22113344, i.e. relative errors −2.39 × 10⁻³ and −1.51 × 10⁻³; relative volume errors are
−6.50 × 10⁻³ and −4.12 × 10⁻³. Both improve with refinement, as does the area-weighted standard
deviation of the interface radius (1.01 × 10⁻⁴ → 5.68 × 10⁻⁵), consistent with near-spherical
evolution. The interface radius is already biased inward at t = 0 by the piecewise-linear cut
surface (relative −1.19 × 10⁻³ → −6.66 × 10⁻⁴, a factor 1.79 for a spacing ratio 1.34), so part
of every final error is geometric reconstruction, not evolution.

Interoperability was demonstrated for diagnostics: the public `CutSurface` (facet
connectivity, vertices, normals, CSR point arrays) was converted to plain Julia geometry and
passed to Curvature.jl's `interfaceCurvature`; mean curvature agreed with 2/R to
2 × 10⁻⁴ relative with zero failed evaluations. Curvature.jl did **not** drive the interface:
it returns one value per facet, and no adapter from facet scalars to a background-node velocity
was identified or validated. The velocity instead uses a finite-difference
nodal curvature with a cosine-tapered narrow band, which is not a normal extension. That missing
adapter, a possible native-library lifetime/state effect seen only when many short-lived
`Field`s were cut in one process, and the absence of strongly nonspherical tests leave the
architecture partially validated rather than proven.

## 2. Mathematical model

Domain, level-set convention, and dynamics:

```text
D = [−0.5, 0.5]³
φ(x, t) < 0 inside Ω(t),  φ > 0 outside,  Γ(t) = {φ = 0}
n = ∇φ / |∇φ|                (outward, i.e. pointing out of the φ < 0 phase)
κ = ∇·n                      (sum of principal curvatures; κ > 0 for a convex phase)
V_n = −μ κ,   μ = 0.01       (V_n measured along n)
```

The equivalent level-set equation is

```text
∂_t φ − μ κ |∇φ| = 0      ⇔      ∂_t φ + V_n |∇φ| = 0
```

Initial data and exact solution:

```text
φ₀(x, y, z) = √(x² + y² + z²) − R₀,   R₀ = 0.27
R(t) = √(R₀² − 4 μ t),   0 ≤ t ≤ 0.6
```

**Sign convention and why the sphere shrinks.** With φ negative inside, ∇φ points outward, so
n is the outward normal. On a sphere of radius R every point is strictly convex and
κ = ∇·n = 2/R > 0. The normal velocity is therefore V_n = −μκ < 0, i.e. opposite to n, hence
inward: the sphere contracts. Equivalently, the φ ≤ 0 phase is a ball whose boundary moves
inward, so each level set {φ = c} is advected with normal speed −μκ. The code preserves this
sign at every level: `dynamics(cfg) = kappa -> -cfg.mu .* kappa`, and the velocity is
`v = rate · n` (code fact, `run.MeanCurvatureSphere3D.jl`).

**Derivation of R(t).** For a sphere φ = r − R(t) with r = |x|:

```text
|∇φ| = 1,     κ = 2/R,     ∂_t φ = −Ṙ
∂_t φ − μ κ |∇φ| = 0   ⇒   −Ṙ − 2μ/R = 0   ⇒   Ṙ = −2μ/R
d/dt (R²) = 2R Ṙ = −4μ   ⇒   R(t)² = R₀² − 4μt      (R₀²/(4μ) = 1.8225 > 0.6)
```

Extinction lies at t = 1.8225; the run stops at t = 0.6, where R ≈ 0.819 R₀. The exact enclosed
volume used for comparison is V(t) = 4πR(t)³/3.

## 3. Mathematical and software decomposition

The implementation is six small, explicit numerical pieces orchestrated by one loop; there is no
curvature-flow solver type. Per time step, for the state φⁿ:

1. **Gradient and normals.** `computeNormals!` reads the level-set array and fills a nodal
   `Array{SVector{3,Float64},3}` with n_h = ∇_h φ / max(|∇_h φ|, ε), ε = 0.5.
2. **Discrete curvature.** `computeCurvature!` fills a nodal `Float64` array with
   κ_h = ∇_h·n_h (central differences of the normal field).
3. **Normal rate.** `dynamics(cfg)` maps κ_h to the geometric rate F = −μκ_h.
4. **Cartesian velocity.** `velocityFromRates!` writes ν = χ(|φ|) F n_h into a flat
   `Vector{SVector{3,Float64}}` of length n³, index `i + (j−1)n + (k−1)n²` for `φ[i,j,k]`
   (first axis fastest, as the transport kernel requires).
5. **Transport.** `advanceGeometry!` calls
   `eulerianTransport!(field, velocity, Δt; level_id=:init, cache=WENO5Cache())`.
6. **Reinitialisation** every 10 steps: `reinitialize!(field, cut(field, :init))`, with the
   interface geometry measured immediately before and after.
7. **Snapshot diagnostics** at t ≈ 0, 0.2, 0.4, 0.6: `cut` → facets → radius, volume,
   axis crossings, and an independent Curvature.jl fit.

Actual loop (verbatim from `run.MeanCurvatureSphere3D.jl`, snapshot block elided):

```julia
while t < cfg.tmax - 1e-12
    dt_step = min(dt, cfg.tmax - t)
    steptime = @elapsed begin
        rate = dynamics(cfg)(workspace.curvature)
        velocityFromRates!(workspace.velocity, phi, rate, workspace.normals, h, cfg)
        advanceGeometry!(field, workspace.velocity, dt_step, cache)
        t += dt_step
        step += 1

        if cfg.reinit_every > 0 && step % cfg.reinit_every == 0
            before = interfaceGeometry(field, h)
            gb = gradientStats(phi, h)
            reinitialize!(field, cut(field, :init))
            after = interfaceGeometry(field, h)
            ga = gradientStats(phi, h)
            push!(reinit_shifts, (radius=abs(after.radius - before.radius),
                                  volume=abs(after.volume - before.volume)))
            push!(reinit_grad_before, gb)
            push!(reinit_grad_after, ga)
        end

        clamped += computeNormals!(workspace.normals, phi, h, cfg.eps_gradient)
        nonfinite += computeCurvature!(workspace.curvature, workspace.normals, h)
    end
    push!(step_times, steptime)
    # ... snapshot diagnostics when t reaches the next recorded time ...
end
```

Note the order: normals and curvature are recomputed *after* the transport, so the next step's
velocity uses the κ of the state at the start of that step; the t = 0 values are computed once
before the loop.

### Division of responsibility

| Operation | Provider |
|---|---|
| Cartesian grid, level → grid mapping, dense value arrays | EvolvingDomains `Grid`, `Discretization`, `Field` |
| WENO5 + SSP-RK3 advection of φ, RK workspace reuse | EvolvingDomains `eulerianTransport!`, `WENO5Cache` |
| Fast-sweeping reinitialisation preserving cut zero crossings | EvolvingDomains `reinitialize!` |
| Piecewise-linear interface triangles, their vertices and unit normals, cell grouping | EvolvingDomains `cut` → `Cut` → `CutSurface` |
| Nodal gradient, clamped normals, discrete mean curvature, rate, narrow-band velocity | validation code (`run.MeanCurvatureSphere3D.jl`) |
| Independent per-facet interface curvature | Curvature.jl `interfaceCurvature` |
| Radius, volume, axis crossings, sphericity, tables, figures | validation code + `visualize.MeanCurvatureSphere3D.jl` |

The level set is stored as an ordinary three-dimensional `Array{Float64,3}`; normals and
velocity are `StaticArrays.SVector{3,Float64}` collections; curvature is an ordinary
`Array{Float64,3}`. No mathematical operator is hidden inside a custom type: the operators are
central differences, a pointwise clamp, a divergence, and a scalar taper, each written out in
the validation file.

## 4. Spatial discretization and time integration

### Normals and curvature

```text
n_h = ∇_h φ / max(|∇_h φ|, ε),        ε = 0.5
κ_h = ∇_h · n_h
```

Both derivatives are centred differences on the uniform grid. The gradient is
(φ[i+1] − φ[i−1]) / (2h) per axis, clamped at the faces (the clamped-index convention the
package also uses for WENO ghost values); at a boundary node this reduces to a half-magnitude
one-sided difference. The box boundary is irrelevant to the physics here: R ≤ 0.27 and the band
ends at R + 4h ≤ 0.355 < 0.5, so no interface velocity is ever computed near a face.

The regularisation is a *magnitude* clamp on the normal, not a smoothing of φ: where
|∇_h φ| < ε the direction is retained but |n_h| ≤ 1. The justification is scale matching
rather than tuning: a reinitialised level set has |∇φ| = 1 at the interface, and ε = 0.5 leaves
the interface untouched while keeping the operator bounded where the level set degenerates (near
|∇φ| = 0, e.g. at the centre of a radially symmetric φ). The clamp is counted and reported:
summed over the initial evaluation and all steps, 2068 node evaluations were clamped in the 48³
run (out of 17.69 million) and 3354 in the 64³ run (out of 74.97 million).

κ_h = ∇_h·n_h is a second-order-consistent approximation of the mean curvature of the
interpolated level set, but only where the normal field is a surface-normal field. Far from Γ,
div n_h is a property of the numerical level set, not of the interface, and is not used: the
velocity is tapered there.

### Narrow-band velocity: the cosine taper

```text
ν(x) = χ(|φ(x)|) · (−μ κ_h(x)) · n_h(x)

χ(d) = 1                                  for d ≤ 1.5 h
χ(d) = cos²( (π/2)(d − 1.5h)/(4h − 1.5h) ) for 1.5 h < d < 4 h
χ(d) = 0                                  for d ≥ 4 h
```

This is a C¹ narrow-band field. On Γ, χ = 1 and |∇φ| ≈ 1, so the transported level set sees the
intended interface law ∂_t φ = μκ|∇φ| exactly (to the accuracy of κ_h); the taper only freezes
the inert part of the level set so that meaningless κ_h values outside the band cannot advect
φ. It is **not** a normal extension of interface data: a normal extension would assign to each
node the interface value carried along ∇φ, i.e. a value independent of distance from Γ; the
taper instead multiplies the *local* nodal rate by a decreasing scalar (≈ 0.90 at |φ| = 2h,
≈ 0.35 at 3h). The two operations answer different questions and are not interchangeable (see
§7).

**Measured interface displacement per step** (from the recorded radius change, 0.0488666 over
the whole run): 3.05 × 10⁻⁴ ≈ 0.0143h at 48³ and 1.71 × 10⁻⁴ ≈ 0.0108h at 64³; the largest
instantaneous value, μκΔt with κ ≈ 9.05, is 1.6 × 10⁻² h and 1.2 × 10⁻² h. The full-speed zone
is therefore ~100 steps wide, and the interface never enters the taper. **Code-fact caveat:** the
docstring of `velocityFromRates!` states "it moves ~0.15 h per step", which overstates the
measured motion by roughly ten; the argument is unaffected (the zone is much wider than the
motion), but the comment is inaccurate and should be corrected in a future edit.

### Time step

```text
Δt = C h² / μ,    C = 1/12
```

| Quantity | 48³ | 64³ |
|---|---:|---:|
| Grid spacing h | 0.02128 | 0.01587 |
| Nominal time step | 3.772e-3 | 2.100e-3 |
| Number of steps | 160 | 286 |
| Final time | 0.6 | 0.6 |

For the linear diffusion operator, the 3D seven-point von Neumann bound is μΔt/h² ≤ 1/6; the
chosen C = 1/12 halves that, and the run records μΔt/h² = 0.0833. This is the correct *scaling*
for a curvature-driven interface, whose parabolic stability scale is Δt ∼ h²/(μ κ² h²) — at
fixed κ, Δt = O(h²/μ). It is **not** a rigorous stability proof for the complete nonlinear
discretisation: WENO5 reconstruction, the SSP-RK3 stages, the nonlinear div n operator and the
periodic reinitialisation all modify the effective linearisation, and no von Neumann analysis or
energy estimate was carried out. What the run supports is weaker and should be stated as such:
with C = 1/12, over 160 and 286 steps, no non-finite curvature, no blow-up, and no spurious
interface mode was observed. The finiteness checks (`count(!isfinite, κ_h)`, plus the transport's
own validation of every velocity component) are recorded as zero.

## 5. Radius, volume, and geometric convergence

| Quantity | 48³ | 64³ |
|---|---:|---:|
| Numerical final radius | 0.22060418 | 0.22080003 |
| Exact final radius | 0.22113344 | 0.22113344 |
| Absolute radius error | −5.293e-4 | −3.334e-4 |
| Relative radius error | −2.39e-3 | −1.51e-3 |
| Relative volume error | −6.50e-3 | −4.12e-3 |
| Standard deviation of interface radius | 1.01e-4 | 5.68e-5 |
| Maximum axis-crossing error | 5.94e-4 | 3.77e-4 |

### How the radius is measured

At each recorded time the run takes the public `cut` of the current `Field`, converts the
boundary triangles to plain geometry, and computes per-facet area a_q and centroid c_q. The
reported radius is the **area-weighted mean of |c_q|**:

```text
R_num = Σ a_q |c_q| / Σ a_q,     σ_R² = Σ a_q (|c_q| − R_num)² / Σ a_q
```

σ_R is the sphericity measure: it is the area-weighted spread of the triangulated interface's
distance from the centre, so uniform shrinkage contributes nothing to it and a triaxial
deformation would show immediately. The enclosed volume uses the divergence theorem on the same
closed triangulation,

```text
V_num = (1/3) Σ a_q (c_q · n_q)
```

with the facet normals supplied by `CutSurface` (already oriented outward from the φ < 0 phase).
The axis radii are a second, independent estimator: the trilinear interpolant of φ on each of the
six rays along the coordinate axes is bracketed by bisection to 80 iterations. All six crossings
were identical at every recorded time (to the printed 10⁻¹⁰), as symmetry requires; their common
value is the "axis" radius.

### Reconstruction bias versus evolution error

The interface is represented by the package's six-tetrahedron piecewise-linear cut. For a
strictly convex body that polyhedron is inscribed, so its area-weighted centroid radius is
*systematically* smaller than the sphere it replaces. The t = 0 row of `timeseries.csv` measures
exactly that, because at t = 0 φ is the exact signed distance of the exact sphere:

| t = 0 measured | 48³ | 64³ |
|---|---:|---:|
| radius | 0.26967748 | 0.26982007 |
| radius relative error | −1.1945e-3 | −6.6640e-4 |
| volume relative error | −3.0976e-3 | −1.7292e-3 |

The ratio between grids is 1.792 (radius) and 1.791 (volume): consistent with the second-order
trend expected from an inscribed reconstruction, since the spacing ratio squared is
(0.0212766/0.0158730)² = 1.797. This is evidence for a reconstruction bias of order h², not a
proof of asymptotic order from two grids.

The recorded final errors therefore mix two contributions. Separating them with the measured
t = 0 bias:

| | 48³ | 64³ |
|---|---:|---:|
| radius error drift over the run (final − initial) | −2.07e-4 | −1.53e-4 |
| as a fraction of the exact shrinkage 0.0488666 | 0.42 % | 0.31 % |
| bias rescaled to the final radius, assuming the relative bias scales as (h/R)² | −3.94e-4 | −2.20e-4 |
| measured final error | −5.29e-4 | −3.33e-4 |
| residual not explained by the rescaled bias | −1.36e-4 | −1.14e-4 |
| residual as a fraction of the shrinkage | 0.28 % | 0.23 % |

The rescaling assumption is itself supported by the t = 0 scaling across the two grids; the
residual shrinks only ~15 % between them, so **no order can be inferred for it** from this data.
What can be said is bounded: the flow reproduces the exact shrinkage to ≈ 0.3 % of its own
magnitude, and roughly three quarters of the final radius error at both resolutions is the
geometric reconstruction bias rather than accumulated dynamics.

### Evidence for near-spherical evolution

σ_R is 1.02 × 10⁻⁴ at t = 0 and 1.01 × 10⁻⁴ at t = 0.6 (48³); 5.80 × 10⁻⁵ and 5.68 × 10⁻⁵
(64³). It is essentially constant and slightly decreasing, i.e. the same order as the static
discretisation pattern of the cut surface and not growing. Over a run in which the radius
decreases by 2.3h (48³), no measurable nonspherical mode develops. This is a statement about
*this* estimator at *these* resolutions: σ_R is dominated by the faceting pattern, so it bounds
deformation above but cannot resolve very small axisymmetric deviations. The axis-radius
estimator is even less sensitive, because all six rays are symmetry-equivalent.

## 6. Curvature accuracy and Curvature.jl interoperability

### The exposed geometry

`cut(field, :init)` returns an owned `Cut`; its `boundary::CutSurface{3,Float64}` carries

- `facet_to_bgcell` — one entry per boundary facet, sorted by background cell (the grouping
  guarantee);
- `point_ids::Vector{Int32}` and `point_pointers::Vector{Int32}` — CSR connectivity: facet `q`
  owns `point_ids[point_pointers[q] : point_pointers[q+1]−1]`, three vertices in 3D;
- `physical_points::Vector{NTuple{3,Float64}}` — mesh-coordinate vertices;
- `reference_points` — the same vertices in the parent cell's [0,1]³ frame (unused here);
- `normals::Vector{NTuple{3,Float64}}` — one unit normal per facet, outward from φ < 0.

The conversion in `interfaceGeometry` is deliberately mechanical and asserts the arity the
producer promises:

```julia
a1 = Int(boundary.point_pointers[q])
a2 = Int(boundary.point_pointers[q + 1]) - 1
a2 - a1 + 1 == 3 || error("facet $q has $(a2 - a1 + 1) vertices, expected 3")
p1 = SVector{3,Float64}(boundary.physical_points[Int(boundary.point_ids[a1])])
p2 = SVector{3,Float64}(boundary.physical_points[Int(boundary.point_ids[a1 + 1])])
p3 = SVector{3,Float64}(boundary.physical_points[Int(boundary.point_ids[a1 + 2])])
facets[q] = (p1, p2, p3);  normals[q] = SVector{3,Float64}(boundary.normals[q])
```

Facet order is preserved, so `interfaceCurvature(...)` returns κ values indexable by the same
`q`. The call is the documented one-shot public entry point:

```text
κ_facet = interfaceCurvature(facets, normals; spacing = h, origin = −0.5, radius = 4h)
```

`spacing` and `origin` are the uniform lattice's length scale and lower corner (a scalar here
because the grid is cubic); `radius = 4h` selects the fit window. The normalisation of the
normals is Curvature.jl's business — it renormalises on entry — and the signs agree by
construction: both libraries orient outward from the φ < 0 phase, so both report κ > 0 for this
convex sphere.

Facet counts per snapshot (from `timeseries.csv`): 18 096 → 12 024 at 48³; 32 568 → 21 840 at
64³.

### Measured curvature

| Quantity | 48³ | 64³ |
|---|---:|---:|
| Mean curvature, Curvature.jl (area-weighted) | 9.04250 | 9.04668 |
| Analytical mean curvature 2/R_exact | 9.04431 | 9.04431 |
| Finite-difference band mean | 9.04874 | 9.04434 |
| Failed Curvature.jl evaluations (NaN) | 0 | 0 |

Relative errors: Curvature.jl −2.0 × 10⁻⁴ (48³) and +2.6 × 10⁻⁴ (64³); finite-difference band
mean +4.9 × 10⁻⁴ and +3 × 10⁻⁶. The sign flip between grids for the Curvature.jl mean is a
reminder that two resolutions cannot establish a convergence rate for this quantity.

Two different numerical objects are being compared, and the difference matters:

- **Curvature.jl** fits a local polynomial through the faceted interface inside a ball of
  radius 4h, then applies its contour correction, and returns one value *per facet*, i.e. a
  quantity attached to Γ. Its area-weighted mean is therefore an estimate of the surface average
  of κ over the discretised interface; its area-weighted spread at the final time is 0.0468
  (48³) and 0.0445 (64³), about 0.5 % of κ.
- **The finite-difference operator** returns κ at every *node* of the background grid from the
  nodal normal field. Its band statistics here are taken over all nodes with |φ| ≤ 0.8h, a
  slightly wider shell than the interface itself, and include nodes where div n_h is not the
  surface curvature; its spread is correspondingly larger, 0.421 (48³) and 0.304 (64³).

Comparing the two estimates *pointwise* is the sharper test, and the code does it by trilinear
interpolation of the nodal κ_h to each facet centroid. The area-weighted RMS discrepancy,
normalised by |2/R_exact|, is

| | 48³ | 64³ |
|---|---:|---:|
| RMS(κ_h − κ_Curvature.jl) at facet centroids / (2/R_exact) | 6.79e-3 | 5.06e-3 |
| in absolute curvature units | 0.061 | 0.046 |

This is an order of magnitude larger than the mean bias of either estimator: the disagreement is
**local and distributional**, not a uniform offset. That is exactly why the mean alone cannot
certify the curvature field — a nodal operator can have a very accurate area-weighted mean while
misplacing curvature locally by several percent of κ.

### Which reference each figure uses

The two references are different and both are legitimate:

```text
2/R_num(0.6)  = 2/0.22080003 = 9.058   (64³)      ← the measured interface's own scale
2/R_exact(0.6) = 9.04431                            ← the analytical sphere at the recorded time
```

`fig_curvature.png` compares two different things in its two panels:

- **left panel** (histogram): the *per-facet* Curvature.jl values at the final snapshot, with a
  dashed line at **2/R_num**, computed in `figureCurvature` as
  `rbar = Σ a_q|c_q| / Σ a_q;  kexact = 2/rbar`. The comparison asks: *is the curvature
  distribution centred on what the measured interface's own scale would imply?* The figure
  therefore shows κ ≈ 9.046 against 9.058 for 64³, i.e. the facet curvature is ~0.13 % below the
  measured interface's 2/R.
- **right panel** and the whole `timeseries.csv` `kappa_*` columns compare the area-weighted
  means with **2/R_exact(t)**. Those ask: *do the curvature estimators track the analytical
  solution?*

Both are meaningful; conflating them would turn a 1.3 × 10⁻³ discrepancy into a 2 × 10⁻⁴ one
or vice versa without any physical change.

## 7. Interface velocity extension: an unresolved integration point

The intended composition, and the one that would make the curvature library part of the
dynamics, is

```text
CutSurface → facet curvature (κ_Γ) → normal speed (−μκ_Γ) → background velocity → φ transport
```

The pieces that exist today:

1. **Mathematical input/output an extension operator needs.** Input: a scalar (or vector)
   sampled on Γ — here κ per facet, together with the facet position and normal — plus the
   background grid and the level set. Output: a value at every background node such that
   (i) on Γ the resulting velocity satisfies v·∇φ = −μκ_Γ, and (ii) away from Γ the field is
   constant along the characteristics of ∇φ (a normal extension), so that the level sets move
   together and the distance property is not destroyed. The two clauses are separate: (i) is a
   boundary/interface condition, (ii) is an extrapolation rule.
2. **What Curvature.jl provides.** `interfaceSamples`/`interfaceCurvature` return per-facet
   data indexed by facet (positions = centroids, normals, measures, curvatures, vertices,
   `cell_to_facets`, fit radius) and the queries `getCurvature`, `curvatureAt(s, x)`,
   `interfaceTangents`, `interfaceTangentPlanes`. `curvatureAt` does evaluate the interface
   curvature nearest an arbitrary query point, which is a pointwise *sampling* of κ_Γ at grid
   nodes; it does not extend along ∇φ, and it gives no guarantee about clause (i) at a node
   whose φ-gradient direction differs from the nearest facet's normal.
3. **What EvolvingDomains requires.** `eulerianTransport!` accepts
   `Vector{SVector{3,Float64}}` with exactly one entry per node, first axis fastest, finite
   components. There is no interface-level or facet-level input form.
4. **The adapter that is missing from this workflow.** A projection from facet-associated
   scalars to nodal values plus a normal extension, or an equivalent formulation in which the
   facet curvature drives the transport without being materialised as a nodal velocity (e.g. a
   ghost-fluid or cut-cell flux using the facet data directly, with the level set advanced by
   the resulting interface motion). Neither was implemented or validated here.
5. **Why the cosine taper is not that adapter.** The taper multiplies the *nodal* rate by
   χ(|φ|) and leaves the direction at n_h. It is defined from the local φ, not from Γ's data;
   its value changes with distance, whereas a normal extension is constant along ∇φ; it cannot
   carry a facet quantity at all. It is a well-posedness device for the nodal operator, not a
   transport of interface information.

**Open question, explicitly not resolved.** During implementation no suitable *public*
extension operator was identified, and the campaign stopped there. That is not evidence that no
such operator exists: the maintainer recalls earlier functionality, it may be reachable under a
different name, may be internal, or may have been removed in the redesign. A concrete, unexamined
pointer exists in this checkout — `validation/archive/legacy_numerics/Transfer/ExtensionOperator.jl`
— which was *not* read for this report and whose API must not be assumed. Determining whether
existing mathematical objects (transport maps, the level-set gradient field, `curvatureAt`
sampling plus a normal-extension sweep) can provide the transformation is investigation A of the
second campaign; only after that determination should anyone conclude that a new public API is
needed.

Consequence for this campaign: the interface is advanced by `computeCurvature!`, a
finite-difference nodal operator written in the validation file, and Curvature.jl serves purely
as an independent estimator of the same surface property (§6). The two roles must not be
conflated in any claim about the architecture.

## 8. Reinitialization and zero-level preservation

`reinitialize!` is called every 10 steps (`reinit_every = 10`; 16 calls at 48³, 28 at 64³), never
every step, and each call is measured: the cut-surface radius and enclosed volume are computed
immediately before and after, and the mean |∇φ| over the band |φ| ≤ 1.5h is recorded on both
sides.

```text
reinitialize!(field, cut(field, :init))
```

Mathematically the routine solves the eikonal reinitialisation to first order by fast sweeping,
seeding the nodes of explicit cut cells so that the represented zero crossings are held fixed.
Its purpose here is not to move Γ but to undo the compression/stretching of φ that advection with
a spatially varying normal speed produces, so that κ_h = div n_h is evaluated on a level set with
|∇φ| ≈ 1.

**Measured effect** (`summary.md`):

| | 48³ | 64³ |
|---|---:|---:|
| calls | 16 | 28 |
| max \|ΔR\| per call | 5.55e-17 | 8.33e-17 |
| max \|ΔV\| per call | 1.39e-17 | 2.78e-17 |
| mean \|∇φ\| in the band, before → after (averaged over calls) | 0.98452 → 0.99692 | 0.99131 → 0.99832 |

The interface displacement is at round-off, which is what the documented "zero crossings remain
fixed" guarantee predicts: the cut vertices are seeded at their current positions, so the
triangulated surface is unchanged as a *set*, and both the radius and the volume estimator are
functionals of exactly that surface. This supports the narrow claim: **reinitialisation does not
move the measured interface**, and it improves the level-set gradient in the band by 1.2 %
(48³) and 0.7 % (64³) on average.

It does **not** support stronger claims. The result does not show that φ is an exact signed
distance function afterwards — `|∇φ| = 0.997` and the routine's own first-order sweeping error
remain — nor that reinitialisation is free of long-term bias (its documented tendency to amplify
a small cut-volume bias was not observable in 16–28 calls at this accuracy, but was not
specifically tested). What *is* observed is that omitting reinitialisation is noticeably worse:
an exploratory run with `reinit_every = 0` (prototype, not part of the recorded artifacts) ended
with a larger radius error and a markedly rougher curvature band than the periodic schedule. The
recorded artifacts contain only the every-10-steps case, so that comparison is reported as a
prototype observation, not as a benchmark result.

## 9. Performance

| Quantity | 48³ | 64³ |
|---|---:|---:|
| Total runtime (as recorded) | 15.3 s | 41.4 s |
| Reported mean time per step | 44.2 ms | 92.2 ms |
| Failed curvature evaluations | 0 | 0 |

**Coverage of those numbers (important).** The recorded `elapsed` spans the time loop, which
includes transport, the reinitialisation cuts every 10 steps, and the *three* in-loop snapshot
diagnostics (each of which performs a `cut`, a Curvature.jl fit and the axis bisections). It
excludes the t = 0 diagnostic, the CSV/raw writes, and the five figure renders. The string in
`summary.md` says "(incl. snapshots and figures)"; the "figures" part is imprecise and should be
read as "post-processing". Splitting the recorded numbers:

- transport loop: 160 × 44.18 ms = 7.1 s (48³); 286 × 92.24 ms = 26.4 s (64³);
- in-loop snapshot diagnostics: ≈ 8.3 s (48³), ≈ 15.0 s (64³), i.e. ≈ 2.8 s and ≈ 5.0 s per
  snapshot, dominated by the Curvature.jl fits at 4h window on 12 000–22 000 facets;
- per-node transport cost: 0.40 µs (48³), 0.35 µs (64³), the finer grid being slightly more
  efficient per node.

Allocations, as code facts rather than measurements: the normals array, the curvature array, the
velocity vector and the `WENO5Cache` workspace are allocated once and reused; each step
materialises the rate field `rate = dynamics(cfg)(workspace.curvature)` as a fresh N³
`Float64` array (0.88 MB at 48³, 2.10 MB at 64³) because `dynamics` returns an allocating
closure. No `@allocated` accounting was recorded, so the total allocation volume is not known.

What is done every step and what is not: per step the code computes normals, curvature, rate and
velocity (four O(N³) passes) and runs one WENO5/SSP-RK3 advection; it does *not* call
Curvature.jl. The expensive facet fit runs four times per simulation (t = 0, 0.2, 0.4, 0.6).
This is the performance shape the benchmark was designed around, and it means the measured step
cost is not representative of a Curvature.jl-driven simulation.

**Estimated cost of the alternative** (estimate, not benchmarked): driving transport from
Curvature.jl would require one facet fit plus an extension per step. At the measured ≈ 2.8 s per
snapshot fit (48³) and at least one extension sweep, a 160-step run would cost minutes rather
than 15 s; the 64³ case likewise grows from tens of seconds to tens of minutes. Those figures are
projections from the snapshot timings, not measurements; the extension cost is unknown because no
extension was implemented.

Plausible improvements, none implemented: precompute a narrow-band node mask once per
reinitialisation instead of testing |φ| at every node every step; extend the velocity by a
fast-sweeping normal extension (cheap, O(N³) with a few sweeps) instead of the cosine taper, which
would make the band edge physical rather than artificial; avoid the per-step rate allocation by
fusing rate and velocity construction; and use `Curvature.jl`'s existing `InterfaceSamples`
object rather than the one-shot wrapper if several window radii must be compared, since window
changes require a rebuild.

## 10. Possible native backend state or lifetime issue

**Status: potential high-priority reproducibility issue; not diagnosed. No cause is asserted.**

Observation from exploratory testing (before the production script was fixed), recorded here
because it bears on how reliable the cut/transport combination is in long-lived processes:

- Two prototype scripts ran several complete simulations *in one Julia process*, each
  simulation constructing a fresh `Grid`, `Discretization` and `Field`, and calling `cut` both
  for reinitialisation and for diagnostics on fields that were subsequently abandoned.
- In both scripts, the third or fourth such inner run — the failing index differed between the
  two invocations — aborted inside `eulerianTransport!` with
  `ArgumentError: Field evolution velocity values must be finite`. The prototype computed no
  further diagnostics at that point, so the non-finite entries were not attributed to a specific
  array element or region.
- A separate script that ran the *same* configuration three times in one process was stable
  and produced bit-identical results each time; every production run — one simulation per
  process — completed, and the 48³ numbers reproduced exactly across separate invocations.
- Production code differs from the prototypes in one relevant respect: it reuses a single
  `Field` for every `cut` (reinitialisation and all diagnostics), whereas the prototypes built a
  new short-lived `Field` on a new `Grid` inside the diagnostic helper for every call.

Candidate explanations — all compatible with the above, none distinguished: (i) an artifact of
the prototype setup (the new-field-per-call diagnostic helper, or stale references held by the
prototype's closures); (ii) Julia object lifetime / garbage collection interacting with arrays
the native cutter returned or adopted; (iii) reuse of native scratch state across cuts on
different mesh identities; (iv) a genuine non-finite value produced by the prototype's own
physics on that particular configuration, unrelated to memory. The observation therefore does
**not** establish a native-memory bug, and it should not be reported as one.

Recommended handling: investigation E of §13 — a minimal reproduction that varies one factor at
a time (number of Fields, whether `cut` is called on each, whether the previous Field is still
referenced), with the failing array dumped before transport. Until that exists, the safe usage
pattern observed here is one simulation per process or, at minimum, one long-lived `Field` with
all cuts taken from it.

## 11. Architectural assessment

EvolvingDomains states two design aims — minimalism and composability — and this benchmark is a
reasonable probe of both, with a clear boundary between what it confirmed and what it did not.

### Minimalism

The package supplied exactly four data structures (`Grid`, `Discretization`, `Field`, and the
`Cut` snapshot) and four operations (`eulerianTransport!`, `cut`, `reinitialize!`, plus the
multilevel transfer API, unused here). No curvature-flow model class, no solver framework, no
configuration object, and no plotting was needed or used. The simulation's own state is three
arrays (φ, n, κ), a flat velocity vector, and a `WENO5Cache`; the mathematical decomposition —
normals, curvature, rate, velocity, transport — is visible as five functions in the validation
file. The loop above reads as the mathematics it implements, with one caveat: node ordering is
implicit in a flat array, so the code must carry an explicit self-check that no type can enforce.

### Composability

The experiment produced genuine cross-package traffic: EvolvingDomains' public `CutSurface` →
plain Julia geometry → Curvature.jl's public `interfaceCurvature`, with no changes to either
package and no access to package internals. This is stronger evidence than comparing two
transport schemes inside one package, because the data crossed a boundary the packages never
agreed on: neither library knows the other's types, and the adapter is ten lines of index
arithmetic in the validation file.

But the claim must be scoped precisely. What was demonstrated:

- **geometric diagnostics** — facet data produced by EvolvingDomains consumed by an independent
  curvature library, at four times, at two resolutions, with zero failed evaluations and mean
  agreement at the 10⁻⁴ level;
- **evolution with an independently implemented nodal operator** — the transport itself ran on
  a curvature field computed in the validation file, not on Curvature.jl output;
- **one-array interchange** — level set in, level set out, with no hidden geometric state.

What was *not* demonstrated: **Curvature.jl driving the interface velocity**. That pipeline still
lacks the facet-to-node adapter of §7, so the architecture is validated as an open
data-exchange architecture, not yet as a "curvature library plugged into a transport engine"
architecture.

### Interface assumptions the glue code had to make

All of these are code facts, and each is a place where a silent mismatch could hide:

1. **Velocity ordering.** One `SVector{3,Float64}` per node, first axis fastest, matching the
   column-major order of `Array{Float64,3}`. Enforced by an exact self-check: a linear φ
   advected by a constant velocity must shift by exactly v·Δt at an interior node; the check
   reports |error| = 0.0 and the program aborts if it exceeds 10⁻¹⁰. A transposed or
   wrong-stride velocity would fail this test.
2. **Facet arity and CSR layout.** Every boundary facet has exactly three vertices; pointers are
   one-based with `nfacets+1` entries; `normals[q]` belongs to facet `q`. Asserted at conversion
   time.
3. **Normal orientation.** Both packages use "outward from φ < 0", so signs agree with no
   conversion; had EvolvingDomains used the opposite convention, every curvature sign and hence
   the direction of motion would flip — this is the single most dangerous unchecked convention,
   and it is documented rather than tested.
4. **Uniform cubic grid.** The trilinear sampler and the axis bisection assume one spacing for
   all three axes and `origin = first(axis)`; the case is cubic, but the code would silently
   mis-sample on an anisotropic grid.
5. **Single level, `:init`.** Every call passes `level_id=:init` explicitly.
6. **`reinitialize!` needs a `Cut` of the same `Field`.** No cross-mesh use is attempted.
7. **Native library availability.** This checkout publishes no cut artifact, so the documented
   developer build (`scripts/build_cxx_bindings.sh`, output under `$TMPDIR`) must be pointed at,
   either through `EVOLVINGDOMAINS_CUTCELLS_LIBRARY` or the fallback search in the script. The
   ABI must match the source: an older build from the same directory produced a different
   library and was superseded for these runs.

### Balanced verdict

The architecture did what this benchmark demanded: ordinary arrays in, geometry and diagnostics
out, with small explicit operators and no framework intrusion, and with cross-package data
exchange that worked on the first attempt once the conventions were read. It did not *by itself*
provide the one adapter that a curvature-library-driven flow needs, and the native-library
lifetime question is open. Calling the architecture "fully validated" would overstate the
evidence; "validated for geometric interoperability and for array-level orchestration, with a
documented gap at facet-to-node extension" is the accurate summary.

## 12. Existing figures

The 64³ run retained the complete figure set; the 48³ run retained two of the five (its plotting
stage was truncated after `fig_slices` — the simulation, `timeseries.csv`, and `summary.md` for
48³ are complete). No figure was regenerated for this report.

- `output/n48/fig_interface3d.png` — 48³ interface point clouds at t = 0.00, 0.41, 0.60.
- `output/n64/fig_interface3d.png` — 64³ interface point clouds at t = 0.00, 0.40, 0.60.
- `output/n64/fig_slices.png` — mid-plane slices with the zero contour, four times, three planes.
- `output/n64/fig_radius.png` — radius history and absolute radius error.
- `output/n64/fig_volume.png` — relative volume error and radius error / sphericity.
- `output/n64/fig_curvature.png` — facet-curvature distribution and mean-curvature history.

![Interface at 48³](output/n48/fig_interface3d.png)

*48³: three snapshots of the cut-surface centroids, coloured by (|c_q| − R̄)/h with a fixed
colour range of one cell. The sphere shrinks visibly while the colour pattern stays at the
±0.01h level of the static faceting pattern; the red/blue banding is the discrete cut geometry
(centroids of Kuhn-tetrahedron facets lie at systematically different distances depending on
orientation), not a growing instability.*

![Interface at 64³](output/n64/fig_interface3d.png)

*64³: the same visualisation at the finer grid, where the faceting pattern is finer and the
sphericity is halved (σ_R = 5.7 × 10⁻⁵). The interface remains visibly spherical at t = 0.6, with
R = 0.2208 against the initial 0.27.*

![Zero-level slices](output/n64/fig_slices.png)

*A dependable 2D view of the same field: φ/h on the z = 0, y = 0 and x = 0 mid-planes, with the
zero contour drawn. The colour scale is saturated at ±2h so the near-interface structure is
visible; the contour is smooth and single-valued at every recorded time, which rules out spurious
secondary interfaces and visible grid-scale wrinkling.*

![Numerical and exact radius](output/n64/fig_radius.png)

*Left: the area-weighted cut-surface radius at the four recorded times lies on the exact curve
√(R₀² − 4μt); the axis-crossing estimator coincides with it at the plotted scale. Right: the
absolute errors on a logarithmic axis, ≈ 2 × 10⁻⁴ → 3 × 10⁻⁴, growing slowly and monotonically —
the signature of a small residual rate error, not of a systematic offset appearing during the
run.*

![Volume error](output/n64/fig_volume.png)

*Left: relative enclosed-volume error, from −0.17 % at t = 0 (pure reconstruction bias of the
inscribed polyhedron) to −0.41 % at t = 0.6. Right: the radius error in cells alongside the
sphericity σ_R/h: the radius error grows from −0.011h at t = 0 to −0.021h at t = 0.6, i.e. by a
hundredth of a cell over the whole run, while σ_R/h stays flat.*

![Curvature](output/n64/fig_curvature.png)

*Left: area-weighted distribution of the per-facet Curvature.jl curvature at t = 0.6, with the
dashed line at 2/R_num = 9.058 (§6 explains why that reference differs from 2/R_exact = 9.044).
The distribution is concentrated within ±1 % of its mean and has a tail to higher values.
Right: area-weighted means of both curvature estimators against 2/R_exact(t); the three curves
are indistinguishable at this scale, and the local RMS disagreement (0.51 %) is much larger than
their mean offsets.*

## 13. Proposed second-campaign investigations

None of the following was performed. For each, the uncertainty it would resolve is stated.

**A. Recover or identify the interface extension operator.** Inspect the public API and the Git
history for a previous or internal extension/projection capability; examine the unread
`validation/archive/legacy_numerics/Transfer/ExtensionOperator.jl` and determine whether it (or a
successor) accepts facet-associated scalars and produces a background-grid array. This resolves
the single architectural question of this campaign: whether the facet→node transformation already
exists as a mathematical object or must be written. Do this *before* proposing new public API.

**B. Curvature.jl-driven evolution.** Given an extension mechanism, replace the finite-difference
curvature input with Curvature.jl facet curvature, without touching package internals, and compare
accuracy and cost against the present run. This resolves whether the diagnostic pathway is also a
viable evolution pathway — the distinction §6 and §7 keep separate — and quantifies the price of
the facet fit per step.

**C. Perturbed geometry.** Evolve an ellipsoid or a low-amplitude spherical-harmonic perturbation
and check that non-uniform curvature produces the expected rounding (axis ratios relaxing, high-κ
regions moving faster). This resolves the central scientific limitation of the present case: a
sphere is the one geometry where "curvature flow" and "constant inward normal speed" are hard to
distinguish, so C is the test that the operator is genuinely κ-driven. It also probes stability
away from spherical symmetry, where the λ ∝ h²/μ argument is weakest.

**D. Higher-resolution convergence.** Add at least two more resolutions (or a fixed-grid
time-refinement study) to separate reconstruction bias (h², measurable at t = 0), curvature
discretisation error, transport error, and reinitialisation error. This resolves the fact that
the present residual (0.23–0.28 % of the shrinkage) is not identifiable from two grids, and would
give the first actual measured order for the evolved solution rather than a two-point trend.

**E. Native backend lifetime investigation.** Build a minimal reproducible example that varies,
one factor at a time: number of `Field`s/`Grid`s, whether `cut` is called on each, whether the
previous `Field` is still reachable, GC pressure. Record the precise failing array and stack, and
determine whether the non-finite velocity comes from stale/adopted native arrays, cache reuse
keyed by mesh identity, or a numerical blow-up in that particular configuration. This resolves the
reproducibility risk of §10 and determines whether the one-Field usage pattern is a requirement or
a coincidence. It must be isolated before any long-running multi-case process is trusted.

**F. Extension performance.** Compare candidate facet→node transformations — nearest-facet
projection (`curvatureAt`-style sampling), weighted local fits, fast-sweeping normal extension,
closest-point evaluation with a spatial index — on cost and on whether they preserve the interface
condition v·∇φ = −μκ_Γ. This resolves the performance side of B: §9's estimates suggest a
per-step facet fit could cost minutes per run, but the extension's share of that is unknown.

## 14. Conclusions

**Demonstrated.** A working three-dimensional level-set mean-curvature-flow simulation in
EvolvingDomains, driven at every step by curvature computed from the evolving level set and not
by the analytic solution; quantitative agreement with R(t) = √(R₀² − 4μt) at the 0.15–0.24 %
level in radius and 0.41–0.65 % in volume after t = 0.6, with both errors decreasing under
refinement; near-spherical evolution (σ_R at or below 10⁻⁴, not growing); explicit, readable
mathematical orchestration in five small operators and one loop; public `CutSurface` geometry
consumed by an independent Curvature.jl fit with zero failures and mean κ agreeing to 2 × 10⁻⁴;
and reinitialisation every ten steps that leaves the measured interface unchanged to round-off
while restoring |∇φ| to ≈ 0.997–0.998.

**Not demonstrated.** Curvature.jl directly driving the moving interface — the facet-to-node
extension was neither found nor implemented; a complete and efficient facet→grid extension
workflow, including any normal-extension step; rigorous convergence of the evolved solution (two
grids, a two-point trend, and a residual whose order is unresolvable); stability or accuracy for
strongly nonspherical interfaces; and any diagnosis of the non-finite-velocity episode seen when
many short-lived fields were cut in one process, which remains a flagged reproducibility risk
with several plausible and undistinguished causes.

**Assessment.** On the scientific question, the evidence is positive but bounded: the scheme
evolves a genuine three-dimensional interface under numerically computed mean curvature, tracks
the analytic sphere within a fraction of a percent over two grid spacings, and preserves
sphericity to the resolution of the estimator. On the architectural question, the evidence is
positive for *data-level* composability — a second library consumed the package's geometry
through its public surface and returned a mathematically independent estimate — and explicitly
incomplete for *operator-level* composition, because the interface quantity produced by that
library cannot yet reach the transport without an adapter that this campaign did not find,
implement, or validate.

That separation — between what was measured, what was computed by which component, and what
remains unknown — is the substantive result of this first campaign.
