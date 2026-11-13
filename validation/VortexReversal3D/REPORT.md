# VortexReversal3D — Scientific and Architectural Report

**Benchmark:** three-dimensional reversible vortex transport of a level set and a passive
scalar on `EvolvingDomains.jl`
**Location:** `validation/VortexReversal3D/`
**Sources:** `run.VortexReversal3D.jl` (397 lines), `visualize.VortexReversal3D.jl` (233 lines)
**Executed command:**

```bash
julia --project="validation/env" \
  "/Users/naoufel/Documents/Scientific stack/EvolvingDomains-Redesign/validation/VortexReversal3D/run.VortexReversal3D.jl"
```

Julia 1.11.7; `validation/env` provides `EvolvingDomains` (path dependency), `StaticArrays`,
`CairoMakie` 0.15.17. No package source, project file, or manifest was modified.

---

## 3.1 Executive summary

A genuine three-dimensional reversible-vortex benchmark was implemented and executed at
$40^3$ and $64^3$ nodes over the full interval $0 \le t \le 1$. The prescribed velocity is
$\mathbf v(\mathbf x,t)=\cos(\pi t)\,\mathbf V(\mathbf x)$ with $\nabla\cdot\mathbf V=0$ and
$\mathbf V\cdot\mathbf n=0$ on $\partial D$, so the exact fields return to their initial
values at $t=1$ and the return error is an analytic reference requiring no trajectory
integration.

Both transport methods were exercised through public `EvolvingDomains` interfaces:
`eulerianTransport!` (WENO5 + SSP-RK3) on the level set and on one copy of the scalar, and
`lagrangianTransport!` (directionally split CIP) on a second copy of the scalar. The two
schemes consumed the same `Grid`/`Discretization`/`Field` representation and the same frozen
velocity array, with the dynamics, time loop and diagnostics owned entirely by the
application.

Level-set and scalar return errors decreased with refinement. CIP produced lower scalar
return errors than WENO5 at both resolutions, in both $L^1$ and $L^2$. Three independent
volume estimators showed small relative drift over the reversal ($\le 5.5\times10^{-4}$ at
$64^3$ for the cut-cell estimate).

Two qualifications are important. First, the chosen flow admits a first integral whose
level sets are approximately spherical near the box centre, so a centred sphere is a
comparatively mild test of geometric deformation. Second, CIP exhibited substantial
per-step allocations despite a supplied `CIPCache`, dominating both total allocation and
garbage-collection time.

No upstream source modifications were required. With only two resolutions, no asymptotic
convergence order is claimed.

---

## 3.2 Mathematical problem

Domain and horizon:

$$D=[0,1]^3,\qquad 0\le t\le 1 .$$

Prescribed velocity:

$$
\mathbf v(x,y,z,t)=\cos(\pi t)\,
\begin{pmatrix}
2\sin^2(\pi x)\sin(2\pi y)\sin(2\pi z)\\
-\sin(2\pi x)\sin^2(\pi y)\sin(2\pi z)\\
-\sin(2\pi x)\sin(2\pi y)\sin^2(\pi z)
\end{pmatrix}
=\cos(\pi t)\,\mathbf V(\mathbf x).
$$

Both fields satisfy the same advection equation

$$\partial_t u+\mathbf v\cdot\nabla u=0,\qquad u\in\{\phi,c\},$$

with initial data

$$\phi_0(\mathbf x)=\|\mathbf x-(0.5,0.5,0.5)\|-0.16,$$

$$c_0(\mathbf x)=\exp\!\left(-\frac{\|\mathbf x-(0.43,0.52,0.49)\|^2}{0.08^2}\right).$$

**Divergence-free.** Writing $s_j=\sin(\pi x_j)$ and using $s_{2j}:=\sin(2\pi x_j)$, all
three divergence terms are the same product $S=s_{2x}s_{2y}s_{2z}$ times a scalar:

$$\partial_x v_x=2\pi\cos(\pi t)\,S,\qquad
\partial_y v_y=-\pi\cos(\pi t)\,S,\qquad
\partial_z v_z=-\pi\cos(\pi t)\,S,$$

so $\nabla\cdot\mathbf v=(2\pi-\pi-\pi)\cos(\pi t)S=0$.

**Tangent to the boundary.** Every component carries a factor $\sin^2(\pi x_j)$ or
$\sin(2\pi x_j)$ for each coordinate, and both vanish at $x_j\in\{0,1\}$; hence
$\mathbf v\equiv\mathbf 0$ on all six faces of $\partial D$, so no mass crosses the boundary.

**Reversible.** $\mathbf v(\mathbf x,1-t)=-\mathbf v(\mathbf x,t)$ for all $t$. If
$\Phi_t$ denotes the flow map, this implies $\Phi_{1-t}=\Phi_t^{-1}$, hence the exact
solution satisfies

$$u(\mathbf x,1)=u(\mathbf x,0)\quad\text{for every initial field }u.$$

The return error therefore measures the total transport error of the numerical scheme
without a manufactured solution, without trajectory bookkeeping, and without a reference
solver. It is a strict test: it aggregates every failure mode along the whole path,
including numerical diffusion, dispersion and phase error.

**Reported error norms.** For a final field $u_N$ and its initial state $u_0$,

$$\varepsilon_{L^1}=\frac{\sum|u_N-u_0|}{\sum|u_0|},\qquad
\varepsilon_{L^2}=\frac{\big(\sum|u_N-u_0|^2\big)^{1/2}}{\big(\sum|u_0|^2\big)^{1/2}},$$

sums taken over the nodal values on the grid (see `relativeErrors` in
`run.VortexReversal3D.jl`).

---

## 3.3 A subtle limitation of the benchmark geometry

The spatial factor admits a first integral. Define

$$I(x,y,z)=\sin^2(\pi x)\sin^2(\pi y)\sin^2(\pi z).$$

A direct computation gives

$$\mathbf V\cdot\nabla I = 2\pi\cos(\pi t)\bigl(8-4-4\bigr)
s_x^3c_x\,s_y^3c_y\,s_z^3c_z = 0 ,$$

so $I$ is constant along every trajectory: the flow maps each level set $I=\text{const}$
onto itself. This is a mathematical property of the prescribed velocity, independent of
the numerical scheme.

Near the centre, with $\mathbf x=(1/2,1/2,1/2)+\boldsymbol\xi$ and
$\sin\bigl(\pi(1/2+\xi_j)\bigr)=\cos(\pi\xi_j)$,

$$\log I=\sum_{j=1}^{3}2\log\cos(\pi\xi_j)=-\pi^2\|\boldsymbol\xi\|^2+O(\|\boldsymbol\xi\|^4),$$

using $\log\cos u=-u^2/2-u^4/12-\cdots$. The invariant surfaces are therefore
**approximately spheres** near the centre, with the anisotropy appearing only at fourth
order.

**Implication.** The initial interface is the sphere $\|\boldsymbol\xi\|=0.16$, centred
exactly on the point about which the invariant surfaces are spherical. To leading order in
$\|\boldsymbol\xi\|$, the initial interface *is* an invariant surface. A flow that maps an
interface onto itself can only rearrange material tangentially along it; radial deformation
is driven solely by the $O(\|\boldsymbol\xi\|^4)$ mismatch, which for
$\|\boldsymbol\xi\|=0.16$ is of relative size $\pi^4\|\boldsymbol\xi\|^4/6\approx1\%$ in
$\log I$.

**Conclusion (as interpretation, not as a solver defect).** This is a meaningful
three-dimensional transport test — a full 3D flow, a 3D interface, genuine tangential
advection across the whole box, and a strict return-error criterion — but it is
*comparatively mild* as a test of strong geometric deformation. The near-circular zero-level
contours visible in all three orthogonal slices at $t=0.5$ (section 3.10,
`output/n64/slices.png`) are consistent with this analysis rather than indicative of a
numerical failure to deform the geometry: at $t=T/2$ the geometry *does* change, but only
modestly (the level-set minimum moves from $-0.1463$ to $-0.1459$ at $64^3$, and the
cut-cell volume dips by $2.7\times10^{-4}$ at $64^3$).

Independent supporting observation: the passive scalar, a compact Gaussian inside the same
region, *does* visibly filament at $t=0.5$ (`output/n64/scalar_schemes.png`). This is
consistent — a smooth scalar field responds to the local shear, whereas a set that nearly
coincides with an invariant surface cannot be strongly stretched by a flow preserving that
surface.

**Recommendation:** an off-centre sphere, or an ellipsoid, whose support avoids the
approximately spherical invariant surfaces, would exercise strong deformation. This was not
implemented here.

---

## 3.4 Software design and actual implementation

### Data model

- `Grid` owns the three uniform node axes and a mesh-identity token; `Discretization` maps
  the stable level symbol `:init` to it; `Field` attaches one dense `Array{Float64,3}` per
  level. Four `Field`s (`geometry`, `scalarEuler`, `scalarCip`, `scalarCipTarget`) share one
  `Discretization`.
- All numerical state is reachable as plain arrays through `field.data[:init]`. Diagnostics,
  snapshots (`copy`), comparisons with the initial data, binary persistence and plotting all
  operate on those arrays directly; nothing in this project needed a solver-specific wrapper
  type.

### The loop actually executed

```julia
for step in 1:config.steps
    v = dynamics(config, config.time[] + config.dt / 2)

    advanceGeometry!(config, v)
    advanceState!(config, :weno5, v)
    advanceState!(config, :cip, v)

    config.time[] = step * config.dt
end
```

`dynamics` fills a pre-allocated `Vector{SVector{3,Float64}}`. `advanceGeometry!` calls
`eulerianTransport!` on the level set. `advanceState!` dispatches on a `Symbol`: `:weno5`
calls `eulerianTransport!` on `scalarEuler`, `:cip` calls
`lagrangianTransport!(target, source, velocity, dt)` on the CIP pair. Snapshots are pushed at
step `steps÷2` and at `steps`, with `steps` forced even so that $t=T/2$ is a step boundary.

### Division of responsibility

| Owned by EvolvingDomains | Owned by the application |
|---|---|
| `Grid` / `Discretization` / `Field` container semantics | the velocity field and its time dependence |
| WENO5 + SSP-RK3 Eulerian update of one level | the global time loop, step size, midpoint evaluation |
| directionally split CIP Lagrangian update between two Fields | buffer management for CIP (distinct source/target) |
| `WENO5Cache`, `CIPCache` workspace types | level-set versus scalar bookkeeping |
| `cut` snapshot: background states, interface facets, normals | all diagnostics (return errors, volumes, extrema) |
| argument validation (aliasing, dimensions, finiteness, uniform spacing) | persistence and visualization |

### Eulerian versus Lagrangian

`eulerianTransport!` is a fixed-stencil, method-of-lines update: WENO5 reconstruction of
$-\mathbf v\cdot\nabla u$ with SSP-RK3, using a frozen nodal velocity, and it updates the
`Field` in place. `lagrangianTransport!` is a departure-point method: at each node it traces
the characteristic backwards by $-\Delta t\,\mathbf v$, builds a cubic Hermite profile from
nodal values and one-sided derivatives, and samples it there, done in directional splits.
Consequently its interface is `(target, source, velocity, dt)` — two Fields — rather than
in-place, and the package raises `ArgumentError` if the two alias. The run therefore keeps
`scalarCip` and `scalarCipTarget` distinct and copies the target back after each step
(`copyto!`), which preserves the documented non-aliasing contract without ping-pong state in
the configuration.

Because $\mathbf V$ vanishes on $\partial D$, no departure point can leave the box, so CIP's
documented refusal of outer-boundary inflow was never triggered over 512 total steps.

### Verified velocity ordering

Both kernels index the velocity vector directly with a linear index whose **first axis is
fastest**, a convention that no type enforces. The run does not assume it: `verifyVelocityOrdering`
transports the affine field $1+2x-\tfrac12y+3z$ with a constant diagonal velocity and
compares against the exact translate on the interior. Measured maximum errors were
$1.78\times10^{-15}$ at $40^3$ and $3.55\times10^{-15}$ at $64^3$; a transposed ordering
would produce $O(1)$ errors. The check is cheap and is repeated per resolution.

### Temporal accuracy

The velocity is evaluated **once per global step at the midpoint** $t+\Delta t/2$ and shared
by the level set and both scalars. Freezing at the midpoint makes the advection term second
order in $\Delta t$ for smooth $\mathbf v$; the overall temporal order is bounded by that
freezing and by SSP-RK3 (level set) and the operator splitting (CIP), and is not separated
from the spatial error in these measurements. No temporal refinement study was performed, so
the achieved temporal order is not claimed — only that the scheme is second-order-consistent
in time.

---

## 3.5 Numerical implementation and stability

Uniform node lattices on $[0,1]^3$; `h` is the minimum node spacing.

| Quantity | $40^3$ | $64^3$ |
|---|---:|---:|
| Number of time steps | 196 | 316 |
| Time step | $5.1020\times10^{-3}$ | $3.1646\times10^{-3}$ |
| Nominal CFL ($\Delta t\,\lvert\mathbf v\rvert_\infty/h$) | 0.4 | 0.4 |
| $h$ | $2.5641\times10^{-2}$ | $1.5873\times10^{-2}$ |
| $\lvert\mathbf V\rvert_\infty$ (sampled, $t=0$) | 1.99440 | 1.99784 |

**CFL estimate.** Since $|\cos(\pi t)|\le1$, the instantaneous speed is bounded by
$|\mathbf V|_\infty$, and the sampled maximum at $t=0$ is the worst-case bound. The step is

$$\Delta t_{\text{CFL}}=\frac{0.4\,h}{\max_k|\mathbf V(\mathbf x_k)|},$$

then rounded down to an even step count:

$$\texttt{steps}=2\left\lceil \frac{T}{2\,\Delta t_{\text{CFL}}}\right\rceil,\qquad
\Delta t = \frac{T}{\texttt{steps}},$$

so $T=1$ is reached exactly by construction and $t=T/2$ falls on a step boundary. Time is
stored as `step*dt` rather than accumulated (the recorded snapshot times are
$0$, $0.49999999999999994$, $0.9999999999999999$ — within one ulp of $T/2$ and $T$). The
step is stability-appropriate rather than merely non-divergent: the package documents that
WENO5 + SSP-RK3 stays bounded at CFL near 1 and collapses near 3, and CIP additionally
rejects any departure point outside the domain, which was never raised.

**No reinitialization.** `reinitialize!` is not called anywhere in the run and $\phi$ is
never reset from the initial sphere; the level set is only ever advanced by
`eulerianTransport!`. This is deliberate: reinitialization projects $\phi$ back towards a
signed-distance function, which would destroy the very thing the benchmark measures — the
accumulated transport error of the level-set values. It would also make the return error
non-comparable between schemes and between resolutions, because the reinitialization
operator introduces its own error and its own grid dependence. Reporting it "separately"
would require a second, differently defined experiment.

---

## 3.6 Numerical errors and grid refinement

Recorded values (`output/n40/errors.csv`, `output/n64/errors.csv`):

| Quantity | $40^3$ | $64^3$ | reduction factor |
|---|---:|---:|---:|
| Level-set relative $L^1$ return error | $4.00\times10^{-4}$ | $7.48\times10^{-5}$ | 5.35 |
| Level-set relative $L^2$ return error | $1.07\times10^{-3}$ | $2.61\times10^{-4}$ | 4.11 |
| WENO5 scalar relative $L^1$ return error | $8.74\times10^{-2}$ | $1.94\times10^{-2}$ | 4.49 |
| WENO5 scalar relative $L^2$ return error | $1.04\times10^{-1}$ | $2.75\times10^{-2}$ | 3.76 |
| CIP scalar relative $L^1$ return error | $3.04\times10^{-2}$ | $1.54\times10^{-2}$ | 1.97 |
| CIP scalar relative $L^2$ return error | $2.67\times10^{-2}$ | $1.32\times10^{-2}$ | 2.03 |

All errors are measured against the exact initial state after the full reversal
($u_N=u(\cdot,1)$, $u_0=u(\cdot,0)$), i.e. they are pure return errors.

**Verified observations.** Every quantity decreases under refinement; the level set returns
with the smallest error by more than two orders of magnitude; CIP is more accurate than
WENO5 for the scalar at both resolutions and in both norms; WENO5's scalar error falls by a
factor of 3.8–4.5 while CIP's falls by a factor of about 2.

**What is *not* claimed.** Two resolutions cannot establish an asymptotic order — a two-grid
estimate would be indistinguishable from pre-asymptotic behaviour. The reduction factors are
reported as measurements, not converted into orders.

**Plausible explanation (hypothesis).** The different trends are consistent with the initial
Gaussian being under-resolved, and with the two schemes degrading differently on a
sharp-gradient profile. With $\sigma_g=0.08$ and $h=1/39$,

$$\sigma_g/h\approx3.12 .$$

If one uses the conventional Gaussian standard deviation from
$\exp(-r^2/(2\sigma^2))$, then $\sigma=\sigma_g/\sqrt2\approx0.0566$ and
$\sigma/h\approx2.21$ at $40^3$ ($\approx3.56$ at $64^3$). The bump is therefore resolved by
only a few cells — correctly described as *narrow*, not as one grid spacing wide. A narrow
initial profile mixes initial-data resolution error into the transport error in a way that
depends on the scheme's dissipation, which could explain both the varying reduction factors
and the narrowing of the CIP advantage from a factor 2.9 ($L^1$) at $40^3$ to 1.26 at
$64^3$. This remains a hypothesis until a wider-Gaussian comparison is run.

**Other recorded extrema** (evidence that nothing diverged and that neither scheme is
positivity preserving):

| | $40^3$ | $64^3$ |
|---|---:|---:|
| $\min\phi$ at $t=0$ / $t=T/2$ / $t=T$ | $-0.1378$ / $-0.1368$ / $-0.1360$ | $-0.1463$ / $-0.1459$ / $-0.1455$ |
| $\max\phi$ at all times | $0.70603$ | $0.70603$ |
| scalar range at $t=T$, WENO5 | $[-1.8\times10^{-6},\,0.948]$ | $[-2.0\times10^{-7},\,0.989]$ |
| scalar range at $t=T$, CIP | $[-1.7\times10^{-4},\,0.982]$ | $[-9.9\times10^{-8},\,0.994]$ |

$\max\phi$ is exactly $\sqrt3/2-0.16$, the initial corner value, i.e. the level set retains
its outer scale. Both scalar schemes produce small negative undershoots ($O(10^{-4})$ at
$40^3$, $O(10^{-7})$ at $64^3$); neither is bound preserving. The initial scalar maximum is
$0.985$ ($40^3$) and $0.997$ ($64^3$) because the Gaussian's peak at $(0.43,0.52,0.49)$ is
not a grid node.

---

## 3.7 Geometry and volume preservation

Three independent estimators of the enclosed volume of $\{\phi<0\}$ were computed. The
analytic reference is

$$V_0=\frac{4\pi R^3}{3}=4\pi(0.16)^3/3=1.71573\times10^{-2}.$$

**1. Voxel counting** — $V=h^3\,\#\{\phi_k<0\}$ (function `voxelVolume`). A zeroth-order
indicator sum: first-order in $h$, biased by roughly half a cell layer of surface area, and
it returns only quantized values.

**2. Smoothed Heaviside** — $V=h^3\sum_k \operatorname{clamp}\bigl(\tfrac12-\phi_k/h,0,1\bigr)$
(function `smearedVolume`). The cell-scale smearing makes the nodal sum exact for a planar
interface and second-order for a curved one; it is still a nodal quadrature, not a geometric
reconstruction.

**3. Cut-cell integration** — the package's own `cut` snapshot is turned into a volume by the
divergence theorem over its interface facets (function `cutCellVolume`):

$$V=\frac13\oint_{\partial\{\phi<0\}}\mathbf x\cdot\mathbf n\,dA
=\sum_{f}\frac{(\mathbf a_f+\mathbf b_f+\mathbf e_f)\cdot\mathbf n_f}{18}
\bigl\lVert(\mathbf b_f-\mathbf a_f)\times(\mathbf e_f-\mathbf a_f)\bigr\rVert .$$

This is genuinely cut-cell based: the facets come from the package's six-tetrahedron
piecewise-linear 3D cut model, and the normals are outward from the $\phi<0$ phase.

Recorded estimates at $64^3$, and their error relative to $V_0$:

| Estimator | $t=0$ | error vs $V_0$ | $t=T$ | drift vs its own $t=0$ |
|---|---:|---:|---:|---:|
| Voxel | $1.70848\times10^{-2}$ | $-0.42\%$ | $1.70848\times10^{-2}$ | $0.0$ |
| Smoothed Heaviside | $1.72029\times10^{-2}$ | $+0.27\%$ | $1.71938\times10^{-2}$ | $-5.26\times10^{-4}$ |
| Cut-cell | $1.70725\times10^{-2}$ | $-0.49\%$ | $1.70631\times10^{-2}$ | $-5.51\times10^{-4}$ |

At $40^3$ the same estimators give $+4.55\%$ (voxel), $+0.91\%$ (smoothed) and
$-1.27\%$ (cut-cell) at $t=0$.

**Two distinct quantities.** *Absolute volume accuracy* (how well the discrete interface
placement reproduces the analytic sphere at a given grid) and *volume drift over time* (how
much a fixed estimator changes as the field is transported) are different measurements. At
$64^3$ the smoothed-Heaviside estimate is the most accurate in absolute terms (+0.27%) while
the cut-cell estimate is the most accurate at $40^3$; all three drift by at most
$5.5\times10^{-4}$ over the reversal. The exact zero drift of the voxel count is not evidence
of good conservation: the indicator is insensitive to sub-cell motion, and it returned the
same integer count at $t=0$ and $t=T$ at both resolutions.

**Validity condition for the divergence estimate.** The closed-surface integral equals the
enclosed volume only if the negative phase does not touch $\partial D$; otherwise the box
walls contribute a surface the snapshot does not represent. The run therefore counts
background cells on $\partial D$ whose state is not `CUT_OUT` (function `wallContacts`).
**Verified:** that count was 0 at all three snapshot times at both resolutions (column
`wall_contacts` in `metrics.csv`), so the estimate is valid as an enclosed volume here.

**Limitation.** Neither whole-domain level-set norms nor global volume characterize interface
*geometry*: a field can have small $L^1$ error and a locally displaced surface, and volume can
be preserved while the interface is distorted (or vice versa). A direct measurement of the
displacement between the exact and reconstructed zero level sets is required and was not
performed.

---

## 3.8 Computational performance

Measured by `@timed` over the time loop only (kernel warm-up, diagnostics, I/O and plotting
excluded — for reference, the whole program, including both resolutions and all figures,
took 2 min 2 s of wall clock):

| Quantity | $40^3$ | $64^3$ |
|---|---:|---:|
| Wall-clock time (loop) | 13.1 s | 87.6 s |
| Allocations | 13.3 GB | 85.5 GB |
| Garbage-collection time | 3.1 s | 17.5 s |
| GC fraction of loop time | 24% | 20% |

Total allocation tracks nodes × steps: $6.4\times$ for $4.1\times$ the nodes and $1.61\times$
the steps, i.e. per-step allocation ($68$ MB at $40^3$, $271$ MB at $64^3$) is close to
proportional to the node count. No scaling study beyond these two points was performed.

**Where it comes from — measured, but not exhaustively attributed.** Cached WENO5 is
comparatively allocation-efficient: an ad-hoc `@allocated` measurement of a single
`eulerianTransport!` call with `WENO5Cache` at $40^3$ after compilation returned about
2.8 kB. The same measurement for a single `lagrangianTransport!` call with `CIPCache`
returned about 71 MB (~1.1 kB per node), and the copy-back `copyto!` returned 0 bytes. That
localizes the allocation to the CIP call. It does **not** identify the responsible internal
arrays — that measurement is outside the saved artifacts, no profiler was run, and the split
between derivative reconstruction, departure-point interpolation, per-axis temporaries and
cache handling is unknown.

**Framing.** This is a performance concern for the next campaign, not a defect claim and not
an invitation to modify the package: `CIPCache` exists precisely to make this path reusable,
and the observed per-step allocation suggests its workspace is not covering everything the
CIP sweep needs. The practical consequence measured here is real: about a fifth of the
runtime, and tens of gigabytes of allocation, are spent on memory management.

---

## 3.9 Public API and architectural evaluation

### Minimalism

The data model is three structs. `Grid` (axes + identity), `Discretization` (symbol →
`Grid`), `Field` (symbol → dense array). In practice this meant the entire application could
be written against `field.data[:init]`: snapshots are `copy`, error norms are `sum(abs, ...)`,
min/max are `minimum`/`maximum`, persistence is `write(io, array)`, and plotting reads the
same arrays back from disk. No abstraction, adapter or solver type was needed anywhere in
397 + 233 lines, and the orchestration reads as the mathematics:

$$\texttt{v = dynamics(t)}\,\to\,\texttt{advanceGeometry!}\,\to\,\texttt{advanceState!}.$$

The minimalism is genuine, with one consequence: because the data model does not carry a
grid for the velocity, the node-ordering convention is untyped and had to be verified
explicitly (section 3.4). That convention is load-bearing and silently wrong is a possible
failure mode — which is why the run asserts it.

### Composability

The same scalar initial data was transported by two different numerical methods supplied by
the package, with a shared velocity representation and no cross-talk:

- both methods accept the *same* `Vector{SVector{3,Float64}}` built once per step;
- both operate on `Field`s attached to the *same* `Grid`/`Discretization`;
- neither knows about the other, about the time loop, or about the diagnostics.

The Eulerian scheme is in-place and the Lagrangian one is a mapping; the application
absorbs that difference with a target buffer and a copy-back, i.e. the divergence between
the two interfaces is handled *outside* the package rather than by a unifying abstraction.
That is exactly the composability claim being tested, and it holds for methods inside the
package.

### What this does *not* demonstrate

It demonstrates interoperability between transport operators **shipped by EvolvingDomains
itself**. No third-party or independently developed transport discretization was integrated;
there is no evidence here about consuming these fields from an external PDE framework (for
example Gridap or Trixi), beyond the observation that the state is an ordinary `Array`.

### Remaining limitations, by kind

*Packaging / environment*
- `Artifacts.toml` is absent from this checkout, so `cut` cannot resolve the native CutCells
  library by design; the run sets the documented `EVOLVINGDOMAINS_CUTCELLS_LIBRARY` override
  to the prebuilt `validation/lib/cutcells/lib/` library. The fallback path — reporting the
  cut-cell volume as `missing` with the error text — exists in the code but was never
  exercised, because the library was always found.

*Performance*
- Substantial CIP per-step allocation despite `CIPCache` (section 3.8).

*Mathematical API limitations*
- The tested `Cut` interface exposes background-cell states, interface facets and facet
  normals, but **no subcell volume weights**, so a cut-cell volume required the
  divergence-theorem workaround over the boundary facets. Note also that
  `boundary.normals` is indexed per facet while `physical_points` is indexed per point.
- Velocity indexing conventions are implicit: node order per axis, and "first axis fastest",
  are documented but not enforced by any type.
- The transport operators carry no CFL or stability diagnostic; choosing `dt` is entirely the
  caller's responsibility, and the documented guidance had to be reproduced manually.

### Overall

The architecture claim holds for this benchmark: a minimal data model plus independent
operators let a user own the dynamics, the time loop, the diagnostics and the visualization
while reusing the package's representation and kernels. The evaluation is scoped — one
application, two operators from the same package, one flow, one geometry family — so it
supports "the design works as advertised for this kind of coupling" and not "the design is
proven interoperable in general".

---

## 3.10 Visualization results

All figures were produced by `visualize.VortexReversal3D.jl` from the saved arrays
(`params.txt`, `fields.bin`, `metrics.csv`, `errors.csv`); the visualization module contains
no transport, no physics and does not modify state. It can be re-run independently with
`julia --project=validation/env validation/VortexReversal3D/visualize.VortexReversal3D.jl`.

**Return error against resolution** — `output/convergence.png`
Level-set and scalar return errors for both schemes, $L^1$ and $L^2$, at $40^3$ and $64^3$.
It shows the monotone decrease with refinement and the CIP-over-WENO5 ordering. With two
points it establishes direction, not order, and the near-vertical appearance reflects the
narrow resolution range (40 → 64).

**Passive scalar, both schemes** — `output/n40/scalar_schemes.png`, `output/n64/scalar_schemes.png`
The $xy$ slice through the bump ($z\approx0.49$) at $t=0$, $T/2$, $T$, identical colour scale
per resolution. It shows the vortex bending the Gaussian into a tail at $T/2$ and both
schemes returning a compact, slightly diffused profile at $T$; CIP retains a visibly higher
peak and a tighter profile. It cannot separate transport error from initial under-resolution,
and it shows one plane only.

![Passive scalar, WENO5 vs CIP, 64³](output/n64/scalar_schemes.png)

**Orthogonal slices of the level set** — `output/n64/slices.png`
$xy$, $xz$, $yz$ slices through the box centre with the $\phi=0$ contour, at the three
snapshot times. The contours remain close to circular, with a small genuine distortion at
$T/2$ — the visual counterpart of the invariant-surface argument of section 3.3. It is a
three-plane sample, so it cannot certify an interface shape, and a slice through the centre
is the least informative plane for a nearly spherical geometry.

![Level-set slices at 64³](output/n64/slices.png)

**Projected 3D geometry** — `output/n64/geometry3d.png`
Two oblique orthographic projections of the zero-crossing point set (the points where $\phi$
changes sign along a mesh edge, linearly interpolated), at the three snapshot times, coloured
by $z$. This is **not** a rendered isosurface: CairoMakie 0.15.17 builds a `Volume` plot for a
3D `contour!`/`volume!` and the Cairo backend drops it with a warning, so neither an
isosurface nor a volume rendering is available. The projection conveys the global shape and
its mild deformation; it cannot establish surface area, curvature or interface displacement.

![Projected zero-level geometry at 64³](output/n64/geometry3d.png)

**Volume drift and return errors** — `output/n64/metrics.png`
Left: relative drift of the three volume estimators over the snapshots. Right: relative
$L^1$ and $L^2$ return errors per quantity on a log scale. It shows the cut-cell and
smoothed-Heaviside tracks nearly coinciding (they differ by about 5% in drift at $64^3$)
and the voxel track's larger, non-monotone excursion at $T/2$ ($-2.3\%$ at $40^3$, $+0.28\%$
at $64^3$). Three time samples only, so no claim about the drift's shape in time.

![Volume drift and return errors at 64³](output/n64/metrics.png)

All figures were checked against the recorded CSV values; no number in this report is taken
from a figure.

---

## 3.11 Recommendations for the second campaign

None of these were executed.

**Experiment A — Stronger geometric deformation.**
Move the initial sphere off the approximately spherical invariant surfaces, or replace it
with an anisotropic ellipsoid, and/or use a flow without a spherical near-centre first
integral. Measure deformation at $t=0.5$, return error at $t=1$, and interface displacement.
*Question answered:* whether the package's transport operators maintain accuracy under strong
stretching and thinning, which the current geometry does not test. It also separates
"scheme is accurate" from "geometry is easy" — the present $O(10^{-4})$ level-set error may
be flattered by the near-invariance of the centred sphere.

**Experiment B — CIP allocation profiling.**
Independently time and measure allocations for one CIP and one WENO5 transport step after
compilation, and profile inside CIP to attribute the bytes to derivative reconstruction,
departure-point interpolation, per-axis temporaries, or cache handling.
*Question answered:* whether `CIPCache` covers the workspace the sweep needs — a performance
question with a concrete architectural consequence (whether the cache contract is complete).

**Experiment C — Improved convergence study.**
At least three resolutions, a wider initial Gaussian, and, where practical, spatial/temporal
separation by refining $\Delta t$ at fixed $h$ and vice versa.
*Question answered:* the actual asymptotic order of WENO5 and CIP in this setting, and
whether the observed CIP/WENO5 gap narrowing is an initial-data artefact.

**Experiment D — Interface geometry error.**
Measure the discrepancy between exact and reconstructed interface directly — a signed
distance sampled on the interface, or a symmetric-difference measure of the enclosed sets —
rather than whole-domain level-set norms.
*Question answered:* whether small global norms actually imply a well-placed interface, which
volume drift alone cannot establish.

**Experiment E — Independent external transport.**
Implement or integrate a transport discretization outside EvolvingDomains that consumes the
same exposed arrays, and compare its output with WENO5 and CIP.
*Question answered:* whether the composability claim extends past the package's own operators
— the strongest architectural question raised by this benchmark, and one that this campaign
deliberately leaves open.

---

## 3.12 Conclusions

**Demonstrated.**
- Three-dimensional operation over the full interval $0\le t\le1$ on genuine $40^3$ and
  $64^3$ lattices, with all saved values finite and every diagnostic computed from the data.
- Both public transport interfaces — `eulerianTransport!` (WENO5 + SSP-RK3) and
  `lagrangianTransport!` (CIP) — exercised over 512 total steps, sharing one velocity array
  and one data representation, with an explicitly verified node ordering.
- An explicit, mathematically readable simulation loop in which the dynamics, time loop,
  diagnostics and visualization remain application code.
- Return errors decreasing with refinement for the level set and both scalar schemes; CIP
  more accurate than WENO5 for the scalar at both resolutions.
- Independent geometric diagnostics: three volume estimators, a cut-cell integration from the
  package's own interface facets, and a verified proof that the enclosed phase never touches
  the computational boundary.

**Not demonstrated.**
- Strong geometric deformation under a demanding flow — the centred sphere is close to an
  invariant surface of the prescribed velocity.
- A rigorous asymptotic convergence order, or a separated temporal order (two resolutions,
  one time step per resolution).
- Allocation-efficient CIP execution, or a full attribution of its allocation.
- Interoperability with independently developed external PDE packages.

**Assessment.** The benchmark is a meaningful first validation of the minimal-and-composable
design: the package behaved as an orchestrator of shared state rather than a monolithic
solver, and no source modification was needed to build a complete three-dimensional
reversible-transport experiment around it. The same campaign exposes specific, actionable
questions — a geometry that under-tests deformation, a CIP path that allocates heavily, an
interface-geometry diagnostic that is not yet measured, and interoperability that is only
demonstrated within the package's own operators. Those are the targets of the second
campaign.

---

### Appendix — saved artifacts

| Path | Content |
|---|---|
| `output/n40/`, `output/n64/` | one run each: `params.txt`, `metrics.csv`, `errors.csv`, `metrics.txt`, `fields.bin`, four figures |
| `output/n{40,64}/fields.bin` | `Float64` records of `nodes³`, column-major (first axis fastest), in this order: `phi0`, `c0`, `phi` at $t=0,T/2,T$, `c_weno5` at $t=0,T/2,T$, `c_cip` at $t=0,T/2,T$ |
| `output/convergence.png` | return error versus resolution, both schemes, both norms |

Re-running `run.VortexReversal3D.jl` regenerates everything deterministically from the same
parameters; the visualization module can be run standalone against `output/`.
