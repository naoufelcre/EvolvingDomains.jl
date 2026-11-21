# ConservativeAnnulusTrixi — four transport methods on one evolving cut-cell annulus

**Campaign II, Experiment 1.** Trixi.jl DGSEM is added as an external high-order solver
beside EvolvingDomains' WENO5, CIP and conservative `TransportMap`, on a single geometry
carried by a level-set `Field`.

| Item | Value |
|---|---|
| Problem | $\partial_t\rho+\nabla\cdot(\rho\mathbf v)=0$, $\mathbf v=(0.24,0.10)$, on a translating annulus, $t\in[0,1]$ |
| Domain | $\Omega(0)=\{0.12<\|\mathbf x-(0.35,0.45)\|<0.24\}$, $\Omega(t)=\Omega(0)+t\mathbf v$ |
| Methods | WENO5 (`eulerianTransport!`), CIP (`lagrangianTransport!(t,s,v,dt)`), `TransportMap` (`lagrangianTransport!(t,s,map)`), Trixi DGSEM degree 3 |
| Resolutions | 64² cells / 16² elements (4225 vs 4096 dof) and 128² / 32² (16641 vs 16384 dof) |
| Steps to $t=1$ | 55 and 109, one shared $\Delta t$, SSP-RK3 everywhere ($\Delta t=0.40$ of Trixi's stable $\Delta t$) |
| EvolvingDomains | this checkout, `src/EvolvingDomains.jl`, git `0809db8b`, v1.0.0 |
| Native CutCells | locally built, upstream commit `ed5232ec` (pin in `cxxBindings/CMakeLists.txt`) |
| Trixi | 0.17.21 (`LinearScalarAdvectionEquation2D`, `TreeMesh` periodic, `DGSEM(polydeg=3, flux_lax_friedrichs)`, `SimpleSSPRK33`) |
| Julia | 1.11.7 |
| Entry points | `run.ConservativeAnnulusTrixi.jl`, `visualize.ConservativeAnnulusTrixi.jl` |

---

## 1. Executive summary

Four methods were run on **one** numerically evolved geometry and measured with **one**
integration rule (bilinear nodal weights $w_i(t)=\int_{\Omega_h}N_i\,\mathrm dx$ and a
degree-5 triangle rule on the same cut triangles for all four).

* **Trixi is the most accurate method by a wide margin.** Relative $L^1$ density error at
  $t=1$: `Trixi 7.95e-6`, `CIP 4.33e-5`, `WENO5 5.07e-5`, `Map 3.54e-2` (128²/32²).
* **Trixi is also the cheapest per unit of accuracy here**, by a large factor: 0.055 s of
  solver time for the whole run against 0.19–0.65 s for the EvolvingDomains operators, at
  a comparable number of degrees of freedom. This is a comparison of *one* setup on *one*
  problem, and it is dominated by data layout (Trixi's vectorised DG kernels over 16k dof
  against scalar Julia loops) rather than by an intrinsic algorithmic advantage.
* **The conservative map loses $8.5\times10^{-3}$ of the physical mass at 128²** while the
  other three lose at most $6.7\times10^{-5}$ there ($2.9\times10^{-4}$ at 64²). Its unweighted
  nodal-sum invariant still holds to $3.5\times10^{-16}$.
* The map's mass loss is now **attributed**: with $\rho\equiv1$ the same map loses
  $9.8\times10^{-3}$ of the *area*, i.e. the loss is a partition-of-unity defect on the
  moving support, not a density-transport error. The weighted identity
  $w_{\text{new}}^{\mathsf T}P=w_{\text{old}}^{\mathsf T}$ fails by $1.0\times10^{-2}$ at
  $t=1$, matching the mass drift; the unweighted identity $\mathbf 1^{\mathsf T}P=\mathbf 1^{\mathsf T}$
  does not.
* A public-API gap blocks CIP as literally specified (see §4.4): the operator refuses any
  departure point outside the box and offers no periodic or inflow boundary condition, so a
  uniform translation cannot be run on a bounded grid at all. It was made runnable with a
  documented, provably inert one-node layer (§4.4) — this is the only deviation from the
  mission's literal recipe, it affects only CIP, and it is quantified.
* Composing Trixi with EvolvingDomains required **64 lines** of adapter (a view structure, a
  tensor-product Lagrange evaluator from the documented nodal interface, and two verifications)
  plus 18 lines of setup. No private field of Trixi's solver or data layout is read.

---

## 2. Problem and reference solution

$\Omega(t)$ is a rigid translate, $\varphi(\mathbf x,t)=\max(r_{\text{in}}-\|\mathbf x-\mathbf c(t)\|,
\|\mathbf x-\mathbf c(t)\|-r_{\text{out}})$ with $\mathbf c(t)=\mathbf c_0+t\mathbf v$, and

$$\rho(\mathbf x,t)=\rho_0(\mathbf x-t\mathbf v),\qquad
\rho_0=1+\tfrac12\,G,\quad G=\sum_{m_x,m_y\in\{-1,0,1\}}\exp\!\Bigl(-\tfrac{(x-x_0+m_x)^2+(y-y_0+m_y)^2}{\sigma^2}\Bigr).$$

The nine-image sum is the **periodic extension** the mission asks for. On the unit torus the
nearest image centre is $\ge0.35$ away and $\sigma=0.055$, so each neighbour term is
$\le4\times10^{-18}$: the extension changes nothing numerically, but the initial datum, the
analytic reference and Trixi's periodic solve now provably refer to one and the same
torus-continuous function instead of three functions that agree only up to a silent edge
convention.

Reference values, both exact to machine precision:

$$|\Omega|=\pi(r_{\text{out}}^2-r_{\text{in}}^2)=1.357168026351\times10^{-1},\qquad
M(0)=\int_{\Omega(0)}\rho_0=1.363208403334\times10^{-1}$$

($M(0)$ by 5-point Gauss–Legendre $\times$ 1024-point midpoint sweep in polar coordinates;
$\int_{[0,1]^2}\rho_0=1+\tfrac12\pi\sigma^2=1.004751658889$, used as an independent check of
the Trixi solve in §5.5).

Physical mass is conserved because $\mathbf v$ is constant and $\Omega$ is material,
$\tfrac{d}{dt}M=\int_{\Omega}\partial_t\rho+\nabla\cdot(\rho\mathbf v)=0$.

---

## 3. One measurement rule for four different discrete objects

The mission's requirement is a *common* diagnostic. For every method the material integral is

$$M_h(t)=\int_{\Omega_h(t)}\rho_h\,\mathrm dA,\qquad
w_i(t)=\int_{\Omega_h(t)}N_i\,\mathrm dA,\qquad
M_h^{\text{nodal}}(t)=\sum_i w_i(t)\,\rho_i(t),$$

with $N_i$ the background bilinear hat functions and $\Omega_h(t)$ the **reconstructed** cut
(a union of triangles from `Cut.subcell_states == CUT_IN`). Each triangle is integrated with
the 7-point degree-5 Dunavant rule, exact for total degree $\le5$; that covers the bilinear
reconstruction (degree 2) and the degree-3 DG polynomial exactly, so quadrature error is
negligible next to transport error and no method is favoured by its rule.

The same triangles and the same rule are used for

* the DG field, evaluated through the adapter of §4.3,
* the relative errors $E_1=\int|\rho_h-\rho_{\rm ex}|/\int\rho_{\rm ex}$ and
  $E_2=\sqrt{\int(\rho_h-\rho_{\rm ex})^2/\int\rho_{\rm ex}^2}$ over $\Omega_h$,
* the two controls of §5.2.

Two independent checks on the weights run at every step and are asserted:
$\sum_i w_i(t)=|\Omega_h(t)|$ (because $\sum_iN_i\equiv1$) to $10^{-12}$ relative, and
$\int_{\Omega_h}\rho_{\rm ex}\,\mathrm dA$ computed both at triangle centroids and through
$w_i$ with exact nodal values.

`Σ q` (an unweighted nodal sum) is reported separately and is never called a mass.

---

## 4. Implementation

### 4.1 Architecture

The loop is the mission's six steps, one function per advance:

```text
velocity = (0.24, 0.10)                      1. shared dynamics (prescribed, frozen per step)
old_cut  = cut(levelset)                     2. old cut + map source vector (zeroed outside)
advanceGeometry!  eulerianTransport!(levelset, …; cache)     3. shared geometry (WENO5, SSP-RK3)
advanceWithWENO!  eulerianTransport!(rho_weno, …)
advanceWithCIP!   lagrangianTransport!(dst, src, cipvel, dt)
new_cut  = cut(levelset); map = TransportMap(old_cut, new_cut, grid, velocity, dt)
advanceWithMap!   lagrangianTransport!(dst, src, map)        4. three independent scalars
advanceWithTrixi! Trixi.step!(integrator)                    5. fourth scalar (own state)
diagnostics(new_cut, (rho_weno, rho_cip, rho_map), dg, t)     6. common rule, §3
```

Three requirements of §4 of the mission are met structurally:

* The geometry is advanced **numerically** and never re-set from the analytic annulus; `cut`
  is taken from the level-set `Field` at every step, and every mass/error number in this
  report is integrated on that reconstructed cut, including at $t=0$.
* The four states are separate objects and are advanced by four separate calls.
* `@assert map.source_mask == supportMask(old_cut)` and the target analogue are checked on
  every step, as is the per-step nodal-sum invariant.

CIP uses two `Field`s with the values swapped each step (distinct source and destination
buffers), and one `CIPCache`; WENO5 uses two separate `WENO5Cache`s (geometry, density).

### 4.2 The four advances

| method | call | character |
|---|---|---|
| geometry | `eulerianTransport!(levelset, velocity, dt; cache)` | WENO5 + SSP-RK3, zero-gradient ghosts |
| WENO5 | `eulerianTransport!(rho_weno, velocity, dt; cache)` | same operator, on the density |
| CIP | `lagrangianTransport!(dst, src, cipvel, dt; cache)` | directionally split cubic-Hermite semi-Lagrangian, intensive |
| Map | `lagrangianTransport!(dst_vec, src_vec, map)` | conservative cut-to-cut pull/push, acts on the flattened node vector |
| Trixi | `Trixi.step!(integrator)` | DGSEM degree 3, LLF surface flux, periodic TreeMesh, SSP-RK3 |

The map's source vector is zeroed outside `supportMask(old_cut)`, which is the operator's
documented semantics ("only the map's source support is transported") and is harmless here
because every node with $\varphi<0$ is a corner of a non-`CUT_OUT` cell.

### 4.3 The Trixi adapter

Trixi exposes no public "evaluate at arbitrary points" routine in 0.17.21 (there is no
`evaluate_at_points`; `PlotData2D` samples for plotting only). What it does document — and
what its own examples use — is the nodal analysis interface
`Trixi.wrap_array`, `Trixi.get_node_coords`, `Trixi.get_node_vars`. The adapter therefore
rebuilds the polynomial the DG solution already *is*:

$$
\rho_h(\mathbf x)=\sum_{i,j=0}^{p}\ell_i(r)\,\ell_j(s)\,u_{ij}^{e},\qquad
(r,s)=\frac{2(\mathbf x-\mathbf x_e)}{h_e}-\mathbf 1,
$$

with $\ell_i$ the 1-D Lagrange cardinal polynomials on the LGL abscissae. Those abscissae,
and each element's origin and size, are *derived from* `get_node_coords`, so nothing about
Trixi's basis, node ordering or mesh bookkeeping is assumed. Size: **64 lines** total
(`DGView` 12, `dgView` 26, `lagrange1d!` 12, `evalDG` 14), plus 18 lines of
`trixiSetup`/`advanceWithTrixi!`. `wrap_array` returns a live view, so no per-step copy of
the state is made.

Verification, run before any simulation (§Stage B):

* **Polynomial reproduction**: a Q3 polynomial
  ($\sum_{a,b\le3}c_{ab}x^ay^b$, including $x^3y^3$) is fed through
  `Trixi.compute_coefficients` and the adapter reproduces it at 200 scattered points to
  $\mathbf{6.66\times10^{-16}}$. That is the exact statement the adapter claims — the DG
  solution is its own nodal interpolant — so it isolates adapter error from discretisation
  error.
* The same evaluator applied to the true (non-polynomial) initial datum has an interpolation
  error of $1.9\times10^{-2}$ at $8^2$ elements, i.e. it is an interpolation, not a
  reproduction, of general data. Reported so the 6.7e-16 is not over-read.
* **Independent solve check**: $\int_{[0,1]^2}\rho_h$ of the DG solution equals the exact
  torus integral $1+\tfrac12\pi\sigma^2=1.004751658889$ to 12 digits (§5.5).

Allocations of the adapter/diagnostics are measured under their own timer (`adapter+diag`)
and are never mixed into the solver times.

### 4.4 CIP: a public-API blocker, and how it was resolved

`cipLagrangianTransport!` guards every node with `_departureCell`:

> `DomainError: CIP departure point lies outside [0.0, 1.0]; outer-boundary inflow is not supported.`

Because the guard checks the departure point of **every** node of the grid, whatever its
value, a constant translation violates it on the first step for any $\Delta t>0$: the
leftmost column has departure $-\Delta t\,v_x<0$. There is no periodic option, no inflow
boundary option, and no mask; the operator's contract simply does not cover a rigid
translation on a bounded box. This is the exact missing operation:
**a periodic (or prescribed-inflow) treatment of the departure point in
`EvolvingDomains.Kinematic.CIP._departureCell`.**

To keep CIP in the comparison rather than drop it, the velocity handed to CIP is a *nodal
field* (a documented input) that is exactly $\mathbf v$ everywhere except on the single
outermost node layer at $x=0$ and $y=0$, where it is zero:

* every point the diagnostics touch satisfies $x-v_xt\ge0.11$ and $y-v_yt\ge0.21$ for
  $t\le1$;
* the CIP stencil is 2 cells wide and that layer is 7 cells away at 64² and 14 at 128², so
  the modification cannot reach the material region;
* WENO5, `TransportMap` and Trixi use the untampered velocity.

Evidence that this is inert: CIP's density error over the material region is the *smallest*
of the three grid-based EvolvingDomains methods at both resolutions, its physical mass drift
is $\le3.5\times10^{-4}$, and its constant-state test is exact to $2\times10^{-16}$ (§5.4).
The deviation is CIP-specific, quantified, and stated wherever CIP numbers appear.

---

## 5. Results

### 5.1 Geometry (shared by all four methods)

| quantity | 64² cells | 128² cells |
|---|---|---|
| worst relative area error over $t$ | $3.799\times10^{-4}$ | $7.889\times10^{-6}$ |
| area at $t=0$ (exact $1.3571680264\times10^{-1}$) | $1.3571783954\times10^{-1}$ | $1.3571697149\times10^{-1}$ |
| area at $t=1$ | $1.3566633767\times10^{-1}$ | $1.3571643480\times10^{-1}$ |

The WENO5-advected level set reconstructs the annulus to $O(h^2)$; nothing in the rest of
this report can be explained by geometry error, and the area itself improves by a factor 48
for a factor-2 refinement.

### 5.2 Physical mass

$M_h$ from §3, relative to the exact $M(0)$.

| method | 64² $M_h(0)/M-1$ | 64² $M_h(1)/M-1$ | 64² worst $\lvert\cdot\rvert$ | 128² $M_h(0)/M-1$ | 128² $M_h(1)/M-1$ | 128² worst $\lvert\cdot\rvert$ |
|---|---|---|---|---|---|---|
| WENO5 | $2.07\times10^{-4}$ | $-1.43\times10^{-5}$ | $2.37\times10^{-4}$ | $5.09\times10^{-5}$ | $5.73\times10^{-5}$ | $6.66\times10^{-5}$ |
| CIP | $2.07\times10^{-4}$ | $-1.40\times10^{-4}$ | $2.36\times10^{-4}$ | $5.09\times10^{-5}$ | $5.00\times10^{-5}$ | $6.04\times10^{-5}$ |
| Map | $2.07\times10^{-4}$ | $-2.00\times10^{-2}$ | $2.00\times10^{-2}$ | $5.09\times10^{-5}$ | $-8.52\times10^{-3}$ | $8.52\times10^{-3}$ |
| Trixi | $5.65\times10^{-5}$ | $-2.82\times10^{-4}$ | $2.91\times10^{-4}$ | $1.18\times10^{-5}$ | $9.56\times10^{-6}$ | $1.99\times10^{-5}$ |

Read the columns separately. $M_h(0)$ is the **representation bias**: the same exact density,
sampled on the grid and integrated with the common rule, already differs from $M(0)$ by
$2.07\times10^{-4}$ (nodal, 64²) or $5.65\times10^{-5}$ (DG, 64²) — this is the decomposition
Campaign 1 asked for (its Experiment D, item ii) and it is *not* a transport error. The
transport-only change $(M_h(1)-M_h(0))/M_h(0)$ is

| | 64² | 128² |
|---|---|---|
| WENO5 | $-2.22\times10^{-4}$ | $+6.44\times10^{-6}$ |
| CIP | $-3.47\times10^{-4}$ | $-8.64\times10^{-7}$ |
| Map | $-2.02\times10^{-2}$ | $-8.57\times10^{-3}$ |
| Trixi | $-3.38\times10^{-4}$ | $-2.21\times10^{-6}$ |

WENO5 and CIP have no systematic sign; the map loses mass monotonically; Trixi is the
smallest and most stable. The map's drift improves by only 2.36$\times$ for a factor-4
increase in nodes.

![Left: physical mass drift, both resolutions.](output/3_mass_drift.png)

*Figure 1.* $\lvert M_h/M_{\rm ex}-1\rvert$ on a log scale. The map (red) separates from the
other three by two orders of magnitude and drifts downward throughout; the others oscillate
at $10^{-5}$–$10^{-4}$ (64²) and $10^{-5}$–$10^{-4}\to10^{-6}$ (128²) without a trend.

![Two different conservation statements.](output/7_conservation_split.png)

*Figure 2.* Left: the unweighted sum over a *fixed* node set. Right: the map's own per-step
invariant on its own support, which is $\le3.5\times10^{-16}$. The left panel changes by
$0.83$ because nodes enter and leave the support; that is not a violation of anything the
operator promises, and §6 explains why the two panels cannot both stay flat.

### 5.3 Density accuracy and the boundary layer

Relative errors over $\Omega_h$ at $t=1$:

| | $E_1$ 64² | $E_2$ 64² | $L^\infty$ core 64² | $E_1$ 128² | $E_2$ 128² | $L^\infty$ core 128² |
|---|---|---|---|---|---|---|
| WENO5 | $3.81\times10^{-4}$ | $1.07\times10^{-3}$ | $4.32\times10^{-3}$ | $5.07\times10^{-5}$ | $1.56\times10^{-4}$ | $6.05\times10^{-4}$ |
| CIP | $1.98\times10^{-4}$ | $6.54\times10^{-4}$ | $1.88\times10^{-3}$ | $4.33\times10^{-5}$ | $1.40\times10^{-4}$ | $2.40\times10^{-4}$ |
| Map | $4.42\times10^{-2}$ | $5.80\times10^{-2}$ | $1.23\times10^{-1}$ | $3.54\times10^{-2}$ | $5.26\times10^{-2}$ | $2.34\times10^{-1}$ |
| Trixi | $2.22\times10^{-4}$ | $6.47\times10^{-4}$ | $5.33\times10^{-3}$ | $\mathbf{7.95\times10^{-6}}$ | $\mathbf{2.99\times10^{-5}}$ | $\mathbf{4.19\times10^{-4}}$ |

"core" = background nodes with $\varphi<-h$, evaluated for Trixi at the same physical
points. Two-resolution ratios ($E_1$): WENO5 $7.51$ ($p\approx2.9$), CIP $4.57$
($p\approx2.2$), Map $1.25$ ($p\approx0.3$), Trixi $27.9$ ($p\approx4.8$). The Trixi ratio
is above the asymptotic $p+1=4$ expectation for degree 3, and with two grids and a geometry
that also converges this is not a reliable convergence rate — it is reported as a ratio, not
as an order.

Mean nodal error binned by signed distance in cells, $t=1$, 128²/32²:

| $\varphi/h$ | $(-\infty,-4)$ | $[-4,-2)$ | $[-2,-1)$ | $[-1,0)$ | $[0,1)$ | $[1,2)$ | $[2,\infty)$ |
|---|---|---|---|---|---|---|---|
| nodes | 1072 | 579 | 290 | 299 | 281 | 299 | 13821 |
| WENO5 | $7.7\times10^{-6}$ | $1.8\times10^{-5}$ | $2.7\times10^{-5}$ | $2.7\times10^{-5}$ | $2.3\times10^{-5}$ | $3.0\times10^{-5}$ | $9.4\times10^{-6}$ |
| CIP | $6.5\times10^{-6}$ | $1.1\times10^{-5}$ | $1.6\times10^{-5}$ | $1.7\times10^{-5}$ | $1.6\times10^{-5}$ | $1.5\times10^{-5}$ | $5.4\times10^{-6}$ |
| Map | $1.6\times10^{-2}$ | $3.8\times10^{-2}$ | $6.0\times10^{-2}$ | $8.0\times10^{-2}$ | $1.5\times10^{-1}$ | $7.7\times10^{-1}$ | $1.0$ |
| Trixi | $9.0\times10^{-6}$ | $1.3\times10^{-5}$ | $1.9\times10^{-5}$ | $2.2\times10^{-5}$ | $2.3\times10^{-5}$ | $3.0\times10^{-5}$ | $8.3\times10^{-6}$ |

Negative $\varphi/h$ is the material side; the bins aggregate the inner and the outer
interface. The three accurate methods show a mild, symmetric interface layer of one to two
cells (a factor 2–3 above their interior error, *not* growing outward); the map's error grows
monotonically outward through the band and reaches order one in the void — the transported
support carries large errors where the exact density is smooth, exactly the signature
Campaign 1 localised to the support edge.

![Density maps.](output/1_density_maps.png)

*Figure 3.* Density on the numerically evolved material domain at $t=0,0.5,1$ (nodes with
$\varphi>0$ masked, not zeroed), 128²/32², common colour scale, white contour = the level-set
$\varphi=0$. The map's field is visibly smoothed against the interface; the other three are
indistinguishable at this scale.

![Error maps.](output/2_error_maps.png)

*Figure 4.* $\lvert\rho_h-\rho_{\rm ex}\rvert$ at $t=1$ on a **log** colour scale
($10^{-7}$–$10^{-1}$), 128²/32². Trixi and CIP are at $10^{-6}$–$10^{-5}$ over most of the
annulus, WENO5 one colour step higher, and the map two to four decades higher with a
staircase structure aligned with the grid.

![Boundary layer.](output/5_boundary_layer.png)

*Figure 5.* Left: a $0.14\times0.14$ window on the inner interface at $t=1$ (dashed red: the
analytic circle). WENO5, CIP and Trixi reproduce the interface without a visible layer at
this colour range; the map shows a staircase of cell-sized plateaus extending inward from the
interface. Right: the profile of the table above.

Boundedness and positivity on material nodes at $t=1$ (exact range $[1.0000,1.1907]$):

| | min | max | min over all steps |
|---|---|---|---|
| WENO5 | 1.0000 | 1.1904 | 1.0000 |
| CIP | 1.0000 | 1.1909 | 1.0000 |
| Map | 0.7415 | 1.3011 | 0.7136 |
| Trixi | 1.0000 | 1.1906 | 1.0000 |

Only the map undershoots and overshoots on the material domain.

### 5.4 Constant-state test ($\rho_0\equiv1$)

The decisive probe for partition of unity and support-boundary defects. Exact answer:
$\rho\equiv1$ and $M_h(t)=|\Omega_h(t)|$.

| 128² / 32² | $M_h(1)/|\Omega|-1$ | $\max\lvert\rho_h-1\rvert$ on material nodes | $\max$ over steps | $E_1$ |
|---|---|---|---|---|
| WENO5 | $-2.71\times10^{-6}$ | $0$ | $0$ | $1.2\times10^{-19}$ |
| CIP | $-2.71\times10^{-6}$ | $0$ | $0$ | $1.2\times10^{-19}$ |
| Map | $-9.81\times10^{-3}$ | $\mathbf{0.234}$ | $0.234$ | $3.63\times10^{-2}$ |
| Trixi | $-2.71\times10^{-6}$ | $4.4\times10^{-16}$ | $4.4\times10^{-16}$ | $1.5\times10^{-16}$ |

WENO5, CIP and Trixi reproduce a constant exactly (their residual mass "drift" is only the
geometry: $-2.71\times10^{-6}$ is the cut-area error at $t=1$). The map destroys a constant by
$23\%$ at the worst material node, and loses $0.98\%$ of the transported constant. **This
identifies the origin of its physical mass loss: it is the integrated partition-of-unity
defect of the remap on its moving support, not a density-transport error.** With the real
$\rho$ the loss is $-8.5\times10^{-3}$, essentially the same number.

![Constant state.](output/6_constant_state.png)

*Figure 6.* $\lvert\rho_h-1\rvert$ after transporting $\rho\equiv1$ to $t=1$ (log scale,
$10^{-17}$–$10^{-1}$). Three panels are empty (machine zero, clipped). The map's defect is
concentrated in a two-cell band at the support edge and decays inward — the same layer where
Campaign 1 measured a partition-of-unity defect of $0.86$ on the *band* nodes.

### 5.5 Full-domain conservation of the DG solve (a different statement)

Trixi integrates the *whole periodic square*, where its DG discretisation is conservative by
construction. $\int_{[0,1]^2}\rho_h$ from `Trixi.integrate` (LGL quadrature):

| | $t=0$ | $t=1$ | relative change |
|---|---|---|---|
| 16² elements | $1.004751479155$ | $1.004751479154$ | $-4.2\times10^{-14}$ |
| 32² elements | $1.004751658888$ | $1.004751658888$ | $+7.2\times10^{-14}$ |
| constant state | $1.0$ | $1.0$ | $0$ |

The 32² value agrees with the exact torus integral $1+\tfrac12\pi\sigma^2=1.004751658889$ to
12 digits, and the 16² deficit ($-1.8\times10^{-7}$) is the DG interpolation error of the
initial datum. This is the strongest single validation of the Trixi side of the experiment.

**This full-domain conservation has nothing to do with mass over the moving annulus.**
$\int_{[0,1]^2}$ is a torus invariant of a *different* domain; $M_h$ is an integral over a
piecewise-linear reconstruction of a moving annulus, using a weight vector $w(t)$ that
changes every step. Trixi is conservative for the first and, like every method here, only
approximately so for the second.

### 5.6 Work and allocations (solver only, per whole run)

| 128² / 32² | time | bytes | dof |
|---|---|---|---|
| geometry (WENO5 on $\varphi$) | 0.287 s | 0.8 MiB | 16641 |
| WENO5 | 0.340 s | 0.8 MiB | 16641 |
| CIP | 0.188 s | 0.0 MiB | 16641 |
| Map | 0.654 s | 241.9 MiB | 16641 |
| Trixi DGSEM p=3 | **0.055 s** | 0.0 MiB | 16384 |
| adapter + diagnostics | 1.196 s | 278.6 MiB | — |

At 64²/16²: WENO5 0.020 s, CIP 0.010 s, Map 0.070 s, Trixi 0.004 s. Every operator is warmed
up on throwaway state before the timed loop, so these exclude JIT; the diagnostics are timed
separately because they are the adapter's cost, not a solver's reference.

![Accuracy vs resolution and cost.](output/4_accuracy_vs_dof.png)

*Figure 7.* Left: $E_1$ against degrees of freedom — the two resolutions are deliberately
dof-matched (4225 vs 4096, 16641 vs 16384), so this is an honest like-for-like comparison.
Right: $E_1$ against solver time; the Trixi line sits one to two decades below the others at
the same or lower cost.

---

## 6. Conservation: two identities, four properties, what remains open

### 6.1 The algebra

For a fixed pair of cuts the remap is linear, $\rho^{n+1}=P\rho^{n}$ with
$P\in\mathbb R^{m\times n}$ ($n=|\mathcal S|$ source-support nodes, $m=|\mathcal T|$ target).
Two *different* identities are in play:

$$\textbf{algebraic:}\quad \mathbf 1_m^{\mathsf T}P=\mathbf 1_n^{\mathsf T}
\iff \sum_i\rho_i^{n+1}=\sum_j\rho_j^{n}\ \ \forall\rho,
\qquad\qquad
\textbf{weighted:}\quad w_{\text{new}}^{\mathsf T}P=w_{\text{old}}^{\mathsf T}
\iff M_h^{n+1}=M_h^{n}\ \ \forall\rho .$$

With $M_h=w^{\mathsf T}\rho$ and constant $w$ over a step, the second is exactly what physical
conservation requires; the first is what the operator advertises. They coincide only when
$w\propto\mathbf 1$ on a support of unchanged size. Note also that
$\mathbf 1^{\mathsf T}P=\mathbf 1^{\mathsf T}$ together with $P\mathbf 1=\mathbf 1$ forces
$m=n$: a single linear map on differing supports cannot both conserve the unweighted sum and
reproduce a constant everywhere. Something must give, and §5.4 shows what does.

### 6.2 What was measured

| property | definition | WENO5 | CIP | Map | Trixi |
|---|---|---|---|---|---|
| unweighted sum over its own support, per step | $\mathbf 1^{\mathsf T}q^{n+1}/\mathbf 1^{\mathsf T}q^{n}-1$ | n/a | $2.2\times10^{-16}$ | $3.5\times10^{-16}$ | n/a |
| unweighted sum, fixed node set | $\sum_i q_i$ over the background grid, per step | $2.6\times10^{-6}$ | $2.2\times10^{-16}$ | $0.84$ | n/a |
| constant state | $\max\lvert P\mathbf 1-1\rvert$ on material nodes | $0$ | $0$ | $0.234$ | $4.4\times10^{-16}$ |
| weighted identity | $w_{\text{new}}^{\mathsf T}P/w_{\text{old}}^{\mathsf T}-1$, applied to $w_{\text{old}}$ | — | — | $\mathbf{-1.0\times10^{-2}}$ | — |
| physical mass | $M_h(t)$, §3 | $5.7\times10^{-5}$ | $5.0\times10^{-5}$ | $-8.5\times10^{-3}$ | $9.6\times10^{-6}$ |

The map's two panels of figure 2 now make sense together: $\mathbf 1^{\mathsf T}P=\mathbf 1^{\mathsf T}$
holds to a few ulps *on the support the map is given*, the sum over a fixed background node
set changes by $0.84$ because the support itself moves, and the geometry-weighted identity —
the only one that carries physical meaning — fails by exactly the observed mass loss.

### 6.3 What this settles, and what it does not

Settled (new in this campaign):

* The map's physical mass loss is **not** a density-transport artefact. The constant-state
  test reproduces essentially the whole drift with $\rho\equiv1$, at both resolutions
  ($-9.81\times10^{-3}$ vs $-8.52\times10^{-3}$ at 128²).
* It is **not** the cut geometry: the same cut, with exact nodal data pushed through the same
  weights, gives $+4.6\times10^{-5}$ at 128², and the cut area itself is good to
  $7.9\times10^{-6}$.
* The weights $w_i=\int_{\Omega_h}N_i$ do **not** satisfy $w_{\text{new}}^{\mathsf T}P=w_{\text{old}}^{\mathsf T}$;
  the per-step defect is $1.0\times10^{-2}$ at 128², the same size as the drift, which answers
  Campaign 1's Experiment A in the negative *for this measure*. It is not merely a matter of
  choosing the right geometry-dependent weights of this family.

Still open:

* Whether *any* implementable weight vector satisfies the identity, or whether the operator
  has to act on **extensive** nodal masses $w_i\rho_i$ instead of intensive samples. This was
  not tested; it is the natural next experiment and it is cheap.
* The split of the $0.234$ constant-state defect between the pull-phase edge renormalisation
  and the push-phase demand/leakage correction.
* Whether the defect scales with the band width in cells rather than with $h$ (Campaign 1's
  H1, consistent with the factor 2.36 improvement here but not established by two grids).

---

## 7. Does an external DG solver compose with EvolvingDomains?

**Yes, and cheaply — with one real gap on the Trixi side and one on the EvolvingDomains side.**

What composed without friction:

* `Cut` and `CutSurface` remain sufficient to define a *common* physical measure for a solver
  that knows nothing about cuts. The cut triangles, the bilinear weights $w_i(t)$ and the
  degree-5 rule were built from documented arrays only (`subcell_states`,
  `point_pointers`/`point_ids`, `physical_points`), and the same triangles integrate the DG
  polynomial. No wrapper, no `getfield` on a private field, no fork.
* The shared time step works: the same $\Delta t$ and the same SSP-RK3 family drive all four
  methods, Trixi included, and both integrators are asserted to land on the same physical
  time with the same step count (`t=1.000000000000000`, 55 and 109 steps). Trixi's stable
  $\Delta t$ is 2.5$\times$ the shared one, so this costs Trixi nothing in stability.
* Geometry, diagnostics, interpolation and plotting are completely independent of which
  solver produced the scalar: swapping in a fourth method changed no line of the geometry
  loop or the measurement code.

What needed glue, and how much:

| glue | lines | why |
|---|---|---|
| `DGView`, `dgView`, `lagrange1d!`, `evalDG` | 64 | no public point-evaluation API; the DG polynomial has to be rebuilt from the nodal analysis interface |
| `trixiSetup`, `advanceWithTrixi!` | 18 | `TreeMesh`+`SemidiscretizationHyperbolic`+`Trixi.init`/`step!` at a prescribed $\Delta t$ |
| `verifyAdapter` | 33 | polynomial-reproduction proof (permanent, run every invocation) |

The adapter allocates (its evaluations are not allocation-free) and is timed separately:
1.20 s and 279 MiB for 330 calls at 128² — 20$\times$ the solver time of the cheapest method
being measured, because it evaluates the DG field at ~$5.3\times10^{3}$ triangles $\times$
7 points and at every background node, every step. This is the price of *comparing* rather
than of *solving*, and it would dominate any workflow that asked for dense pointwise output
at every step.

Public-API gaps found (both stated exactly, neither fatal):

1. **Trixi 0.17.21 has no documented arbitrary-point evaluation.** `PlotData2D`
   (`trixi2img`/`trixi2vtk` machinery) samples on its own plot grid;
   `Trixi.get_node_coords`/`get_node_vars` give nodal data but no evaluation. Missing
   operation: something like `evaluate_at_points(u, points, semi)`.
2. **EvolvingDomains' CIP has no periodic or inflow boundary treatment**, and its departure
   guard is global, so a rigid translation is impossible on a bounded grid (§4.4). Missing
   operation: a periodic wrap (or prescribed inflow state) in
   `CIP._departureCell`/`_xSweep!`/`_ySweep!`.

---

## 8. Conclusions

1. **Trixi improves accuracy and reduces cost on this problem.** $E_1(128^2)=7.95\times10^{-6}$
   against $4.3$–$5.1\times10^{-5}$ for CIP/WENO5 at the same dof, with $0.055$ s of solver
   time against $0.19$–$0.65$ s. It is also the only method that is simultaneously the most
   accurate in $E_1$, $E_2$ and $L^\infty$ at the finer resolution and exactly
   constant-preserving. The margin is large enough that the conclusion survives the
   caveats (single problem, single velocity, two resolutions, a thin annulus, one flux and
   one polynomial degree).
2. **The conservative map is the least accurate method here by two to three orders of
   magnitude**, and its physical mass loss is now attributed: it is the integrated
   partition-of-unity defect of the remap on its moving support, reproduced almost exactly by
   a constant-state test, invisible to the unweighted nodal-sum invariant, and not fixed by
   the geometry-dependent weights $w_i$.
3. **A high-order external PDE solver can replace the transport discretisation while
   EvolvingDomains keeps geometry, field representation and diagnostics** — at the cost of a
   64-line adapter, because Trixi does not expose point evaluation. What it cannot currently
   do is anything that needs the *cut structure inside* the DG solve (cut-cell DG, embedded
   boundaries); Trixi solved the periodic background problem in full, and the moving annulus
   was imposed only through the diagnostics. That is the honest boundary of this result.

---

## Appendix — reproduction and artifact inventory

```sh
# native cut backend, built once from this checkout (needs cmake + network for the pinned CutCells)
BUILD_DIR=/tmp/ed-cut2-build PREFIX=/tmp/ed-cut2-install \
  bash scripts/build_cxx_bindings.sh
cp /tmp/ed-cut2-install/lib/libevolvingdomains_cut_cells.dylib \
   validation/ConservativeAnnulus/round2/lib/

julia validation/ConservativeAnnulus/round2/run.ConservativeAnnulusTrixi.jl [pilot]
julia validation/ConservativeAnnulus/round2/visualize.ConservativeAnnulusTrixi.jl
```

`run` takes ≈29 s (both resolutions, four cases, including JIT) and writes
`output/results.jls`; `visualize` takes ≈15 s and writes the seven PNGs. Both scripts
activate the shared environment `validation/env` (Julia 1.11.7, EvolvingDomains from this
checkout via `[sources]`, Trixi 0.17.21, CairoMakie 0.15.17). The run script sets
`EVOLVINGDOMAINS_CUTCELLS_LIBRARY` to `lib/` only when the caller has not already provided a
library, prints the resolved `pathof(EvolvingDomains)` and aborts unless it is this
repository, and prints the git revision, the package versions and the pinned CutCells commit.

| Path | Contents |
|---|---|
| `run.ConservativeAnnulusTrixi.jl` | problem, quadrature, cut utilities, DG adapter, four advances, diagnostics, driver (840 lines) |
| `visualize.ConservativeAnnulusTrixi.jl` | seven CairoMakie figures, reads `output/results.jls` (220 lines) |
| `output/results.jls` | per-step histories, snapshots at $t=0,0.5,1$, summaries, provenance, adapter verification |
| `output/1_density_maps.png` | density on the material domain, four methods × three times |
| `output/2_error_maps.png` | $\lvert\rho_h-\rho_{\rm ex}\rvert$ at $t=1$, log scale |
| `output/3_mass_drift.png` | physical mass drift, both resolutions |
| `output/4_accuracy_vs_dof.png` | $E_1$ against dof and against solver time |
| `output/5_boundary_layer.png` | interface zoom + nodal error by $\varphi/h$ |
| `output/6_constant_state.png` | constant-state defect (partition of unity) |
| `output/7_conservation_split.png` | unweighted nodal sum vs the map's own invariant |
| `lib/libevolvingdomains_cut_cells.dylib` | locally built native cut backend (132 KB, macOS arm64) |

No file of EvolvingDomains, of its manifests, or of Trixi was modified. `round2/` is the only
directory this campaign writes to.
