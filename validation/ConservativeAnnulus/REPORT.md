# ConservativeAnnulus — Conservative cut-to-cut transport on a translating annular domain

**Campaign 1 report.** Experiment under `validation/ConservativeAnnulus/` against
EvolvingDomains.jl (this checkout) in the pinned environment `validation/env`.

| Item | Value |
|---|---|
| Problem | $\partial_t\rho+\nabla\cdot(\rho\mathbf v)=0$ on a translating annulus, $t\in[0,1]$ |
| Domain | $\Omega(0)=\{0.12<\|\mathbf x-(0.35,0.45)\|<0.24\}$, $\mathbf v=(0.24,0.10)$ |
| Tested operator | `TransportMap` + `lagrangianTransport!` (Lentine–Grétarsson–Fedkiw conservative cut-to-cut map) |
| Resolutions | 65² nodes (64² cells) and 129² nodes (128² cells), nested by `refine!` |
| Steps to $t=1$ | 55 and 109 (CFL 0.40), plus a CFL sweep {0.2, 0.4, 0.8} |
| Entry point | `run.ConservativeAnnulus.jl`, figures by `visualize.ConservativeAnnulus.jl` |
| Julia | 1.11.7 |

---

## 1. Executive summary

The experiment transports a smooth density on a translating annular material domain, whose
geometry is carried by a level-set `Field` and reconstructed per step by `cut`. The
conservative remap `TransportMap` / `lagrangianTransport!` is exercised on the flattened
background nodal arrays, and physical mass is measured independently by exact per-triangle
quadrature over the cut subcells.

The geometric reconstruction is accurate and converges strongly: the relative area error
falls from $3.80\times10^{-4}$ to $7.89\times10^{-6}$ and the absolute centre error from
$2.02\times10^{-5}$ to $1.89\times10^{-6}$ between the two grids.

The transport map preserves its documented invariant — the unweighted nodal sum over its
source support — to at most $4.4\times10^{-16}$ relative, i.e. a few units in the last place,
at every step of both runs. This does not imply conservation of physical mass. Integrating
$\rho$ over the reconstructed geometry gives a drift of $-1.999\times10^{-2}$ on the quick grid
and $-8.512\times10^{-3}$ on the fine grid, while the same quadrature applied to the exact
density on the *same* cut geometry drifts by only $-3.78\times10^{-4}$ and
$-6.59\times10^{-6}$. The density error is concentrated in the two-node edge layer of the
transported support, where the map's own partition-of-unity defect reaches $0.86$; deep inside
the annulus it is an order of magnitude smaller.

Every operation was performed through exported names and documented data: `Grid`,
`Discretization`, `Field`, `refine!`, the hierarchy queries, `transfer!`, `cut`,
`eulerianTransport!`, `TransportMap`, `lagrangianTransport!`, and the owned arrays of `Cut`
and `CutSurface`. **The numerical origin of the physical mass loss has not been conclusively
established.** Three candidate mechanisms are identified and a controlled decomposition is
specified for Campaign 2.

---

## 2. Mathematical formulation

### 2.1 Domain, velocity, exact solution

Let $D=[0,1]^{2}$ and $t\in[0,1]$. The material domain is the rigid translate

$$\Omega(t)=\Omega(0)+t\mathbf v,\qquad \mathbf v=(0.24,\,0.10)^{\mathsf T},$$

$$\Omega(0)=\bigl\{\mathbf x\in D:\; r_{\mathrm{in}}<\|\mathbf x-\mathbf c_0\|<r_{\mathrm{out}}\bigr\},
\qquad \mathbf c_0=(0.35,0.45),\quad r_{\mathrm{in}}=0.12,\quad r_{\mathrm{out}}=0.24 .$$

$\Omega(t)$ stays strictly inside $D$ for $t\le 1$ (the outer radius reaches
$0.79$ in $y$ and $0.83$ in $x$). The level-set representation used to cut the grid is

$$\phi(\mathbf x,t)=\max\bigl(r_{\mathrm{in}}-\|\mathbf x-\mathbf c_0-t\mathbf v\|,\;
\|\mathbf x-\mathbf c_0-t\mathbf v\|-r_{\mathrm{out}}\bigr),$$

negative on the material, with $\phi_0(\mathbf x)=\phi(\mathbf x,0)$.

The density obeys

$$\partial_t\rho+\nabla\cdot(\rho\mathbf v)=0\quad\text{in }\Omega(t),
\qquad
\rho_0(\mathbf x)=1+\tfrac12\exp\!\Bigl(-\frac{\|\mathbf x-(0.41,0.48)\|^{2}}{0.055^{2}}\Bigr).$$

Because $\mathbf v$ is constant, the exact solution is the translate of the initial datum,

$$\rho(\mathbf x,t)=\rho_0(\mathbf x-t\mathbf v),$$

and, since $\mathbf v\cdot\mathbf n=0$ on $\partial\Omega(t)$ in the frame moving with the
domain, the physical mass is constant:

$$M(t)=\int_{\Omega(t)}\rho(\mathbf x,t)\,d\mathbf x=M(0)=M_{\mathrm{exact}},
\qquad M_{\mathrm{exact}}=1.36320840\times10^{-1},$$

$$|\Omega(t)|=\pi\bigl(r_{\mathrm{out}}^{2}-r_{\mathrm{in}}^{2}\bigr)=1.35716803\times10^{-1},
\qquad \mathbf x_{c}(t)=\mathbf c_0+t\mathbf v .$$

$M_{\mathrm{exact}}$ is evaluated in polar coordinates about $\mathbf c_0$: 5-point
Gauss–Legendre on 20 radial subintervals of $[r_{\mathrm{in}},r_{\mathrm{out}}]$ times a
1024-point midpoint sweep in $\theta$ (`exactMass` in `run.ConservativeAnnulus.jl`). The
integrand is analytic; both rules converge spectrally, so the reference value is converged far
beyond the accuracy of the diagnostics. The peak of $\rho_0$ sits at radius $0.067$ from
$\mathbf c_0$, i.e. **inside the hole**, so the density variation actually present on the
material is moderate: its exact range on the in-phase nodes is $[1.0000,1.1848]$ at $t=0$
and $[1.0000,1.1907]$ at $t=1$ on the fine grid.

### 2.2 Discrete representation

The interface is carried by a nodal level-set field and reconstructed per step as a cut-cell
decomposition: each background cell is classified `CUT_IN` ($\phi<0$ at all four corners),
`CUT_OUT` ($\phi>0$ at all four corners) or `CUT_CUT` (mixed), and each cut cell is
subdivided into simplices — triangles in 2D — whose union is the in-phase part of the cell.
`Cut` stores these as owned flat arrays: `subcell_states` (per subcell, `CUT_IN` or
`CUT_OUT`), `cell_to_bgcell` (parent cell, sorted), `point_pointers`/`point_ids` (CSR
connectivity into `physical_points`), plus a `CutSurface` with the interface facets, their
normals and their parent cells.

A nodal density $q_i$ is carried on the background nodes. In this report the physical mass is
estimated from the cut subcells, not from the nodes:

$$M_h(t)=\sum_{K\,:\,\texttt{CUT\_IN}}|K|\;\rho_h(\mathbf x_K),\qquad
\mathbf x_K=\tfrac13\bigl(\mathbf p_1^{K}+\mathbf p_2^{K}+\mathbf p_3^{K}\bigr),$$

with $|K|$ the exact area of the physical triangle and $\rho_h$ the multilinear interpolant of
the nodal samples (`linearInterpolator`). The rule is exact for a linear density and second
order for a smooth one. It is a quadrature **weighted by physical cell areas**; the transport
operator's nodal sum is a different functional and is measured separately everywhere in this
report.

---

## 3. Software architecture and implementation

### 3.1 Objects used, and what they are

| Object | Kind | Role in the experiment |
|---|---|---|
| `Grid((xaxis, yaxis))` | immutable struct, coordinate vectors + `mesh_id` | the uniform $65^2/129^2/257^2$ Cartesian grids |
| `Discretization` | `Dict{Symbol,Grid}` with a shared-box invariant | the multilevel container |
| `Field` | `Dict{Symbol,Array{T,2}}` + its `Discretization` | level set $\phi$, density $\rho$, transfer study |
| `Cut` / `CutSurface` | owned flat `Vector`s (§2.2) | geometry snapshot at $t^-$ and $t^+$ |
| `TransportMap` | immutable, geometry-only | the discrete flow between two cuts |
| `WENO5Cache`, vectors, matrices | plain Julia | workspace, nodal vectors, diagnostics |

Levels are created and navigated exactly as documented:

```julia
d = Discretization(Dict(:coarse => Grid((axis, copy(axis)))))   # 65² nodes
refine!(d, :coarse, :fine;  factor = 2)                          # 129²
refine!(d, :fine,   :finer; factor = 2)                          # 257²
coarsest(d), finest(d)                                           # :coarse, :finer
immediateFiner(d, :coarse), immediateCoarser(d, :fine)           # :fine, :coarse
```

`refine!` preserves the physical box and produces exactly nested nodes; that nesting is what
`transfer!` requires and what makes the coarse/fine comparison in §4.4 well posed.

### 3.2 Discrete layout: how the cut relates to the arrays

Both `Field`s live on the same `Discretization`, so the level-set field and the density field
share the same grid identity (required by `TransportMap` and by the `Field` overload of
`lagrangianTransport!`). The map does **not** act on "cut degrees of freedom" in the sense of
the cut's own point ordering. It acts on the background node vector

$$q=\operatorname{vec}\bigl(\texttt{field.data[level]}\bigr)\in\mathbb R^{N},\qquad N=n_xn_y,$$

whose index is $i+(j-1)n_x$ — first axis fastest, which is exactly the node order
`eulerianTransport!` documents for its velocity vector ("one `SVector` per node, first axis
fastest") and the cell order `CutSurface` documents for `facet_to_bgcell`. The cut contributes
only a **support**: the mask of nodes that are corners of some `CUT_IN` or `CUT_CUT`
background cell. That mask is computed here from `bg_states` alone (`supportMask`, 12 lines),
and its equivalence with the map's own supports is asserted on every step of every run:

```julia
@assert map.source_mask == old_support
@assert map.target_mask == new_support
```

No step of any run (55 + 109 + 28 + 2 control steps, plus the warm-up and sensitivity runs)
failed these assertions. The cut's point ordering is therefore never assumed: it is used only
for geometry, via `point_pointers`/`point_ids`/`physical_points` (§2.2).

The source vector handed to the operator is the current nodal density zeroed outside the old
support. This is required by the documented semantics of the map — "only the map's source
support is transported, and everything else is set to zero, so a field defined everywhere
loses whatever lies outside the cut" (`TransportMap` docstring) — and it is harmless here
because every node with $\phi<0$ is by construction a corner of a cell that is not `CUT_OUT`,
hence inside the support. The dropped entries lie in the void.

### 3.3 The simulation loop

Per step, the mathematical content of §2 is visible line by line (`simulate` in
`run.ConservativeAnnulus.jl`):

```text
levelset ← Field(d, φ₀ sampled on nodes)          density ← Field(d, ρ₀ sampled on nodes)
velocity ← vec([SVector(vx, vy) for each node])   # 1. prescribed velocity, frozen per step
c_old    ← cut(levelset, level)                   # 2. old cut
S_old    ← supportMask(c_old)                     #    source support
q_old    ← vec(density.data[level]);  q_old[complement(S_old)] .= 0.0
eulerianTransport!(levelset, velocity, dt)        # 3. move the level set (WENO5 + SSP-RK3)
c_new    ← cut(levelset, level)                   # 4. new cut
S_new    ← supportMask(c_new)
map      ← TransportMap(c_old, c_new, grid, velocity, dt)   # 5. conservative map
q_new    ← zeros(N); lagrangianTransport!(q_new, q_old, map) # 6. remap
density.data[level] .= reshape(q_new, nₓ, n_y)
diagnostics(c_old, c_new, q_old, q_new, map, t)   # 7. every number in §4
```

Two diagnostics reuse the *same* map, so they cost one extra transport each and no extra
geometry work:

* a constant-state probe, `lagrangianTransport!(·, Float64.(S_old), map)`, which measures
  $P\mathbf 1$ (§5);
* the per-step "exact density on the same cut" integral, which isolates the geometry.

Both cuts are constructed every step (rather than carrying the previous step's cut) so that
the loop reads like §3 of the mission; the redundant `cut` roughly doubles the geometric cost
and is noted as such.

### 3.4 Adaptation code, and what was *not* needed

Everything numerical in the experiment is ordinary Julia on plain arrays; EvolvingDomains
supplies containers, geometry and operators, and no wrapper layer was written around them. The
project-specific code is:

| Piece | Lines | Nature |
|---|---|---|
| `supportMask` (cut cells → nodal support) | 12 | bookkeeping the package does not expose as a function |
| `subcellVertices`, `triangleArea`, `cutIntegrals` | 28 | physical quadrature over `Cut` arrays |
| `cutTriangles`, `cutSegments` | 19 | geometry export for plotting |
| `exactMass`, `layerProfile` | 40 | reference integral, error binning |
| `simulate` loop + records | ≈130 | the experiment itself |
| `transferStudy`, report, entry point | ≈170 | diagnostics and I/O |

The run script is 584 lines including comments and the printed report; no source file of
EvolvingDomains, and no manifest, was modified. The only nontrivial piece of "glue" is
`supportMask`, and it is a mapping from *documented cell states* to node indices — not a
reconstruction of anything hidden. Confirming that it agrees with the operator required no
inspection of package internals: the map publishes `source_mask`/`target_mask` as fields of an
exported struct, and those are used exclusively in the assertions above.

---

## 4. Numerical results

Error definitions used throughout: **area error** $=|\hat A|/|\Omega(t)|-1$ (relative);
**centre error** $=\|\hat{\mathbf x}_c-\mathbf x_c(t)\|$ (absolute, physical units);
**mass drift** $=M_h(t)/M_{\mathrm{exact}}-1$ (relative to the exact constant);
**nodal drift** $=\sum q^{n+1}/\sum q^{n}-1$ (relative, per step);
**density error** against $\rho_0(\mathbf x-t\mathbf v)$ at the nodes, split by signed
distance in node units $\phi/h$; **core** means nodes of the transported support with
$\phi<-h$.

### 4.1 Geometry

| quantity (worst over all steps) | 64² cells | 128² cells |
|---|---|---|
| relative area error | $3.799\times10^{-4}$ | $7.889\times10^{-6}$ |
| absolute centre error | $2.024\times10^{-5}$ | $1.889\times10^{-6}$ |
| area at $t=1$ (exact $1.35716803\times10^{-1}$) | $1.35666338\times10^{-1}$ | $1.35716435\times10^{-1}$ |
| centre at $t=1$ (exact $(0.590,0.550)$) | $(0.589982,\,0.549998)$ | $(0.590000,\,0.549999)$ |

Both error measures decrease by roughly two orders of magnitude for a factor-two refinement.
The coarse area error is *not* a random walk: it decreases systematically over the interval,
with small oscillations (figure 3, left), i.e. the level-set advection erodes the domain slightly, and the sign and
size of that erosion is itself resolution dependent. The centre error stays at both
resolutions far below a cell ($h=1.56\times10^{-2}$ and $7.81\times10^{-3}$), so the
translated geometry is captured essentially exactly.

![Cut geometry at t = 0, 0.49, 1.00. Left: interface cut cells in a 0.1×0.1 window on the outer interface, with the Cartesian grid and the in-phase subcell triangles; right: full box. Solid blue is the reconstruction from `Cut`, dashed orange the analytic circles.](output/1_cut_geometry.png)

*Figure 1.* The cut-cell subtriangulation tiles the annulus exactly, and the reconstructed
interface is a piecewise-linear interpolant of the level-set zero contour that lies on the
exact circle to plotting accuracy at $t=0$, is very slightly outside at $t=0.49$ and inside at
$t=1$ — visible as the sub-pixel doubling of the blue and orange curves, and consistent with
the sign change of the area error in figure 3. The number of triangles is large because every
wholly in-phase cell also emits two triangles (8560 subcells for the $65^2$ grid at $t=0$, of
which 1288 are `CUT_IN`).

![Density on the material domain at t = 0, 0.49, 1.00, masked where φ > 0. Colour range [0.90, 1.25]; colorbar.](output/2_density.png)

*Figure 2.* The remapped density inside the annulus (nodes with $\phi>0$ are masked, not
zeroed). The initial bump — centred inside the hole, so only its tail is material — is
progressively smoothed, and near both interfaces a band of low values ($\approx0.75$–$0.9$)
appears; this is the visible signature of the support-edge layer analysed in §4.3 and §6.

![Left: relative area error vs t for both resolutions. Right: geometric centre error, log scale.](output/3_geometry_errors.png)

*Figure 3.* Cut-reconstruction accuracy against the analytic annulus. The fine-grid errors are
flat and small ($\sim10^{-5}$ relative area, $\sim10^{-6}$ centre); the coarse-grid area error
drifts to $-3.7\times10^{-4}$ over the interval. No monotone convergence study beyond two
resolutions was performed.

### 4.2 Conservation of the nodal sum, and physical mass

| quantity | 64² cells | 128² cells |
|---|---|---|
| nodal sum $\sum q$ at $t=0$ | $7.415965899051855\times10^{2}$ | $2.611025810780435\times10^{3}$ |
| nodal sum at $t=1$ | $7.415965899051821\times10^{2}$ | $2.611025810780482\times10^{3}$ |
| worst per-step relative drift | $3.33\times10^{-16}$ | $4.44\times10^{-16}$ |
| $h^{2}\sum q$ at $t=0$ vs $M_h(0)$ | $+32.8\%$ | $+16.9\%$ |
| $M_h(0)$ (quadrature) | $1.363481\times10^{-1}$ | $1.36328\times10^{-1}$ |
| $M_h(1)$ | $1.33595\times10^{-1}$ | $1.35160\times10^{-1}$ |
| **mass drift** $M_h(1)/M_{\mathrm{exact}}-1$ | $\mathbf{-1.999\times10^{-2}}$ | $\mathbf{-8.512\times10^{-3}}$ |
| same quadrature, exact $\rho$ at the centroids | $-3.777\times10^{-4}$ | $-6.585\times10^{-6}$ |
| same quadrature, void-side nodes taken from $\rho_{\mathrm{exact}}$ | $-2.155\times10^{-2}$ | $-9.291\times10^{-3}$ |

Three independent facts should be read separately.

1. The operator's invariant holds to floating-point precision: at most 4 units in the last
   place of the nodal sum, at every step of both runs (figure 4, left). Taken at face value
   this is exactly what the docstring promises.
2. The invariant is *not* a mass statement. Even at $t=0$, with the exact initial datum, the
   uniform-weight Riemann sum $h^{2}\sum q$ exceeds the physical mass by $32.8\%$ (64²) and
   $16.9\%$ (128²): the support contains about 1.5 times as many nodes as the material region
   holds ($741$ vs $485$ at $64^2$, $2591$ vs $2070$ at $128^2$, measured at $t=0$), so a
   uniform weight counts halo nodes that integrate almost nothing. This mismatch is $O(h)$,
   unlike the $O(h^{2})$ geometry error.
3. The physical mass does drift, by $2.0\%$ and $0.85\%$. The drift is not a geometry error:
   the same quadrature with the exact density *on the same cut* is $50$–$1300$ times smaller.

![Left: nodal-sum drift in units of eps(Σq), linear scale, for both resolutions. Right: physical mass drift M_h/M_exact−1 for both resolutions, with the exact-density-on-the-same-cut control (dotted).](output/4_conservation.png)

*Figure 4.* The two conservation statements side by side. Left: the map's invariant is exact
to a few ulps. Right: the physical integral is not, and its error has a fine structure — a
mean downward trend plus a short-period oscillation, consistent with the support gaining and
losing nodes as the annulus slides through the lattice. The dotted controls integrate the analytic
density over the *same* reconstructed triangles and stay at $10^{-4}$ (coarse) and $10^{-6}$
(fine).

### 4.3 Density accuracy, positivity, boundedness

| quantity | 64² cells | 128² cells |
|---|---|---|
| $L^\infty$ on core nodes ($\phi<-h$) | $1.228\times10^{-1}$ | $2.345\times10^{-1}$ |
| $L^{1}$ on core nodes | $3.671\times10^{-2}$ | $2.900\times10^{-2}$ |
| $L^\infty$ over all nodes (incl. halo) | $1.493$ | $1.499$ |
| range on in-phase nodes at $t=1$ | $[0.7705,\,1.1627]$ | $[0.7415,\,1.3011]$ |
| exact range on the same nodes | $[1.0000,\,1.1848]$ | $[1.0000,\,1.1907]$ |
| overshoot / undershoot on in-phase nodes | $-2.20\times10^{-2}$ / $2.30\times10^{-1}$ | $+1.10\times10^{-1}$ / $2.59\times10^{-1}$ |
| range on the whole support (halo included) | $[0.6555,\,1.9260]$ | $[0.6643,\,1.7911]$ |
| positivity: min over steps, in-phase nodes | $0.6997$ | $0.7136$ |
| positivity: min over steps, support nodes | $0.4834$ | $0.4805$ |
| constant state: $\max\|\cdot-1\|$ on target support | $0.8622$ | $0.8611$ |
| constant state: $\max\|\cdot-1\|$ on core nodes | $0.0396$ | $0.0394$ |

Positivity survives in the material (the field never goes negative) but boundedness fails:
values reach $1.79$–$1.93$ on void-side nodes where the exact density is at most $1.19$, and
the transported field overshoots its exact maximum by $11\%$ on in-phase nodes of the fine
grid. The maximum over *all* nodes ($\approx1.49$) is attained on those void-side nodes and is
therefore not a statement about the material field; it is reported because the same
quantities are what `Σq` and the pull stencils contain.

The error is not distributed uniformly. Binning the final nodal error by $\phi/h$:

| $\phi/h$ bin | 64² nodes | 64² $L^\infty$ | 64² $L^{1}$ | 128² nodes | 128² $L^\infty$ | 128² $L^{1}$ |
|---|---|---|---|---|---|---|
| $[-6,-4)$ | — | — | — | 649 | $7.904\times10^{-2}$ | $1.843\times10^{-2}$ |
| $[-4,-2)$ | 254 | $9.734\times10^{-2}$ | $3.081\times10^{-2}$ | 579 | $1.677\times10^{-1}$ | $3.794\times10^{-2}$ |
| $[-2,-1)$ | 158 | $1.228\times10^{-1}$ | $4.619\times10^{-2}$ | 290 | $2.345\times10^{-1}$ | $6.008\times10^{-2}$ |
| $[-1,0)$ | 149 | $2.295\times10^{-1}$ | $6.674\times10^{-2}$ | 299 | $3.011\times10^{-1}$ | $7.973\times10^{-2}$ |
| $[0,1)$ (void side) | 145 | $9.260\times10^{-1}$ | $1.440\times10^{-1}$ | 281 | $7.911\times10^{-1}$ | $1.502\times10^{-1}$ |
| $[1,2)$ (void side) | 143 | $1.401$ | $7.895\times10^{-1}$ | 299 | $1.295$ | $7.699\times10^{-1}$ |

Two features are worth separating. The **deep** bins improve with refinement
($L^{1}$ at $[-4,-2)$: $3.08\times10^{-2}\to1.84\times10^{-2}$), while the **innermost
material bin** gets worse pointwise and in the mean ($L^\infty$
$2.30\times10^{-1}\to3.01\times10^{-1}$, $L^{1}$ $6.67\times10^{-2}\to7.97\times10^{-2}$). A
boundary layer whose thickness is fixed in *cells* rather than in physical length would
produce exactly this: refining does not move the layer, it only resolves it better, and the
node adjacent to the interface is closer to the interface. This annulus is thin — 7.7 cells
wide at 64², 15 at 128² — so the layer occupies a large share of the material.

Per-step versus accumulated error, from one-step control runs started from exact data:

| level | one-step $L^{1}$ core | one-step $L^\infty$ core | after $n$ steps | ratio |
|---|---|---|---|---|
| 64² | $4.571\times10^{-4}$ | $1.834\times10^{-2}$ | $3.671\times10^{-2}$ (55 steps) | 80.3 |
| 128² | $4.326\times10^{-4}$ | $2.612\times10^{-2}$ | $2.900\times10^{-2}$ (109 steps) | 67.0 |

The single-step error does not decrease between the two grids, and the accumulated error is
sublinear in the number of steps. Both observations point to a resolution-independent
per-step defect (the support-edge layer) superimposed on a per-step error that does converge
in the interior.

### 4.4 Multilevel transfer

Prolongating the analytic field $g(x,y)=\exp(\sin 2\pi x+\cos 2\pi y)$ with `transfer!`
(default `linearInterpolator`):

| direction | $h$ | $L^\infty$ | $L^{1}$ |
|---|---|---|---|
| coarse → fine (65² → 129²) | $1.563\times10^{-2}\to7.813\times10^{-3}$ | $1.767\times10^{-2}$ | $1.161\times10^{-3}$ |
| fine → finer (129² → 257²) | $7.813\times10^{-3}\to3.906\times10^{-3}$ | $4.443\times10^{-3}$ | $2.905\times10^{-4}$ |
| observed order in $L^\infty$ | | $1.99$ | |

The restriction direction returns $L^{\infty}=0$ exactly, which is the expected node sampling
on nested grids and is *not* a conservative averaging restriction — the docstring's warning
("coarser→finer is prolongation, finer→coarser is point-sampling, which is not conservative")
is reproduced numerically. Prolongating the *quick run's transported density* to 129² and
comparing with the fine run gives $L^\infty=1.497\times10^{-1}$ and
$L^{1}=1.861\times10^{-2}$ on fine core nodes: the coarse truncation error is carried upward,
not removed, again as documented.

![Left: error of transfer! coarse→fine against the analytic field. Right: prolongated coarse density minus fine density at t = 1 on fine core nodes.](output/5_transfer.png)

*Figure 5.* Coarse/fine interpolation. The analytic transfer error is a smooth $O(h^{2})$
pattern with the largest errors where the field has the most curvature. The run-to-run
difference of the transported density is an order of magnitude larger, dominated by the same
support-edge structure as figure 2.

### 4.5 Timestep sensitivity, runtime, allocations

Quick grid, CFL sweep (`dt = \mathrm{CFL}\cdot h/\|\mathbf v\|_1`, $\|\mathbf v\|_1=0.34$):

| CFL | steps | worst nodal drift | mass drift at $t=1$ | $L^{1}$ core | area error at $t=1$ |
|---|---|---|---|---|---|
| 0.2 | 109 | $3.33\times10^{-16}$ | $-1.90\times10^{-2}$ | $3.85\times10^{-2}$ | $-3.97\times10^{-4}$ |
| 0.4 | 55 | $3.33\times10^{-16}$ | $-2.00\times10^{-2}$ | $3.67\times10^{-2}$ | $-3.72\times10^{-4}$ |
| 0.8 | 28 | $2.22\times10^{-16}$ | $-2.15\times10^{-2}$ | $2.66\times10^{-2}$ | $-1.53\times10^{-4}$ |

The invariant is timestep independent; the mass drift is nearly so (it changes by about
12% of its own value across a factor-four change in $\Delta t$), while the density error improves with *fewer* steps — the
expected signature of a scheme whose per-step dissipation is the dominant error and whose
total dissipation scales with the number of steps. CFL 0.4 was chosen as the reporting
configuration and is the value in all other tables.

| level | loop wall time | bytes allocated | GC time | steps |
|---|---|---|---|---|
| 64² | 0.30 s | 545 MiB | 0.07 s | 55 |
| 128² | 1.99 s | 4285 MiB | 0.46 s | 109 |

Measured with `@timed` around the loop only, after a warm-up run that compiles the cut and
transport path; the allocations are dominated by cut construction (two cuts per step, each
allocating fresh point, subcell and facet arrays) rather than by the remap. Excluded from
these numbers are the transfer study (which builds a 257² level), the two single-step control
runs and the CFL sweep; a full invocation of `run.ConservativeAnnulus.jl` was observed at
about 10 s wall clock, including Julia start-up and compilation.

---

## 5. The meaning of conservation

Four properties are routinely conflated. This section states them precisely and separates what
was measured.

### 5.1 The map as a linear operator

For a fixed pair of cuts, every scale in the operator depends only on geometry: the pulled
stencil weights, the supersample share $2^{-N}$, the "demand" denominator
$d_j=\sum_{\text{rays}}w_{j}$, and the leakage fractions $1-d_j$. Hence the remap is an
$m\times n$ matrix for fixed geometry,

$$q^{n+1}=P\,q^{n},\qquad P\in\mathbb R^{m\times n},$$

where $n=|\mathcal S|$ is the number of source-support nodes and $m=|\mathcal T|$ the number of
target-support nodes. (In this implementation $m=n=N$ on the embedding space, but rows and
columns outside the supports are identically zero and the source is zeroed outside
$\mathcal S$, so the *effective* map is the restriction displayed above. $|\mathcal S|$ and
$|\mathcal T|$ differ as soon as the band gains or loses nodes.)

**Algebraic conservation** of the unweighted nodal quantity is the row-vector identity

$$\mathbf 1_m^{\mathsf T}P=\mathbf 1_n^{\mathsf T}
\qquad\Longleftrightarrow\qquad
\sum_i q_i^{n+1}=\sum_j q_j^{n}\ \ \forall q .$$

**Constant-state preservation** (partition of unity on the target) is the column identity

$$P\mathbf 1_n=\mathbf 1_m .$$

**Observation.** If both hold, then

$$n=\mathbf 1_n^{\mathsf T}\mathbf 1_n=\mathbf 1_m^{\mathsf T}P\mathbf 1_n=\mathbf 1_m^{\mathsf T}\mathbf 1_m=m .$$

So a single linear map cannot simultaneously conserve the unweighted sum *and* reproduce a
unit constant on the whole of $\mathbb R^n$ and $\mathbb R^m$ when $m\ne n$. This is an
algebraic statement about linear maps, not a defect claim: it says that as soon as the source
and target supports differ in size, at least one of the two properties must fail somewhere,
and a scheme has to choose (or restrict the domain of each property — e.g. reproduce constants
only where the reconstruction stencil is not truncated).

The measurements in §4.2/§4.3 are exactly this trade-off: the row identity holds to a few ulps
by construction, and the constant-state defect is *not* distributed uniformly but concentrated
at the support edge ($0.86$ on the band, $0.039$ one node in), vanishing algebraically inward.
That is the compromise the implementation makes, and it is visible in the data rather than
assumed.

### 5.2 What weighted conservation would require

Let $w\in\mathbb R^{m}$ and $\tilde w\in\mathbb R^{n}$ be geometry-dependent integration
weights for the target and source, so that the physical mass is approximated by
$M_h=w^{\mathsf T}\rho_h$. With $\rho^{n+1}=P\rho^{n}$ and constant weights over the step,

$$M_h^{n+1}=w^{\mathsf T}P\,\rho^{n}\overset{!}{=}\tilde w^{\mathsf T}\rho^{n}
\qquad\Longleftrightarrow\qquad
w^{\mathsf T}P=\tilde w^{\mathsf T}.$$

Since the operator instead supplies $\mathbf 1^{\mathsf T}P=\mathbf 1^{\mathsf T}$, weighted
conservation would be automatic only if $w=\tilde w=\text{const}$ on a support of unchanged
size — precisely the case excluded when nodes enter and leave the band. Two routes are open:
supply weights to the map, or transport **extensive** nodal masses $w_i\rho_i$ rather than
intensive samples, in which case the unweighted invariant of the extensive quantity *is* the
mass. The choice is a modelling decision that the current public interface does not make for
the user, and the experiment does not assume one.

### 5.3 The documented contract

The package's own statement, in the docstring of the conservative overload
(`src/Kinematic/TransportOperators/SemiLagrangian.jl`, `lagrangianTransport!`):

> "Apply conservative CCISL redistribution defined by `map`. For a valid map, this preserves the
> unweighted nodal sum over its source support (not a physical-volume integral); use the
> velocity-based `lagrangianTransport!` overload for intensive CIP."

Read literally, the implementation keeps its promise and explicitly disclaims the physical
integral. The experiment's contribution is not to catch a violated contract; it is to
*quantify* how far the two notions diverge on a concrete moving-domain problem, and to show
that the divergence is large enough ($2\%$ of the mass, versus $4\times10^{-4}$ for the
geometry) that a user could reasonably mistake the first property for the second. The
`TransportMap` docstring also warns that the operator is "conservative, not shape-preserving"
and gives measured smoothing figures; §4.5's timestep behaviour (fewer, larger steps are more
accurate in $L^{1}$) is consistent with that characterisation.

### 5.4 Four properties, tested separately

| property | definition | result |
|---|---|---|
| 1. Algebraic conservation | $\mathbf 1^{\mathsf T}q^{n+1}=\mathbf 1^{\mathsf T}q^{n}$ | holds, $\le4.44\times10^{-16}$ relative |
| 2. Constant-state preservation | $P\mathbf 1=\mathbf 1$ on the target support | fails at the support edge: $\max|P\mathbf 1-1|=0.86$; $0.039$ at one node in |
| 3. Geometry-weighted physical mass | $w^{\mathsf T}q$ with cut-cell weights | fails: $-2.0\%$ (64²), $-0.85\%$ (128²), against $-3.8\times10^{-4}$/$-6.6\times10^{-6}$ for the exact density on the same geometry |
| 4. Pointwise accuracy and boundedness | $\|q-\rho_{\mathrm{exact}}\|$; range constraints | interior $L^{1}$ improves with refinement; the innermost material layer does not; boundedness fails on void-side nodes ($\le1.93$ vs exact $\le1.19$) |

---

## 6. Analysis of the physical mass loss

### 6.1 Observations (measured)

* **O1.** The fine-grid mass drift is $-8.51\times10^{-3}$, against $-6.59\times10^{-6}$ for the
  exact density integrated over the *same* reconstructed triangles. The geometry is therefore
  an unlikely explanation for the larger part of the drift.
* **O2.** Mass drift improves with refinement by a factor $2.35$ (from $-2.00\times10^{-2}$ to
  $-8.51\times10^{-3}$) for a factor-two refinement in $h$. That is between first and second
  order.
* **O3.** Mass drift is nearly timestep independent ($-1.90,-2.00,-2.15\times10^{-2}$ for CFL
  0.2/0.4/0.8), while the $L^{1}$ density error improves as the number of steps decreases.
  Whatever removes mass is therefore not a per-step time-truncation effect in the usual sense.
* **O4.** The drift has a mean downward trend plus an oscillation of a few steps
  (figure 4, right). Its period is comparable to the time the annulus needs to advance one cell
  in $x$, $h/0.24=6.5\times10^{-2}$, i.e. $\approx3.6$ steps at $\Delta t=1.81\times10^{-2}$ —
  consistent with, but not proof of, a mechanism driven by nodes crossing the interface. The
  exact-density control shows no such oscillation at its own scale, so the oscillation belongs
  to the transported field and to the discrete change of the quadrature.
* **O5.** Replacing the void-side nodal values by their exact counterparts does *not* reduce
  the drift (it changes $-1.999\times10^{-2}\to-2.155\times10^{-2}$ and
  $-8.512\times10^{-3}\to-9.291\times10^{-3}$). Void-side pollution is not the cause.
* **O6.** The map's constant-state defect is $0.862$ on the support and $0.0396$ one node in,
  and is essentially unchanged by refinement ($0.8622$ / $0.8611$; $0.0396$ / $0.0394$). It is
  a property of the reconstruction at the truncation edge, not of $h$.
* **O7.** The density error is concentrated at the support edge (§4.3), and its innermost
  material bin degrades with refinement while the deeper bins improve.
* **O8.** The support contains noticeably more nodes than the material: $741$ vs $485$ at
  $64^{2}$, $2591$ vs $2070$ at $128^{2}$ ($t=0$).

### 6.2 Reasonable hypotheses (not established)

* **H1 — edge renormalisation.** The pull stencil is 4 nodes wide while the support band is
  about one cell wide, so rows near the edge are built from a truncated stencil and renormalised
  to unit row sum. Renormalisation preserves $\sum q$ but *moves* mass between nodes with
  different effective quadrature weights. Because the layer occupies a large fraction of a
  7.7-cell-wide annulus, a small per-node bias translates into a per-mille-level mass change
  per step. Consistent with O5–O7.
* **H2 — demand rescaling.** Columns whose demand exceeds one are scaled by $1/d_j$, and the
  remainder is re-injected by the forward leakage pass through a different stencil. Both are
  geometry-determined linear maps that conserve the column sum but not the integral.
  Consistent with O3 (a geometry effect, insensitive to $\Delta t$) and with the sawtooth in
  O4, since $d_j$ changes discontinuously when nodes enter or leave the support.
* **H3 — dissipation/quadrature interaction.** The piecewise-quadratic reconstruction with
  supersampled rays is dissipative; smoothing near the interface lowers nodal values whose
  quadrature weight is mostly inside the material, so the integral falls even though the nodal
  sum is fixed. Consistent with O2 (partial order) and with the CFL trend in O3.
* **H4 — genuinely distinct, competing effects.** The observations are compatible with H1–H3
  acting together at comparable magnitude; nothing in the present data separates them.

### 6.3 Unverified, and one confounded control

The following are **not** established by this campaign:

* the split of the drift between H1, H2 and H3, and its sign in general;
* whether the drift persists for a different velocity, a thicker annulus, or a different
  quadrature rule;
* whether an extensive (mass-carrying) nodal representation would conserve the integral with
  the same operator.

A specific qualification applies to the control used throughout. The "exact density on the
same cut" integral evaluates $\rho_{\mathrm{exact}}$ **analytically at the triangle
centroids**, whereas the numerical mass evaluates $\rho_h$ by **bilinear interpolation of nodal
samples** at those same centroids. The two differ by the interpolation error of the nodal
representation, which is precisely the quantity under investigation in part. That control
therefore bounds the *geometry and quadrature* contribution under an idealised field; it does
**not** isolate transport error alone, and it cannot be read as "the transport is responsible
for the remaining $8.5\times10^{-3}$" without the additional step below — which has not been
performed:

> Initialise exact density samples on the Cartesian grid, reconstruct them with the same
> `linearInterpolator`, and integrate with the identical cut-cell centroid quadrature. The
> difference between that integral and $M_{\mathrm{exact}}$ is the pure
> interpolation-plus-quadrature bias of the chosen representation at $t=0$; only what remains
> after removing it can be attributed to transport.

Campaign 2, Experiment D is designed for exactly this decomposition.

---

## 7. Assessment of architectural composability

### 7.1 Demonstrated

* **Geometry is inspectable and independently integrable.** `Cut` and `CutSurface` publish
  their arrays as owned byproducts "an external consumer can use directly (curvature, surface
  integrals, plotting) without touching `Cut` internals" (`src/CutData.jl`). This experiment
  used them for exactly that: areas, centroids and mass from `subcell_states` +
  `point_pointers`/`point_ids` + `physical_points`, and interface segments for the figures. No
  accessor, wrapper or `getfield` on an undocumented field was needed.
* **The algebraic invariant is directly measurable.** The map exposes `source_mask` and
  `target_mask`; the experiment recomputed the supports itself from `bg_states` and asserted
  equality on every step. That is a genuine, independent check of the operator's support
  semantics rather than a trust relationship. (This is the one place where a field of an
  exported struct that is not in the export list is read — used only in assertions, never in
  the computation.)
* **Multilevel data composes.** Three named levels from `refine!`, hierarchy navigation through
  `coarsest`/`finest`/`immediateCoarser`/`immediateFiner`, and `transfer!` in both directions
  produced a clean second-order prolongation error and an exact nodal restriction; the coarse
  and fine integrals were then compared on the same quadrature.
* **The loop stays explicit.** §3.3 is the whole numerics of the experiment; nothing is hidden
  in a solver abstraction, and the two non-invasive diagnostics (constant-state probe, exact-
  density control) reuse the same `TransportMap` object, which is only possible because the map
  is value-independent, as documented.
* **The numerical limitation was investigated without touching upstream source.** Every
  hypothesis in §6 was formulated and partially tested using public data only.

### 7.2 Not demonstrated

* **Interoperation with an independent solver.** No WENO5/CIP/DG transport of $\rho$ was run
  against the cut map in this campaign (the level set *is* advanced by `eulerianTransport!`,
  but that is a different quantity, not a conservation reference). Experiments C and D are the
  natural continuation.
* **General compatibility with external FEM/DG stacks.** `aggregateCutCells` (AgFEM roots) and
  the `Cut`-to-Gridap/GridapEmbedded route were *not* exercised here; nothing in this report
  should be read as evidence about them.
* **Physical conservation of the chosen density representation** — see §5 and §6.
* **Independence from native packaging.** `cut` requires the native CutCells library. This
  checkout publishes no `Artifacts.toml` entry, so `cut` cannot start until a library is
  provided:

  > "CutCells platform artifacts have not been published in Artifacts.toml. Either point
  > EVOLVINGDOMAINS_CUTCELLS_LIBRARY at an existing library, or build one locally …"

  The experiment built it locally with the repository's own `scripts/build_cxx_bindings.sh`
  (cmake fetches CutCells at a pinned commit; the build writes to a temporary prefix) and ships
  the resulting 132 KB library in `lib/`, which `run.ConservativeAnnulus.jl` points the
  environment variable at when the caller has not. This is a **packaging and deployment**
  concern — reproduction requires a compiler, cmake, network access to the pinned CutCells
  source, and platform-specific binaries — not evidence of an API design failure: the Julia
  surface of `cut` was stable and sufficient once the library was present, and the failure mode
  was a clear, actionable error rather than a silent degradation.

One experiment on one problem is not evidence of universal composability. What it does support
is narrower and still valuable: the geometric and multilevel data are transparent enough that
an external script can reconstruct the discrete operators' inputs, measure its invariants, and
diagnose its accuracy without forking the package.

---

## 8. Proposed experiments for the second campaign

None of the following has been performed. Each states the question and the evidence that would
decide it.

**Experiment A — Weighted conservation.**
*Question:* is there a choice of geometry-dependent weights $w$ such that
$w^{\mathsf T}P=\tilde w^{\mathsf T}$ approximately holds, and does the public interface intend
intensive density or extensive nodal mass?
*Method:* build $w$ from the cut (e.g. the in-phase fraction of each node's dual cell, or the
`CUT_IN` subcell areas attributed to nodes), apply the same maps, and measure
$w^{\mathsf T}q^{n+1}/\tilde w^{\mathsf T}q^{n}-1$; compare with the uniform-weight and the
extensive ($w\odot q$) formulations.
*Decisive evidence:* a weight family giving drift at the $10^{-12}$ level for all steps would
confirm that the operator is weighted-conservative for a suitable discrete measure; a residual
$O(h)$ drift independent of the weights would falsify the hypothesis that the map merely needs
the right measure, and point back to H1/H2.

**Experiment B — Constant-state preservation versus support size.**
*Question:* is the constant-state defect fully explained by $|\mathcal S|\ne|\mathcal T|$ and by
stencil truncation at the support edge?
*Method:* initialise $\rho\equiv1$ on the support, remap, and correlate
$\max|q^{n+1}-1|$ with $\Delta|\mathcal S|=|\mathcal T|-|\mathcal S|$, with the number of
truncated rows, and with the depth of each node in the band (binned by $\phi/h$).
*Decisive evidence:* a defect that vanishes where $\mathcal S=\mathcal T$ and scales with the
count of truncated rows confirms the §5.1 algebra as the operative mechanism (and confirms that
the observed drift is a geometric bookkeeping effect, not an implementation error). A defect
that persists with identical supports would falsify it and indicate a bug.

**Experiment C — Independent transport comparison.**
*Question:* how much of the observed mass and density error belongs to the conservative map
specifically?
*Method:* transport $\rho$ on the same geometry with `eulerianTransport!` (WENO5) and with CIP
(the intensive overload), driving the level set *identically* in all runs; evaluate physical
mass with the same cut-cell quadrature.
*Decisive evidence:* if WENO5/CIP reproduce mass to much better than $10^{-3}$ while the
conservative map loses $10^{-2}$, the map's edge handling is implicated; if they lose
comparable mass, the representation (nodal sampling + cut quadrature) is implicated instead.

**Experiment D — Controlled spatial error decomposition.**
*Question:* which term dominates the budget — level-set geometry, nodal interpolation,
cut-cell quadrature, or transport?
*Method:* four integrals of the same analytic density, all on the same cut: (i) analytic $\rho$
at centroids (already available); (ii) analytic nodal samples reconstructed with
`linearInterpolator` at the centroids, integrated identically (the control §6.3 says is
missing); (iii) the transported field, same quadrature; (iv) a high-order reference on the
exact geometry. Differences (i)→(ii)→(iii) separate interpolation+quadrature bias from
transport error; a comparison of (i) with (iv) isolates the geometry.
*Decisive evidence:* this yields the missing attribution directly, and any later claim about
"transport mass loss" must be expressed as the (ii)→(iii) difference.

**Experiment E — Resolution and support-width study.**
*Question:* do the boundary errors scale with $h$, with the band width in cells, or with the
reconstruction stencil size?
*Method:* fix the geometry and vary (a) $h$ at fixed annulus width (grids 65², 129², 257²,
513²); (b) the annulus width at fixed $h$ (0.06, 0.12, 0.24 wide); (c) the CFL. Record the
mass drift, the constant-state defect and the layer-resolved $L^{1}$ error in each case.
*Decisive evidence:* if the mass drift per unit mass scales with $1/(\text{width in cells})$
at fixed $h$ and with $h$ at fixed width, the layer hypothesis (H1) is quantified; if the deep
interior error converges at second order in every case while the layer error does not, then the
operator is interior-accurate and edge-limited, which is a much sharper statement than "the
scheme is diffusive".

A useful by-product of A–E would be a decision about what the *package* should expose: a
documented node-weight constructor (or a `restrict!`-style conservative restriction), and an
explicit statement in the `TransportMap` docstring of the edge-layer behaviour with its
measured magnitude, would let users distinguish §5's four properties without re-deriving them.

---

## 9. Conclusions

The geometric half of the experiment succeeded clearly. A translating annulus cut from a
level-set field is reconstructed with relative area error $3.8\times10^{-4}$ at 64² cells and
$7.9\times10^{-6}$ at 128² cells, centre error at or below $2\times10^{-5}$, and an interface
that is visually indistinguishable from the analytic circles. The multilevel machinery behaves
as documented: nested refinement, hierarchy navigation, second-order prolongation (observed
order 1.99) and exact nodal restriction.

The algebraic claim also holds: the conservative map preserves the unweighted nodal sum over
its support to a few units in the last place, at every step, independent of the timestep, on
both grids.

The physical claim does not follow from it. Physical mass, measured by independent quadrature
over the cut subcells, decreases by $2.0\%$ at 64² cells and $0.85\%$ at 128² cells, against
$3.8\times10^{-4}$ and $6.6\times10^{-6}$ for the exact density on the same geometry. The
density error is concentrated in the two-node edge layer of the transported support, where the
map's partition-of-unity defect is $0.86$ — the same layer whose existence the §5.1 algebra
predicts once source and target supports differ in size. Refinement improves the mass drift
(by a factor 2.35) and the deep-interior density error, but not the maximum pointwise error at
the interface-adjacent nodes.

On the software side, the experiment supports the package's minimalist and composable design
goals within its scope: geometry, supports, invariants and multilevel data were all accessed
through exported names and documented structures, the loop remained explicit, and the
diagnostics — including the ones that show the operator in a less favourable light — required
no modification of, or dependence on, package internals. What remains unexercised is
interoperation with an independent solver, general compatibility with external FEM/DG stacks,
and independence from native-binary packaging.

The result is therefore neither a success nor a failure, but a precise boundary: **the operator
behaves as documented, and the documentation's notion of conservation is not the notion most
users will bring to a moving-domain transport problem.** The origin of the mass drift has been
localised to the support-edge layer and quantified, but not conclusively attributed; the
controlled decomposition of Campaign 2 (Experiments A–E) is required before drawing further
conclusions.

---

## Appendix — reproduction and artifact inventory

```sh
julia validation/ConservativeAnnulus/run.ConservativeAnnulus.jl        # ≈10 s, writes output/results.jls
julia validation/ConservativeAnnulus/visualize.ConservativeAnnulus.jl # writes the five PNGs
```

Both scripts activate `validation/env` (Julia 1.11.7; EvolvingDomains from this checkout via
`[sources]`, CairoMakie 0.15.17, StaticArrays, Trixi present but unused).
`EVOLVINGDOMAINS_CUTCELLS_LIBRARY` is set by the run script to `lib/` unless the caller
already provided a library; see §7.2 for how `lib/` was produced.

| Path | Contents |
|---|---|
| `run.ConservativeAnnulus.jl` | model, `transferStudy`, `simulate`, diagnostics, report, entry point (584 lines) |
| `visualize.ConservativeAnnulus.jl` | five CairoMakie figures (225 lines) |
| `output/results.jls` | per-step records, snapshots, transfer study, sensitivity runs (Julia `Serialization`) |
| `output/1_cut_geometry.png` | figure 1 |
| `output/2_density.png` | figure 2 |
| `output/3_geometry_errors.png` | figure 3 |
| `output/4_conservation.png` | figure 4 |
| `output/5_transfer.png` | figure 5 |
| `lib/libevolvingdomains_cut_cells.dylib` | locally built native cut backend (132 KB, macOS arm64) |
| `README.md` | short campaign summary |

No EvolvingDomains source file, dependency manifest or project configuration was modified in
this campaign. The archive of the run summarized in §4 is the one stored in
`output/results.jls`; the report does not rerun or extend it.
