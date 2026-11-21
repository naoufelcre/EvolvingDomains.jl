# PlanarStefan — Scientific and Architectural Report

**Case:** one-phase Stefan melting in a two-dimensional rectangle, moving interface represented by a level set, liquid temperature computed by an independent embedded (cut) finite-element solver.

**Components under test:** EvolvingDomains.jl (moving geometry and kinematics), Gridap 0.20.10 and GridapEmbedded 0.9.12 (thermal finite-element discretization and cut geometry), Julia sparse linear algebra (implicit time-step solve), application code (coupling, flux reconstruction, Stefan dynamics).

**Deliverables in this directory:** `run.PlanarStefan.jl` (606 lines), `visualize.PlanarStefan.jl` (150 lines), `Project.toml` / `Manifest.toml` for a self-contained environment, `build/libevolvingdomains_cut_cells.dylib` (locally built native CutCells backend), `output/` (numerical outputs, figures, convergence study), and this report.

**Executed environment:** Julia 1.11.7, EvolvingDomains 1.0.0 (path `../..`), Gridap 0.20.10, GridapEmbedded 0.9.12, CairoMakie 0.15.17, SpecialFunctions 2.9.2, StaticArrays 1.9.24, macOS (arm64). All numerical claims below are taken from `output/summary.txt`, `output/history.csv`, `output/profiles.csv`, `output/snapshots.csv`, `output/operators.csv` and `output/convergence.txt` of this directory.

---

## 3.1 Executive summary

A one-phase Stefan melting problem was solved from $t_0=0.05$ to $t_f=0.20$ on a $65^2$-node Cartesian grid with 75 implicit time steps of size $\Delta t=0.002$. The moving liquid domain was represented by an EvolvingDomains `Field` holding the level set $\phi=x-s(t)$, advanced with the package's WENO5/SSP-RK3 `eulerianTransport!`; the same nodal level-set array was handed to GridapEmbedded, which constructed an independent cut of the background mesh; Gridap assembled the mass and stiffness operators of $\rho c_p\partial_tT=k\Delta T$ on the cut liquid domain, with weak (symmetric Nitsche) enforcement of both Dirichlet boundaries and a ghost penalty on the cut-cell skeleton. The partitioned time loop, the interface heat-flux reconstruction, the Stefan velocity and the geometry update are explicit application code; the thermal operators exist as ordinary sparse matrices ($M$, $K$ in `output/operators.csv`).

Against the analytical similarity solution, the final numerical interface position is $s_h(t_f)=0.394739$ versus $s_{\text{exact}}(t_f)=0.392162$, a relative error of $+0.657\,\%$ ($+2.58\times10^{-3}$). The liquid-temperature RMS error at the final state is $4.3\times10^{-4}$ (maximum nodal error over the run $1.13\times10^{-2}$, occurring at the first step), the mean relative interface heat-flux error is $1.35\,\%$ (decaying to $0.40\,\%$ at $t_f$), and the linear systems solve to a relative residual of $3.4\times10^{-16}$. The temperature field remains numerically one-dimensional ($y$-variation $\le 3.7\times10^{-15}$), and the two independent cuttings of the same level-set array classify all 4096 background cells identically. Accuracy is sensitive to the Nitsche penalty (relative interface error $0.93\,\%$ at $\gamma_d=10$ versus $0.65\,\%$ at $\gamma_d=300$), which points at weak Dirichlet enforcement and the gradient-based flux recovery as the accuracy-limiting mechanisms rather than at the algebraic solve.

This is strong evidence of cross-package interoperability on a nontrivial coupled free-boundary problem: two unrelated geometric constructions agreed cell-by-cell on data exchanged only as arrays of values and coordinates, and the PDE discretization stayed accessible as explicit matrices. It is **not** proof that the moving-domain PDE problem class is validated: the interface stayed exactly planar and axis-aligned, one finite-element order and two mesh resolutions were used, the interface flux is recovered from a background-mesh gradient, the saved visualization snapshots are staggered by one time step, and possible cut-cell/polygonal "apparent convergence" of the flux recovery is uncontrolled. Those gaps define the second campaign.

---

## 3.2 Mathematical formulation

The computational box is $D=[0,1]^2$. The liquid occupies

$$
\Omega(t)=\{(x,y)\in D:\ 0<x<s(t)\},
$$

and the solid occupies $\{x>s(t)\}$. The interface is $\Gamma(t)=\{x=s(t)\}$, with outward (liquid-to-solid) normal $\mathbf n=(1,0)$.

Within the liquid the temperature satisfies the heat equation

$$
\rho c_p\,\partial_t T = k\,\Delta T
\qquad \text{in } \Omega(t),
$$

subject to

$$
T(0,y,t)=T_h,\qquad
T(s(t),y,t)=T_m,\qquad
\partial_y T=0 \ \text{at } y=0,1 .
$$

The interface motion is the one-phase Stefan condition

$$
\rho L\,\frac{ds}{dt}
=
-k\,\frac{\partial T}{\partial x}\bigl(s(t)^{-},t\bigr),
$$

with the positive $x$ direction being the outward normal of the liquid. The solid phase is assumed isothermal at $T_m$ (one-phase model); no equation is solved there.

**Parameters used** (all values read from `output/summary.txt`):

| symbol | value |
|---|---|
| $\rho$ | 1 |
| $c_p$ | 1 |
| $k$ | 0.5 |
| $L$ | 1 |
| $T_h$ | 1 |
| $T_m$ | 0 |

Consequently

$$
\alpha=\frac{k}{\rho c_p}=0.5,
\qquad
\mathrm{Ste}=\frac{c_p(T_h-T_m)}{L}=1 .
$$

### Analytical reference

The classical similarity solution of this problem is

$$
s_{\text{exact}}(t)=2\lambda\sqrt{\alpha t},
\qquad
T_{\text{exact}}(x,t)
=
T_h-(T_h-T_m)\frac{\operatorname{erf}\!\left(\dfrac{x}{2\sqrt{\alpha t}}\right)}{\operatorname{erf}(\lambda)},
$$

where $\lambda>0$ solves

$$
\mathrm{Ste}=\sqrt{\pi}\,\lambda\,e^{\lambda^2}\operatorname{erf}(\lambda).
$$

The scalar equation is solved once at initialization by bisection on $[0,2]$ (`similarityRoot`), giving (recorded)

$$
\lambda = 0.620062633314,
\qquad
\sqrt{\pi}\lambda e^{\lambda^2}\operatorname{erf}(\lambda)-\mathrm{Ste}
= 6.66\times10^{-16},
\qquad
s(t_0)=0.196081,\quad s(t_f)=0.392162 .
$$

**Why the simulation starts at $t_0=0.05$.** Both $T_{\text{exact}}$ and $s_{\text{exact}}$ contain $\sqrt{t}$ and $\sqrt{\alpha t}$ in denominators; $T_{\text{exact}}$ is singular as $t\to0^{+}$, and the similarity length scale is $2\sqrt{\alpha t}$, i.e. $0.316$ at $t_0$ but zero at $t=0$. At $t_0$ that length scale resolves on $h=1/64=0.015625$ ($\approx 20$ cells), whereas at $t=0$ the boundary layer and the interface are unresolved. Starting at $t_0$ therefore tests the discretization rather than the initial singularity, at the price that no conclusion can be drawn about start-up behaviour.

### Role of the analytical solution

The analytical solution is used in three strictly diagnostic ways:

1. **Initialization** at $t_0$: $\phi(x,y,t_0)=x-s(t_0)$ and $T(x,y,t_0)=T_{\text{exact}}(x,t_0)$ for $x<s(t_0)$, $T_m$ elsewhere. This is the only place where analytical data enter the state.
2. **Reference values** for interface position, temperature, heat flux and Nitsche penalty.
3. **Nothing else.** During the loop the interface is advanced only by `eulerianTransport!` applied to the numerical velocity field; the temperature is advanced only by the finite-element implicit step on the numerical cut geometry. The analytical interface is never written into the level set, and the analytical temperature is never written into the temperature field after initialization. This is verifiable in the loop excerpt of §3.3: the sole writes to `state.level` and `state.T` are `advanceGeometry!` and the solution of $(M+\Delta tK)T^{n+1}=MT^n+\Delta t\,r$.

---

## 3.3 Software architecture and component decomposition

### Roles

| Layer | Package | Responsibility in this case |
|---|---|---|
| Data model | EvolvingDomains | `Grid` (uniform Cartesian axes, mesh identity), `Discretization` (level $\to$ grid), `Field` (dense nodal arrays for $\phi$ and $T$); native cut-cell snapshot of $\phi$ used as an independent geometric reference; `eulerianTransport!` (WENO5 + SSP-RK3) for the level-set update |
| Cut geometry | GridapEmbedded | `DiscreteGeometry(nodal values, node coordinates)`; `cut(background model, geometry)`; `Triangulation(cutgeo, PHYSICAL)` (liquid), `EmbeddedBoundary(cutgeo)` (interface), `GhostSkeleton(cutgeo)` (cut-cell facets), `compute_bgcell_to_inoutcut` (per-cell IN/CUT/OUT state) |
| Finite elements | Gridap | Q1 background space on the fixed mesh, `Measure` on the cut triangulations and on the wall facets, `assemble_matrix` / `assemble_vector`, symmetric Nitsche terms, sparse `SparseMatrixCSC` operators |
| Linear algebra | Julia `LinearAlgebra`/`SparseArrays` | LU factorisation and solve of the implicit step; residual and condition-number diagnostics |
| Application | `run.PlanarStefan.jl` | flux reconstruction, Stefan velocity, velocity extension, geometry update, time loop, diagnostics and output |
| Presentation | `visualize.PlanarStefan.jl` (CairoMakie) | reads only CSV files; no PDE solve, no field update |

### The adapter: an exchange of arrays, not of solver objects

The complete geometric coupling is three statements (`buildGeometry`):

```julia
function buildGeometry(state, model, X, nodeindex)
    φnode = [state.level.data[:init][nodeindex[q]...] for q in eachindex(X)]
    geo = DiscreteGeometry(φnode, X)
    cutgeo = cut(model, geo)
    return (geo=geo, cutgeo=cutgeo,
        bgstate=compute_bgcell_to_inoutcut(cutgeo, geo),
        Ω=Triangulation(cutgeo, PHYSICAL),
        Γ=EmbeddedBoundary(cutgeo),
        S=GhostSkeleton(cutgeo))
end
```

`φnode` is a plain `Vector{Float64}` of level-set samples in Gridap node order; `X` is the list of node coordinates. GridapEmbedded is told nothing about EvolvingDomains: it receives values and coordinates and builds its own geometry. The reverse direction is equally shallow — the solved dof vector is scattered into the EvolvingDomains `Field` array (`state.T.data[:init]`); nothing else crosses the interface.

### Node and dof ordering: the only real adapter work

Two consistency conditions must be established and are checked at startup with hard errors rather than assumptions.

*(a) Grid node coordinates ↔ EvolvingDomains grid indices* (`nodeIndexMap`):

```julia
map = [(round(Int, pt[1] / h) + 1, round(Int, pt[2] / h) + 1) for pt in X]
dev = maximum(max(abs(X[q][1] - state.x[map[q][1]]), abs(X[q][2] - state.x[map[q][2]]))
              for q in eachindex(X))
dev == 0.0 || error("Gridap nodes and the EvolvingDomains grid disagree by $dev")
```

The recorded deviation is exactly $0.0$: the two grids have identical node coordinates, including ordering conventions, so the mapping is an exact index bijection (no interpolation is needed anywhere).

*(b) Finte-element dof ↔ background node* (`dofNodeMap`): because the Q1 space is built on the *background* model, dofs are background nodes, and the mapping is recovered from `get_cell_dof_ids(V)` and `get_cell_node_ids(model)` and then verified against `get_free_dof_coordinates(V)`:

```julia
dev = maximum(max(abs(Xg[d][1] - state.x[nodeindex[dofnode[d]][1]]),
    abs(Xg[d][2] - state.x[nodeindex[dofnode[d]][2]])) for d in eachindex(Xg))
dev == 0.0 || error("dof/node mismatch of $dev")
```

Recorded deviation: $0.0$. The permutation is computed once (it cannot change: the background mesh is fixed), which is what makes the temperature transfer in §3.4 a pure copy.

### The orchestration loop (verbatim structure, comments abridged)

```julia
for step in 0:(nsteps - 1)
    # 1. current geometry and cut
    geom = buildGeometry(state, model, X, nodeindex)
    snapshot = ED.cut(state.level, :init)                 # native cut, independent reference
    cmp = compareCuttings(state, model, geom, snapshot)
    s_area = sum(get_cell_measure(geom.Ω))                # liquid area = interface position

    # 2./3. temperature on Ω(t);  interface flux;  4. Stefan speed
    therm = solveTemperature(state, p, geom, dis, dis.dofnode, nodeindex; dt=p.dt)
    q, ifacearea, nmean = interfaceHeatFlux(therm, geom, p, dis)
    Vstefan = dynamics(q, p)                              # V = -q / (rho L)

    ... diagnostics ...

    # 5. extend the speed;  6. advance the level set;  7. advance time
    velocity = velocityFromRates(Vstefan, state)          # (V, 0) at every node
    advanceGeometry!(state, velocity, p.dt)               # ED.eulerianTransport! (WENO5+SSP-RK3)
    prev = (area=s_area, line=s_line, q=q, V=Vstefan)
    t = tT
end
```

with

```julia
interfaceHeatFlux(therm, geom, p, dis) = begin
    uh = FEFunction(dis.U, therm.u)
    nΓ = get_normal_vector(geom.Γ)
    area = sum(∫(1.0) * therm.ops.dΓ)
    q = sum(∫(p.k * (∇(uh) ⋅ nΓ)) * therm.ops.dΓ) / area
    (q, area, nmean)
end

dynamics(q, p) = -q / (p.ρ * p.L)
velocityFromRates(V, state) = vec([SVector(V, 0.0) for _ in state.x, _ in state.x])
advanceGeometry!(state, velocity, dt) =
    ED.eulerianTransport!(state.level, velocity, dt; level_id=:init, cache=ED.WENO5Cache())
```

The loop is the application's own; no framework type hides the sequence, and each of the required application-level operations exists as a named function with a mathematical signature.

---

## 3.4 Independent finite-element formulation

### Spaces and operators

The background model is the fixed tensor mesh of $(n-1)\times(n-1)$ quads, $n=65$. The finite-element space is $Q_1$ (order 1) on that background mesh, i.e. trilinear/bilinear hats with the *background nodes* as dofs (4225 dofs). Integrals are restricted to the physical (liquid) part of the cut mesh, so active dofs are exactly the nodes touched by cells that intersect the liquid.

The heat-capacity mass matrix and the bulk conductivity stiffness matrix are

$$
M_{ij}=\int_{\Omega(t)}\rho c_p\,\varphi_i\varphi_j\,dx,
\qquad
K^{\text{bulk}}_{ij}=\int_{\Omega(t)}k\,\nabla\varphi_i\cdot\nabla\varphi_j\,dx .
$$

Both Dirichlet data are imposed **weakly** by symmetric Nitsche, with $h$ the background mesh size, $\mathbf n_\Gamma$ the outward normal of the liquid on the embedded boundary and $\mathbf n_\Lambda$ the outward normal of the heated wall:

$$
K_{ij}=K^{\text{bulk}}_{ij}
+\int_{\Gamma(t)}\Bigl[-k(\mathbf n_\Gamma\!\cdot\!\nabla\varphi_i)\varphi_j
-k(\mathbf n_\Gamma\!\cdot\!\nabla\varphi_j)\varphi_i
+\frac{\gamma_d k}{h}\varphi_i\varphi_j\Bigr]ds
+\int_{\Lambda}\Bigl[-k(\mathbf n_\Lambda\!\cdot\!\nabla\varphi_i)\varphi_j
-k(\mathbf n_\Lambda\!\cdot\!\nabla\varphi_j)\varphi_i
+\frac{\gamma_d k}{h}\varphi_i\varphi_j\Bigr]ds
+\int_{\mathcal S}\gamma_g k h\,\bigl[\!\bigl[\mathbf n_{\mathcal S}\!\cdot\!\nabla\varphi_i\bigr]\!\bigr]
\bigl[\!\bigl[\mathbf n_{\mathcal S}\!\cdot\!\nabla\varphi_j\bigr]\!\bigr]\,ds,
$$

$$
r_i=\int_{\Lambda}\Bigl[\frac{\gamma_d k}{h}T_h\,\varphi_i-k(\mathbf n_\Lambda\!\cdot\!\nabla\varphi_i)T_h\Bigr]ds,
$$

with the interface data $T_m=0$ (hence no interface contribution to $r$). $\Lambda$ is the $x=0$ facet set (64 facets) selected by a coordinate mask `wallFacetMask`, and $\mathcal S$ is `GhostSkeleton(cutgeo)`. The Nitsche form is the standard consistent, symmetric and (for $\gamma_d$ large enough) coercive one; consistency follows from $-\Delta T=0$ in the liquid and integration by parts.

**Stabilization parameters actually used** (recorded in `summary.txt`):

$$
\gamma_d=100,\qquad \gamma_g=0.1 .
$$

The ghost term penalizes the jump of the normal derivative across background facets adjacent to at least one cut cell; it is what makes the equations of dofs living on the far side of the interface well-posed and keeps arbitrarily small cut cells from destroying the conditioning of the implicit operator.

### Time step

$$
\boxed{\;(M+\Delta t\,K)\,T^{n+1}=M\,T^{n}+\Delta t\,r\;}
\qquad \text{(implicit Euler, } \Delta t=0.002\text{)} .
$$

Here $T^n$ and $T^{n+1}$ are the dof vectors restricted to the *same* geometry $\Omega(t^n)$: the thermal problem is solved with the interface frozen at the beginning-of-step position.

### Algebraic versus physical objects

Two distinct things must not be conflated:

* The **physical operators** are the bilinear forms above: integrals over the current liquid domain of mass, conductivity, Nitsche and ghost contributions. They change every step because $\Omega(t)$, $\Gamma(t)$ and $\mathcal S$ change.
* Their **sparse representations** are `SparseMatrixCSC` objects returned by `assemble_matrix` (`ops.M`, `ops.K`). They are ordinary algebraic objects: `output/operators.csv` contains their upper-triangle triplets `(matrix,dof_i,dof_j,value)` for the last step (607 kB), so a reviewer can re-check symmetry, the partition of unity, or re-assemble the solve independently.

Solving the assembled system correctly says nothing about whether the assembled integrals approximate the PDE well; §3.6 keeps the two apart.

### State transfer between successive active domains

Because the background mesh is fixed and the Q1 dofs *are* the background nodes, no geometric interpolation is required. Concretely (`solveTemperature`):

```julia
uold = [Tnodal[nodeindex[dofnode[d]]...] for d in eachindex(dofnode)]   # previous state
b    = ops.M * uold + dt * ops.r                                       # mass-matrix transfer
...
active    = findall(>(1e-12 * maximum(diagA)), diagA)                  # dofs with a non-empty equation
breg[inactmask] .= p.Tm                                                # solid nodes: T = Tm
Areg = A + spdiagm(0 => Float64.(inactmask))                           # identity row for those dofs
u = Areg \ breg
```

Three points:

1. A dof whose cells contain no liquid has an empty equation in $M+\Delta tK$. Rather than leaving it singular, its row is replaced by the identity with right-hand side $T_m$. This is the physically consistent value, not padding: the one-phase model states that the solid is isothermal at $T_m$, which is also the interface value.
2. A node that the liquid gains in step $n+1$ was, at the end of step $n$, one of those inactive dofs and therefore held exactly $T_m$. Its initialization is thus *identical* to the solid state, and the recorded diagnostic confirms it: over the 845 dof activations of the run (max 65 per step) the initialization deviation $\max|T_{\text{init}}-T_m|$ is exactly $0.0$.
3. Because the mesh is fixed, the transfer is the identity on shared dofs. Nothing is lost, and no interpolation error enters. This is a genuine simplification of the planar case, not a general capability.

**Limitations of this treatment.** (i) It relies on a fixed background mesh: on a moving or adaptively refined mesh, or with a Lagrangian/ALE geometry, the transfer would need a real projection (and would introduce its own error, e.g. the interpolation error documented for multilevel `transfer!`). (ii) It relies on Q1-isoparametric dof = node identity; for $Q_2$ or higher-order spaces the dofs are not grid nodes and the scatter above would have to be replaced by pointwise evaluation or an $L^2$ projection onto the EvolvingDomains samples. (iii) It assumes an exactly one-phase problem with an isothermal solid; two-phase or temperature-dependent solid models would make the inactive-dof value $T_m$ a modelling error rather than an assumption.

---

## 3.5 Temporal coupling and splitting

### Executed sequence

Per step $n$ (from $t^n$ to $t^{n+1}=t^n+\Delta t$):

1. build the cut geometry from the current level set $\phi^n$ (`buildGeometry`), and take the native EvolvingDomains `Cut` snapshot for the cross-check;
2. assemble $M$, $K$, $r$ on $\Omega(t^n)$;
3. advance the temperature implicitly, $T^{n+1}=$ solution of $(M+\Delta tK)T^{n+1}=MT^n+\Delta t\,r$, which is $T$ at time $t^{n+1}$ on the *frozen* domain $\Omega(t^n)$;
4. reconstruct the interface heat flux $q_h=\langle k\nabla T^{n+1}\cdot\mathbf n_\Gamma\rangle_\Gamma$;
5. apply the Stefan condition, $V_n=-q_h/(\rho L)$;
6. extend $V_n$ to the whole grid as $\mathbf u=(V_n,0)$;
7. advance the level set, $\phi^{n+1}=\text{WENO5}(\phi^n;\mathbf u,\Delta t)$, hence $\Omega(t^{n+1})$;
8. advance time, $t\mapsto t+\Delta t$.

### Accuracy consequences

* **First-order partitioned coupling.** The temperature is advanced on the old geometry and the geometry is advanced with the flux from the *new* temperature, but the two are never solved together; the interface position used by the thermal solve lags the interface position known by the geometry by exactly one step, $\Delta s=V_n\Delta t$ (derived from recorded values: $0.0044$ at $t=0.05$, $0.0020$ at $t=0.20$, i.e. $0.28$ to $0.13$ cells). Formally the splitting is $O(\Delta t)$.
* **First-order in time inside the thermal solve.** Implicit Euler is first order; its damping of the resolved frequencies exceeds that of the exact evolution, and this matters most for the sharp similarity profile (curvature concentrated in a boundary layer of width $2\sqrt{\alpha t}\approx0.32$).
* **Behaviour near the imposed initial profile.** Because $T(\cdot,t_0)$ is the *exact* solution and the geometry is frozen during the first thermal step, $T(\cdot,t_0+\Delta t)$ on $\Omega(t_0)$ necessarily differs from $T_{\text{exact}}(\cdot,t_0+\Delta t)$ on $\Omega(t_0+\Delta t)$ near the interface by roughly $|\partial_xT|\,\Delta s\approx4\cdot0.0044\approx0.018$. The measured consequence is exactly of this size: the maximum nodal temperature error over the whole run is $1.13\times10^{-2}$ and occurs at the first step, and the interface flux error peaks at $12.65\,\%$ at the first step, decaying to $0.40\,\%$ by $t_f$. This transient is a property of the initialization + splitting, not of the steady quality of the cut-FEM operators.
* The recorded refined-in-time runs quantify the total temporal contribution: at $65^2$ the interface error falls $0.657\%\to0.577\%\to0.539\%$ for $\Delta t=0.002,\,0.001,\,0.0005$ (§3.8), i.e. the temporal part at $\Delta t=0.002$ is a few tenths of a percent and the remainder is spatial.

### Saved-snapshot time-staggering (verified issue)

In `simulate`, snapshots and profiles are recorded **before** the geometry update:

```julia
        (step % snap_every == 0 || step == nsteps - 1) && recordSnapshot!(tT)
        (step % snap_every == 0 || step == nsteps - 1) && recordProfile!(tT)

        # 5. extend the speed;  6. advance the level set;  7. the state moves with the dofs
        velocity = velocityFromRates(Vstefan, state)
        advanceGeometry!(state, velocity, p.dt)
```

with `tT = t + p.dt`. Hence every saved pair is labelled $t+\Delta t$ while, inside that pair, the temperature `state.T` is the completed implicit step (time $t+\Delta t$) but the level set `state.level` is still the *preceding* geometric state (time $t$). The saved fields of one record are therefore out of phase by one step.

Measured effect in `output/snapshots.csv` (zero crossing of the saved $\phi$ at mid-height, compared with the analytical interface at the record's label):

| label $t$ | $\phi$-interface in snapshot | $s_{\text{exact}}(t)$ | difference | in cells ($h=0.015625$) |
|---:|---:|---:|---:|---:|
| 0.050 (initial record) | 0.196081 | 0.196081 | $-1.3\times10^{-13}$ | $-0.00$ |
| 0.052 | 0.196081 | 0.199964 | $-3.88\times10^{-3}$ | $-0.25$ |
| 0.072 | 0.233465 | 0.235297 | $-1.83\times10^{-3}$ | $-0.12$ |
| 0.112 | 0.292885 | 0.293467 | $-5.8\times10^{-4}$ | $-0.04$ |
| 0.152 | 0.341969 | 0.341879 | $+9.0\times10^{-5}$ | $+0.01$ |
| 0.192 | 0.384760 | 0.384239 | $+5.2\times10^{-4}$ | $+0.03$ |
| 0.200 | 0.392773 | 0.392162 | $+6.1\times10^{-4}$ | $+0.04$ |

The initial record is consistent (both fields at $t_0$). For the later records the drawn interface lags the drawn temperature by one step; the lag competes with the numerical interface error (which is positive and growing), so the *net* offset of the contour in the figures passes through zero around $t\approx0.14$ rather than showing a constant shift.

**Distinct from the reported final interface position.** After the loop, `simulate` builds the geometry from the final level set and reports $s_h(t_f)=0.394739$ from `get_cell_measure` of that geometry. That number is *not* the interface visible in the snapshot labelled $0.200$ (whose $\phi$ is at $t=0.198$ and yields $0.392773$). The $\Delta t$-consistent, post-update value is the one used for the headline error ($+2.577\times10^{-3}$); the snapshot is one step behind. The code is left unmodified here; synchronizing saved states is Recommendation A in §3.13.

---

## 3.6 Numerical results

Default run: $n=65$ nodes ($h=1/64$), $Q_1$, $\Delta t=0.002$, 75 steps, $\gamma_d=100$, $\gamma_g=0.1$, from `output/summary.txt`.

| Quantity | Measured value |
|---|---:|
| Grid | $65^2$ nodes (4225 dofs, 4096 cells) |
| Time step / steps | $0.002$ / $75$ |
| Total run time | 33.56 s (447.5 ms/step) |
| Similarity root $\lambda$ | 0.620062633314 |
| Final numerical interface $s_h(t_f)$ | 0.394739 |
| Exact interface $s_{\text{exact}}(t_f)$ | 0.392162 |
| Absolute interface error | $+2.577\times10^{-3}$ |
| Relative interface error | $+0.657\,\%$ |
| Max $\lvert s$ error $\rvert$ over run | $2.577\times10^{-3}$ |
| Final temperature RMS error (liquid nodes) | $4.291\times10^{-4}$ |
| Max nodal temperature error over run | $1.134\times10^{-2}$ (first step) |
| Max RMS temperature error over run | $5.194\times10^{-3}$ (at $t=0.060$) |
| Mean relative interface-flux error | $1.349\,\%$ |
| Max relative interface-flux error | $12.651\,\%$ (first step) |
| Relative flux error at last step ($t=0.198$) | $-0.404\,\%$ |
| Max relative linear-system residual | $3.366\times10^{-16}$ |
| Max temperature variation in $y$ | $3.664\times10^{-15}$ |
| Cut-cell classification mismatches (vs EvolvingDomains `Cut`) | $0$ of $4096$ |
| Stefan residual $\rho LV_n+k\partial_nT$ | mean $1.55\times10^{-7}$, max $7.68\times10^{-7}$ |
| Measured vs flux-predicted interface motion | max $1.54\times10^{-9}$ |
| Active dofs / total | 910 → 1755 / 4225 |
| Smallest cut-cell volume fraction | $5.168\times10^{-5}$ |
| Condition estimate of the active system | $3.511\times10^{2}$ (at $t=0.050$, 910 dofs) |
| Mass identity $\lvert \mathbf 1^\top M\mathbf 1/(\rho c_p\,|\Omega|) - 1\rvert$ | $1.110\times10^{-15}$ |
| Symmetry $\lvert A-A^\top\rvert_{\max}/\lvert A\rvert_{\max}$ | $M$: $0$; $K$: $1.07\times10^{-16}$ |
| Interface area $\int_\Gamma ds$ | 1.000000 |
| Interface flux range | $[-2.1982,\,-0.9832]$ |

Representative history rows (`output/history.csv`), showing the decay of the initialization/splitting transient:

| step | $t$ | $s_h$ | $s_{\text{err}}$ | $q_h$ | $q$ rel. error | $T$ max err | $T$ RMS | active dofs |
|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| 0 | 0.050 | 0.196081 | $2.8\times10^{-17}$ | $-2.19824$ | $-12.65\,\%$ | $1.13\times10^{-2}$ | $4.02\times10^{-3}$ | 910 |
| 1 | 0.052 | 0.200478 | $5.13\times10^{-4}$ | $-2.03229$ | $-6.39\,\%$ | $1.06\times10^{-2}$ | $4.57\times10^{-3}$ | 910 |
| 5 | 0.060 | 0.215914 | $1.12\times10^{-3}$ | $-1.82870$ | $-2.97\,\%$ | $8.28\times10^{-3}$ | $5.19\times10^{-3}$ | 975 |
| 19 | 0.088 | 0.261931 | $1.80\times10^{-3}$ | $-1.48863$ | $-1.52\,\%$ | $4.12\times10^{-3}$ | $3.15\times10^{-3}$ | 1170 |
| 39 | 0.128 | 0.315941 | $2.21\times10^{-3}$ | $-1.22274$ | $-0.50\,\%$ | $1.64\times10^{-3}$ | $1.20\times10^{-3}$ | 1430 |
| 59 | 0.168 | 0.361873 | $2.45\times10^{-3}$ | $-1.06719$ | $-0.42\,\%$ | $6.80\times10^{-4}$ | $4.06\times10^{-4}$ | 1625 |
| 74 | 0.198 | 0.392773 | $2.58\times10^{-3}$ | $-0.98318$ | $-0.40\,\%$ | $1.20\times10^{-3}$ | $4.29\times10^{-4}$ | 1755 |

### What each quantity establishes — and what it does not

* **Interface position error $0.657\,\%$** — accuracy of the *coupled* free-boundary solution against an exact solution of the same PDE. It is a single-number global measure; it mixes the thermal discretization error (through the flux) with the splitting error, and cannot separate them.
* **Temperature error (RMS $4.3\times10^{-4}$, max $1.13\times10^{-2}$)** — accuracy of the thermal field on the numerical domain, at the nodes currently inside the liquid. The maximum sits at the first step and near the interface (§3.5); the final-state RMS is the informative number, the maximum is dominated by the start-up transient.
* **Interface flux error ($1.35\,\%$ mean, $0.40\,\%$ final)** — accuracy of the *boundary data of the coupling*, i.e. of the recovered normal gradient. It is the most direct measure of the embedded-discretization quality and is the quantity that feeds the Stefan law; §3.8 shows it is the slowest to converge.
* **Linear-system residual $3.4\times10^{-16}$** — algebraic accuracy only. It states that the discrete linear system was solved, not that the discrete equations approximate the PDE.
* **Stefan residual $10^{-7}$** — consistency between the prescribed coupling law and the realized geometric motion of the level set (see §3.7). Not a physical closure test.
* **$y$-variation $3.7\times10^{-15}$** — the discrete solution preserves the one-dimensional symmetry of the exact problem exactly to roundoff, i.e. no spurious $y$-dependence is introduced by the cut geometry, the ghost penalty or the transport.
* **Mass identity $1.1\times10^{-15}$ and symmetry $\sim10^{-16}$** — the assembled discrete operators reproduce their defining integrals (Nitsche/ghost terms are symmetric and do not disturb the partition of unity), so the sparse matrices are faithful representations of the forms.
* **Cut-classification mismatch $0/4096$** — agreement of two *independent* geometric constructions on the same data (§3.10). It says nothing about the accuracy of the PDE solve.
* **Interface area $1.000000$** — the embedded boundary represents the interface without gaps or overlaps at this resolution (the exact value is the height, 1).

---

## 3.7 Interpretation of the Stefan residual

The code records

$$
r_S=\rho L V_n + k\,\partial_n T
$$

with $V_n$ reconstructed from successive *measured* interface positions, $V_n=\bigl(s_h(t^{n})-s_h(t^{n-1})\bigr)/\Delta t$, and $k\partial_nT$ replaced by the reconstructed flux of the previous step. Recorded values: mean $|r_S|=1.55\times10^{-7}$, maximum $|r_S|=7.68\times10^{-7}$ (vs $|q|\approx1$–$2$).

The reason it is so small is structural: the interface speed that is fed to the transport is *defined* by $V_n=-q_h/(\rho L)$, so the residual only measures how faithfully `eulerianTransport!` realizes the commanded velocity. That check is not vacuous — the measured motion agrees with the commanded motion to $1.5\times10^{-9}$, so the level-set update is essentially exact for this (linear, constant-velocity) configuration — but it is **not an independent measure of thermal accuracy**. A completely wrong conductivity or a mis-imposed interface temperature would still produce $r_S\approx0$, because the flux would be wrong consistently with the velocity it induces.

The signal that does carry independent information is the comparison with the analytical flux, $q_{\text{exact}}=k\,\partial_xT_{\text{exact}}(s_h(t),t)$: it measures the recovered gradient against the true gradient, independent of the coupling. That is why the flux error ($1.35\,\%$ mean, $12.65\,\%$ worst, $0.40\,\%$ final) is reported as a primary accuracy metric while the Stefan residual is reported as a coupling-consistency diagnostic.

Likewise, the residual can expose *numerical inconsistencies of the geometric update*: had the transport distorted the level set (it does not here — the interface stays linear to $5.6\times10^{-17}$), or had the velocity extension been wrong, the measured displacement would deviate from $V_n\Delta t$ and $r_S$ would grow. Its non-zero, slowly growing value ($10^{-10}$ early to $\sim10^{-6}$ late) reflects the increasing magnitude of $s_h$ and the finite-difference nature of the measured speed, not a physical defect.

---

## 3.8 Convergence and parameter sensitivity

### Time-step refinement at $65^2$ (recorded in `output/convergence.txt`)

| $\Delta t$ | $s(t_f)$ error | relative | mean $\lvert q\rvert$ rel. error | max $\lvert q\rvert$ rel. error | runtime |
|---:|---:|---:|---:|---:|---:|
| 0.002 | $2.577\times10^{-3}$ | $0.657\,\%$ | $1.349\,\%$ | $12.65\,\%$ | 33.6 s |
| 0.001 | $2.265\times10^{-3}$ | $0.577\,\%$ | $1.292\,\%$ | $9.71\,\%$ | 40.3 s |
| 0.0005 | $2.115\times10^{-3}$ | $0.539\,\%$ | $1.267\,\%$ | $7.82\,\%$ | 58.1 s |

The error decreases monotonically but with a strongly decreasing rate ($-0.080\%$, then $-0.038\%$), i.e. it is converging to a floor of about $0.5\,\%$ as $\Delta t\to0$ at this resolution. That floor must be **spatial** (cut geometry, weak Dirichlet, flux recovery) — or a coupling error tied to the interface flux, which also barely moves ($1.35\%\to1.27\%$ mean). Extrapolating the two-term behaviour would put the spatial limit near $0.50$–$0.52\,\%$; this is an extrapolation, not a measurement.

### Spatial refinement at $\Delta t=0.001$

| grid | $h$ | $s(t_f)$ error | relative | mean $\lvert q\rvert$ rel. error |
|---|---:|---:|---:|---:|
| $65^2$ | $1/64$ | $2.265\times10^{-3}$ | $0.577\,\%$ | $1.292\,\%$ |
| $129^2$ | $1/128$ | $1.276\times10^{-3}$ | $0.325\,\%$ | $0.667\,\%$ |

Both the interface error (ratio $1.78$) and the flux error (ratio $1.94$) improve under refinement, consistent with a first-order-like rate ($p\approx0.8$–$0.9$ by these two points). **This does not establish an asymptotic order**: two resolutions cannot distinguish $O(h)$ from $O(h^{1/2})$ with a favourable constant, cut-cell geometry changes discontinuously as the interface crosses nodes (the smallest cut-cell fraction over the run is $5.2\times10^{-5}$), and no self-similar refinement sequence was run. A proper order study needs at least three grids per family and control of the interface's position relative to the nodes.

### Nitsche penalty sensitivity at $65^2$, $\Delta t=0.002$

| $\gamma_d$ | $s(t_f)$ relative error | mean $\lvert q\rvert$ rel. error | max $\lvert q\rvert$ rel. error |
|---:|---:|---:|---:|
| 10 | $0.928\,\%$ | $2.046\,\%$ | $22.01\,\%$ |
| 30 | $0.687\,\%$ | $1.423\,\%$ | $12.62\,\%$ |
| 100 | $0.657\,\%$ | $1.349\,\%$ | $12.65\,\%$ |
| 300 | $0.650\,\%$ | $1.329\,\%$ | $12.66\,\%$ |

The dependence is real and one-sided: weak enforcement is only *accurate* when the penalty is stiff enough that the discrete interface value is close to $T_m$. With $\gamma_d=10$ the mean interface value is measurably above $T_m$ (an earlier diagnostic on this configuration, computed during development, gave $\langle T_h\rangle_\Gamma=1.4\times10^{-2}$ at $\gamma_d=10$ versus $T_h$ scale 1), the discretely imposed data are then effectively offset, and the recovered gradient inherits the error. Above $\gamma_d\approx100$ the result saturates ($0.657\%\to0.650\%$), so further stiffening buys nothing at this resolution.

**But the penalty alone is not the explanation of the residual error.** Reasons, all supported by the recorded numbers: (i) even at $\gamma_d=300$ the interface error stays at $0.65\%$ and the flux error at $1.33\%$; (ii) the flux error is largest ($12.7\%$) at the first step, where the geometry is frozen and the *splitting* lag dominates; (iii) the flux error is $0.40\%$ at $t_f$ while the interface error is $0.657\%$, so the accumulated interface error is not simply the instantaneous flux error integrated. The plausible irreducible contributors are the Q1 background space (a $C^0$ interpolant cannot represent the kink of the true solution at the interface, so the recovered normal gradient on the liquid side of a cut cell is intrinsically $O(h)$-accurate at best) and the interface representation ($\phi$ linear per cell, interface located by linear interpolation).

---

## 3.9 Stability, geometry, and conditioning

Recorded indicators over the 75-step run:

* **Active dofs** grow from 910 (14 node columns × 65) to 1755 as the liquid grows; total dofs 4225. The active set is detected from the diagonal of $M+\Delta tK$ (dofs whose cells hold no liquid have an empty equation and are set to $T_m$ by an identity row).
* **Smallest cut-cell volume fraction** $5.168\times10^{-5}$ of a background cell (i.e. $0.005\,\%$). Such slivers appear whenever the interface passes within $5\times10^{-5}h$ of a node; 64 cut cells are present at every step.
* **Condition estimate** of the active system: $\kappa_2=3.511\times10^{2}$ at $t=0.050$ (910 dofs), computed by dense SVD of the reduced active block. For the last step (1755 dofs) and for both $129^2$ runs the estimate is deliberately reported as `NaN`: `conditioningEstimate` refuses to form a dense $>1200\times1200$ block. So the conditioning statement is *one measurement on the smallest active set*, not a history; the minimum diagonal ratio over the run ($9.7\times10^{-4}$) is the only per-step conditioning proxy recorded.
* **Algebraic robustness at slivers**: the largest relative residual of any step is $3.4\times10^{-16}$, including steps with cut fractions near $10^{-4}$; the implicit (mass-carrying) operator plus the ghost penalty prevented any observed ill-conditioning failure or iteration breakdown, and no step needed a special treatment.
* **Operator structure**: mass identity to $1.1\times10^{-15}$ (partition of unity on the cut subcells), symmetry to roundoff ($M$ exactly symmetric, $K$ to $1.07\times10^{-16}$ — the residual asymmetry is assembly-order roundoff, checked by an assertion at every step).
* **Geometric cross-check**: 0 mismatches of the per-cell IN/CUT/OUT state between the native EvolvingDomains cutter and GridapEmbedded, on all 4096 background cells, at every step.
* **Level set**: remains exactly linear ($\phi$ zero crossing at mid-height equals the liquid area to $5.6\times10^{-17}$), and the transport realizes the commanded displacement to $1.5\times10^{-9}$. Richardson-type "CFL" of the transport, derived from recorded values: $\max(V_n)\Delta t/h=2.198\cdot0.002/0.015625=0.28$ at $65^2$ ($0.56$ at $129^2$), below the order-one threshold quoted in the package README.

**What this does and does not establish.** These numbers show that *this* configuration (fixed uniform mesh, Q1, $\gamma_d=100$, $\gamma_g=0.1$, a straight interface) is stable and algebraically well behaved, including for cut fractions down to $5\times10^{-5}$. They do **not** establish uniform stability in the cut-cell measure, for other polynomial degrees, or for strongly varying cut fractions; the classical cut-FEM theory guarantees the method is stable independently of the cut position only when the ghost term is present *and* the penalty is scaled with the local cut configuration, and no such scaling study was run. Also, no test probed a cut fraction below $5\times10^{-5}$ deliberately (the value is incidental), and no eigenvalue or condition-number sweep as a function of the interface's distance to nodes exists yet (Recommendation F).

**Role of the two stabilizers.** The Nitsche terms supply the Dirichlet data and the corresponding consistency flux; without a sufficiently large $\gamma_d$ the discrete boundary value drifts (§3.8). The ghost penalty is not a boundary condition but a mesh-dependent bilinear form that couples neighbouring cells across cut facets, restoring control of the dofs that have (almost) no support in the liquid; it is what allows slivers not to pollute the solution, and it is the reason the algebraic residual stays at machine level.

---

## 3.10 Cross-package geometric interoperability

### The same array consumed twice

At every step the level-set samples $\phi_h$ (a `Matrix{Float64}` inside `state.level.data[:init]`) are consumed by two constructions that share no code:

$$
\phi_h \;\longrightarrow\;
\begin{cases}
\text{EvolvingDomains native `Cut` (C++ CutCells backend via C ABI)}, & \text{used as the reference and for comparison},\\[2pt]
\text{GridapEmbedded `DiscreteGeometry`} \to \texttt{cut(bgmodel, geo)}, & \text{used by the PDE solver}.
\end{cases}
$$

The comparison is made per background cell, matched by cell centre so that the two cell orderings are irrelevant (`compareCuttings`): the state codes coincide convention-wise (`IN`$=-1$, `CUT`$=0$, `OUT`$=1$ in both packages), and the recorded mismatch count is **0 of 4096 cells at every one of the 75 steps**.

### Why this is stronger than passing a solver object around

Exchanging a solver object between kindred packages proves only that the object survived. Here the exchanged quantity is the *mathematical content*: nodal values of a scalar field plus node coordinates. Neither package knows the other exists; there is no common abstract interface, no shared geometry type, no plugin mechanism, no adapter registered in either package. Agreement on the classification of all 4096 cells therefore demonstrates agreement on the meaning of the data — sign convention ($\phi<0$ = material), the nodal basis (piecewise-linear interpolation of the nodal values), the cell subdivision rule, and the background mesh indexing — rather than agreement of implementation details.

That is the architectural point: interoperability of this kind requires the components to agree on a *mathematical representation*, not on an API. EvolvingDomains' public surface admits this because a `Field` is essentially an array plus a grid description, so its content can be handed to any consumer that accepts arrays.

### What the comparison does *not* establish

Cell classification is a coarse invariant. The experiment does **not** show that the two packages produce the same interface quadrature, the same facet parametrisation, the same normals, or the same treatment of cells cut by more than two planes, nor does it compare areas, centroids or normals of the interface facets. Two independent consistency signals were recorded, and they are weaker than a full geometric equivalence:

* the interface area from GridapEmbedded integrates to $1.000000$ (exact value for a straight interface of height 1), so its facet set has neither gaps nor double counting at this resolution;
* the liquid area equals the level-set mid-height zero crossing to $5.6\times10^{-17}$, so the cut volume and the interpolated contour agree for this geometry.

For a *curved* interface, volume, facet area, normal orientation and the polygonisation error would all differ between implementations at the $O(h^2)$ level, and cell classification would no longer be an adequate comparison (Recommendation C). In three dimensions the package itself documents unwelded interface vertices, which would make quadrature comparisons essential.

### Cost of the adapter

The adapter is ~20 lines of index bookkeeping (`nodeIndexMap`, `dofNodeMap`) plus the three-line value hand-off, and it is written by the application, not by either package. Both index maps are validated against coordinates at startup and *fail loudly* if the conventions differ. In this configuration the two deviations are exactly zero, so the maps degenerate into identity permutations; that is fortunate, not guaranteed — a different node ordering (e.g. first-axis-fastest vs last-axis-fastest) or a different axis order would require a permutation, which the same code computes from coordinates without knowing either convention.

The honest assessment is therefore two-sided. EvolvingDomains' minimalism means there is no opaque object standing between the application and the consumers of its data — the coupling is possible at all only because the state is exposed as arrays. But minimalism also means the *responsibility* for coordinate conventions, dof ordering and state transfer is entirely the application's: the adapter does not disappear, it moves into user code (here, into `run.PlanarStefan.jl`). Whether that is a gain depends on whether the application needs its own time loop and its own PDE — which is exactly the case this case study targets.

---

## 3.11 Architectural assessment

### Minimalism

The public data model is three types ($\texttt{Grid}$, $\texttt{Discretization}$, $\texttt{Field}$) that organize axes, a level dictionary and dense arrays, plus a `Cut` snapshot and the transport operators. In this experiment that was **sufficient and never obstructive**:

* the level set and the temperature live as plain arrays with a documented indexing convention, so they could be handed to Gridap/GridapEmbedded and to the linear solvers without any conversion layer;
* the `Cut` snapshot was usable as an independent geometric reference without a solver attached to it;
* the transport operator could be called with an application-built velocity field, one call per step, with the timestep chosen by the application (the package deliberately does not choose one).

Nothing in the workflow needed a "simulation", "problem" or "solver" type from EvolvingDomains. The absence of such a type is what allowed the thermal PDE to be discretized by an unrelated package: no plug-in point had to be respected, and no internal data structure had to be reconstructed.

The price of minimalism showed up in the same places that are usually attributed to it: the package does not offer the level-set reinitialization, velocity extension, or mesh transfer needed for a general moving-domain PDE, so the application has to know when each is needed and what it means. In this planar case the omissions were harmless (constant extension is exact; reinitialization is unnecessary for an exactly linear $\phi$), but that is a property of the test, not of the architecture.

### Composability

The demonstration supports the composability claim in a concrete, checkable way:

* the thermal PDE was discretized **entirely outside** EvolvingDomains (Gridap/GridapEmbedded, no EvolvingDomains type in the assembly);
* the geometry and kinematics were handled **entirely outside** Gridap (EvolvingDomains, no Gridap type in the transport);
* the coupling was mediated by an array of values plus coordinates (`φnode`, `X`), verified to be an exact match;
* the assembled operators remained ordinary sparse matrices, written to `output/operators.csv` and re-usable; the thermal solve used Julia's own LU factorisation, not a package-hidden solver;
* the time loop, the flux functional, the Stefan law and the geometry update are application code (~40 lines of physics plus diagnostics);
* **no package source file was modified** to make this work, and no private API was used (see the boundary note below).

The most informative negative detail is what the adapter *did* need: coordinate matching and dof-ordering verification (§3.10). Those are unavoidable whenever two independent discretizations meet, and they are cheap here only because both sides use uniform Cartesian meshes with a nodal Q1 space.

### Practical limitations observed

1. **Adapter is application code.** ~20 lines of index bookkeeping plus verification; it would grow (permutations, interpolation, dof layouts) for non-matching orderings or higher-order spaces.
2. **Native backend deployment.** The `cut` snapshot needs the C++ CutCells shared library. This checkout has no published `Artifacts.toml`, so the library was built locally (`scripts/build_cxx_bindings.sh`, install prefix outside the source tree) and is loaded via `EVOLVINGDOMAINS_CUTCELLS_LIBRARY` (the run script falls back to `./build/`). This is a deployment requirement, not a numerical one, but it must not be overlooked when reproducing the case.
3. **Uniform Cartesian grids only.** Both cutters rely on uniform spacing; the adapter's power-of-two coordinate matching would need to be replaced for graded meshes.
4. **The moving interface is planar and axis-aligned**, so the geometry never exercises curvature, the normals are constant, the velocity extension is exact, and the level set never needs reinitialization.
5. **The thermal solver is untested on curved or topology-changing domains**; nothing here validates cut-cell aggregation, facet welding or normal orientation for a genuinely curved front.
6. **The flux recovery is sensitive** to spatial resolution ($1.27\%$ at $65^2$, $0.67\%$ at $129^2$) and to the boundary treatment ($2.0\%$ at $\gamma_d=10$, $1.33\%$ at $\gamma_d\ge100$), and it is the error source that propagates into the interface motion.
7. **Implicit-interface caveat**: the discretized operators depend on the geometry, and the geometry is updated explicitly; the coupled system is not solved monolithically, so discrete conservation of energy across the interface is not enforced to machine precision (it was not measured).

Assessment: within its stated scope — *arrays and geometry as the contract* — EvolvingDomains composed with an unrelated embedded-FEM stack without friction and without modification. Outside that scope (curvature, reinitialization, topology change, ALE motion, non-Cartesian meshes) nothing was demonstrated, and the applications-level work (adapters, conventions, stabilization choices) remains exactly where minimalism leaves it: with the user, who must also be the one to verify it.

---

## 3.12 Existing visualizations

Generated by `visualize.PlanarStefan.jl` from the saved CSVs only (no PDE solve, no field update). Paths verified relative to this file:

**Temperature heatmaps with the moving interface** — `output/fig_heatmaps.png`

![Temperature heatmaps with numerical and analytical interface](output/fig_heatmaps.png)

Four snapshots (labels $t=0.05,0.09(2),0.15(2),0.20$). Colour: temperature; gray: solid phase ($\phi>0$, blanked); white contour: numerical interface $\phi=0$; black dashed: analytical interface $x=s_{\text{exact}}(t)$. The temperature decays monotonically from $T_h=1$ at the wall to $T_m=0$ at the interface, the boundary layer thickens in time, and the two interfaces coincide within a pixel. **Caveat from §3.5:** in these panels the temperature is at the record's label while the level set (and hence the white contour) is one step behind; the measured offset is $-0.25$ cells at $t=0.052$, passing through zero near $t\approx0.14$ and reaching $+0.04$ cells at $t=0.200$. The final panel therefore shows the interface at $0.392773$, not the reported post-update value $0.394739$.

**Interface position, error, flux, Stefan residual** — `output/fig_history.png`

![Interface history](output/fig_history.png)

(a) $s_h(t)$ against the similarity solution — the two curves are visually indistinguishable, the numerical one slightly above after $t\approx0.1$; (b) $s_h-s_{\text{exact}}$, growing monotonically to $2.58\times10^{-3}$; (c) interface flux, numerical against analytical, showing the $12.7\%$ first-step overshoot and its decay to $0.4\%$; (d) $|r_S|$ on a log scale, $10^{-10}$–$10^{-6}$, with the sawtooth caused by differencing successive measured interface positions (value grows with $s_h$, not with a physical defect; see §3.7).

**Temperature profiles at mid-height** — `output/fig_profiles.png`

![Temperature profiles](output/fig_profiles.png)

Markers: numerical nodal values; lines: analytical profile (drawn only where $\phi<0$, since the analytical continuation is unphysical inside the solid); dash-dotted: analytical and numerical interface positions. The numerical markers lie on the analytical curves inside the liquid to within the errors of §3.6. Two artefacts are visible and are genuine, not plotting errors: (i) the first node on the solid side of the interface carries a small negative temperature (about $-3\times10^{-3}$) because the background Q1 function has a dof outside the liquid whose value is determined by the cut-cell equations — the numerical "extension" is much flatter than the analytical continuation ($-2.7\times10^{-2}$), and this is precisely the mechanism limiting the recovered flux; (ii) the plateau at $T_m=0$ in the solid is the imposed solid state. The staggering of §3.5 also applies to the profile records: the drawn temperature field is one step ahead of the saved level set.

---

## 3.13 Recommendations for the second campaign

Documented, not executed.

**A. Synchronize saved fields.** Record snapshots and profiles *after* `advanceGeometry!` (or record both states with their own time labels), then re-plot and quantify the change. *Uncertainty addressed:* reporting integrity — how much of the visible interface offset in the figures is lag rather than error (measured upper bound here: $0.25$ cells early, $0.04$ cells late).

**B. Improved interface heat-flux reconstruction.** Compare the present background-gradient functional $k\langle\nabla T_h\cdot\mathbf n\rangle_\Gamma$ with consistent alternatives: the Nitsche multiplier (the boundary functional that the weak form itself defines), a one-sided polynomial recovery on the liquid side of the interface, and (if feasible) an $h$-refined or $p$-enriched reconstruction near $\Gamma$. *Uncertainty addressed:* whether the $1.3\,\%$–$0.7\,\%$ flux error is a property of the cut-FEM space or of the recovery operator — i.e. how much of the interface error is inherited boundary-data error.

**C. Curved melting interface.** Replace the straight front by a circular (or otherwise curved) melting front, e.g. by a radially symmetric similarity problem or by imposing a curved initial $\phi$. Compare EvolvingDomains and GridapEmbedded beyond cell classification: interface length, enclosed area, facet centroids, normals and orientation, plus the resulting thermal solution. *Uncertainty addressed:* whether the demonstrated interoperability extends to nontrivial geometry, where normals vary and the two polygonisations must agree quantitatively; also whether the flux recovery degrades faster with curvature.

**D. Higher-order or refined spatial discretization.** At least three grids ($65^2,129^2,257^2$) at fixed $\Delta t$ to estimate an order, and if feasible $Q_2$ on the background mesh with a dof transfer based on pointwise evaluation instead of node identity. *Uncertainty addressed:* whether the observed $0.5\,\%$ floor is an $O(h)$ cut-FEM effect that will disappear or an $O(h^{1/2})$ effect that will not; and whether the reduced $y$-symmetry and mass identities survive non-nodal dof layouts.

**E. Temporal coupling schemes.** Compare the present first-order partitioned update with (i) a predictor–corrector that moves the interface with the previous speed before the thermal solve and corrects afterwards, and (ii) a Gauss–Seidel/monolithic iteration between the thermal solve and the interface position within each step. *Uncertainty addressed:* how much of the $0.657\,\%$ interface error and of the $12.7\,\%$ first-step flux error is splitting rather than discretization; whether the observed $0.5\,\%$ floor moves.

**F. Cut-cell stabilization.** Sweep $\gamma_d$ and $\gamma_g$ jointly, and place the interface deliberately at $10^{-1},10^{-2},\dots,10^{-6}$ of a cell away from nodes, recording condition estimates, residual and flux error. *Uncertainty addressed:* whether the ghost-penalty/nitsche combination is uniformly stable as the cut fraction tends to zero, and what the smallest reliable cut fraction is — the present run only *happens* to reach $5.2\times10^{-5}$ and reports a single $\kappa$ estimate at 910 dofs.

Cross-cutting: record $\kappa$ (or a Krylov-based estimate) at every step or at least at a representative sample, so that conditioning statements are not single-measurement anecdotes.

---

## 3.14 Conclusions

**Demonstrated.**

* EvolvingDomains supplied the moving-domain representation — grid, level-set field, sampled temperature field, native cut snapshot — and advanced the interface with its own transport operator, at a cost of one call per step, with the time step chosen by the application. Its public data model was sufficient for the whole workflow; no internal structure had to be reconstructed and no source file of any dependency was modified.
* An independently developed embedded finite-element stack (Gridap + GridapEmbedded) constructed the cut liquid domain, assembled the mass and conductivity operators with Nitsche boundary conditions and ghost stabilization, and solved the thermal problem implicitly, using only the exposed level-set samples and node coordinates as input.
* The assembled mathematical operators remained explicit sparse matrices, re-usable outside the packages (`output/operators.csv`), and satisfied their defining invariants (partition of unity to $1.1\times10^{-15}$, symmetry to roundoff).
* The coupled solution tracks an analytically known melting front: $+0.657\,\%$ interface error at $t_f$, temperature RMS $4.3\times10^{-4}$, $y$-variation at roundoff, flux error $1.35\,\%$ mean ($0.40\,\%$ final).
* Two independent cut-cell constructions agreed on the phase state of all 4096 background cells at every step, on data exchanged as arrays rather than as objects.
* The time loop is application code and reads as the mathematical sequence (geometry, thermal solve, flux, Stefan speed, velocity extension, level-set update, time advance), which was an explicit design objective.

**Not demonstrated.**

* Synchronized reporting: the saved temperature and level-set fields are one step out of phase, and the interface drawn in the figures is not the post-update value used for the headline error.
* High-order or geometry-independent convergence: two grids and two time steps are not enough for an asymptotic order, and the cut-cell geometry changes discontinuously with resolution.
* Robustness under curved or topology-changing interfaces; the manifest of the experiment is a straight, axis-aligned front, for which velocity extension and level-set regularity are trivially exact.
* Independence from native backend deployment: the CutCells library had to be built locally in this checkout.
* Universally reliable cut-FEM conditioning or interface flux recovery: one condition estimate ($3.5\times10^{2}$ at 910 dofs), no sweep in the cut fraction, and a flux functional whose error depends on the Nitsche penalty and on $h$ with no established order.

**Why the interoperability result is nevertheless meaningful.** The failure modes that make cross-package coupling expensive in practice — opaque solver objects, private geometry caches, hidden conventions, mandatory plug-in points — are all absent here. The two packages met at the lowest possible level of abstraction (values, coordinates, cell classification) and, for this problem class, that level was sufficient to run a genuine coupled free-boundary simulation with an explicit, auditable time loop and accessible operators. The experiment therefore supports the architectural hypothesis that a minimal shared data model plus ordinary numerical arrays is a workable substrate for composing independently developed PDE software — while making equally clear that generality (curvature, $p$-adaptivity, non-Cartesian meshes, topology change) would move work, not remove it, into the application's adapters and diagnostics.

---

### Appendix: reproduction

```sh
# environment (self-contained; does not touch the package environment)
julia --project=validation/PlanarStefan -e 'using Pkg; Pkg.instantiate()'

# native cut-cell backend, if not already present in ./build
BUILD_DIR=/tmp/edbuild PREFIX=/tmp/edinstall bash scripts/build_cxx_bindings.sh
cp /tmp/edinstall/lib/libevolvingdomains_cut_cells.* validation/PlanarStefan/build/

# default run (65^2, dt = 0.002) + figures
julia --project=validation/PlanarStefan validation/PlanarStefan/run.PlanarStefan.jl

# figures only, from saved data
julia --project=validation/PlanarStefan validation/PlanarStefan/visualize.PlanarStefan.jl \
      validation/PlanarStefan/output

# refinement / sensitivity runs recorded in output/convergence.txt
julia --project=validation/PlanarStefan validation/PlanarStefan/run.PlanarStefan.jl \
      --no-viz --dt 0.001 --out validation/PlanarStefan/output/dt0.001
julia --project=validation/PlanarStefan validation/PlanarStefan/run.PlanarStefan.jl \
      --no-viz --n 129 --out validation/PlanarStefan/output/n129dt0.002
julia --project=validation/PlanarStefan validation/PlanarStefan/run.PlanarStefan.jl \
      --no-viz --gd 10 --out validation/PlanarStefan/output/gd10
```

Output files per run: `history.csv` (per-step diagnostics), `profiles.csv` (mid-height profiles), `snapshots.csv` (full nodal fields, see the staggering note in §3.5), `operators.csv` (final-step $M$ and $K$ triplets), `summary.txt`, and — for the default run — the three figures.
