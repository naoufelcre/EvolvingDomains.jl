# Zalesak slotted disk / sphere

**Verdict: accepted** for every variant except the 96³ tilted-axis stress run, which is an
exploratory stress finding. One case, several variants: dimension, solver
and resolution change, the object under test does not — advection of a non-smooth interface
through a closed incompressible velocity field, without reinitialization.

## Structure

Three small files hold what is shared, and each variant is a plain script that runs the
simulation itself:

| File | Contents |
|---|---|
| `zalesak.jl` | the analytic problem only: `disk`/`sphere`/`tilted_sphere`, their rotated exact solutions, the velocity fields, the exact area/volume |
| `measure.jl` | `phase_measure`, `field_norms`, `interface_points`, `interface_distance`, and `within` for one gate |
| `artifacts.jl` | `figure_2d`, `figure_3d`, `figure_comparison_2d/3d`, `save_figure`, `report` |
| `runZalesak2D.jl` … `runZalesak3DTenRotations.jl` | one file per variant: initial data → velocity → discretization → time loop → measurements → gates → figure |
| `TrixiZalesakSolver.jl`, `TrixiZalesak3D.jl` | the Trixi DGSEM implementation and its exchange self-check |

The scripts share no driver, no configuration object and no scheme abstraction: each one
shows the five steps in order, so a reader sees exactly which call does what. Two variants of
the same dimension therefore repeat the time loop — deliberately, because that loop is the
thing the case demonstrates. All variants sit side by side in this directory: the file name
says which one it is.

There is **one refinement level per variant**. Acceptance gates are stated numbers attached
to that level, checked with `within`. There is no refinement study: monotone decrease over a
few resolutions is weak evidence and is not what these cases claim.

Each script activates the pinned environment itself (`validation/env`), so everything runs
from the repository root with no flags:

```sh
julia validation/zalesak/runZalesak2D.jl [NEW_output]
```

## Variants

| Script | What it adds | Verdict |
|---|---|---|
| `runZalesak2D.jl` | Field WENO5, 100 cells/axis, ten rotations | accepted |
| `runZalesak2DTrixi.jl` | independent DGSEM solver, same ten rotations | accepted, independent comparison |
| `runZalesak3D.jl` | canonical 3D sphere, 50³, one revolution | accepted |
| `runZalesak3DTrixi.jl` | independent DGSEM solver, same one revolution | accepted, independent comparison |
| `runZalesak3DTenRotations.jl` | tilted-axis stress, 96³, ten rotations | accepted at 96³; exploratory |

## 2D, Field — `runZalesak2D.jl`

Recorded ten-rotation result (100 cells): sampled measure drift **+0.9202%**, symmetric-
difference measure **0.00453125**, interior L1 **0.00096009**, edge-sampled interface distance
**0.02474645** (not exact Hausdorff).

## 2D, Trixi — `runZalesak2DTrixi.jl`

Field WENO5 and Trixi DGSEM start from the same analytic slotted-disk function and are
sampled at the same ten rotation times. Recorded result: measure drift **+0.9202% /
−1.4648%**, symmetric-difference measure **0.00453125 / 0.0017**, interior L1
**0.00096009 / 0.00033357**, L1 difference between the two solvers **0.000997**
(Field / Trixi). Resolutions are 10,201 versus 16,384 scalar DOFs.

**Limits:** different boundary discretizations, different resolutions, sampled phase measure,
no exact Hausdorff bound, no matched-cost claim. The comparison runs in the pinned
`validation/env` project (Julia 1.11.7, Trixi 0.17.21, CairoMakie 0.15.17); both runners
record the versions and source hashes in `run.txt`. Field clamps off-grid stencil reads while
Trixi DG uses the interior boundary trace, so whole-box error norms reflect boundary
differences, not only the schemes.

## 3D canonical, Field — `runZalesak3D.jl`

Domain `[-1,1]^3`, sphere center `(-1/3,0,0)`, radius `0.5`, slot `|x+1/3|<1/6`, `y∈[0,1]`
through all `z`, velocity `π(y,-x,0)`, one revolution `T=2`.

Recorded result: sampled measure drift **−1.53%**, symmetric-difference measure **0.018640**
(analytic `0.397547`), interior L1/L∞ **0.006384/0.060490**. The slot remains visible in the
surface and the slice.

**Limits:** one revolution only; six-tet piecewise-linear cut; unwelded topology; sampled
measure; boundary differs from a true Neumann implementation. The 3D figures need the native
`CutCells` library: the scripts use `EVOLVINGDOMAINS_CUTCELLS_LIBRARY` when set, otherwise the
local build from `scripts/build_cxx_bindings.sh`, and skip the figure with a message if neither
exists. The numbers never depend on it.

## 3D, Trixi — `runZalesak3DTrixi.jl`

Trixi DGSEM (16³ elements, degree 3, scalar-only Lax–Friedrichs, `SimpleSSPRK33`,
interior-trace boundary) sampled onto the common 50³ Field grid. Recorded result: measure
drift **+0.443%**, symmetric-difference measure **0.009328**, interior L1/L∞
**0.0029566/0.025252**. 262,144 scalar DOFs versus Field's 132,651.

**Limits:** analytic Trixi initialization, not a 3D Field→DG exchange; one revolution; no
convergence claim; sampled phase diagnostics.

## 3D tilted-axis stress — `runZalesak3DTenRotations.jl`

Sphere centered at `(0.5,0.75,0.5)`, radius `0.15`, through-slot 0.05 wide, rotated about
`(1,1,1)/√3` for ten rotations at 96³, where the slot is about **4.8 cells wide**.

| Resolution | Steps | Solver time | Volume drift | Symmetric difference | L1 | L∞ |
|---|---:|---:|---:|---:|---:|---:|
| 96³ | 17,440 | 792 s (8 threads) | **−3.36%** | 0.001358 | 0.000345 | 0.02459 |

An under-resolved 48³ arm was tried earlier and lost the slot (drift −37%); it is not part of
the case, because an under-resolved run answers no question the case asks. The 96³ result is a
stress finding, not a convergence proof.

## Adding a variant

Copy the nearest `runZalesak*.jl`, change the initial data, velocity and resolution, keep the
five steps and the `within` gates, add a row here. Shared mathematics goes to `zalesak.jl`, shared
measurement to `measure.jl`.
