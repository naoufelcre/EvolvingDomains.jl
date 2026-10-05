# ![visuel](concept.svg)  EvolvingDomains.jl

**A Julia package for solving PDEs on moving domains.**

This package provide a set of utilities to write 2D multiphysics moving domain problems.

The paradigm the package is built on is a decoupling between kinematics and dynamics. 

The package provides tools to handle the kinematics side by dedicated geometric structures. The dynamics part is intended to be handled on the active mesh, the package remain at the data-structure level, so you can proceed any way you like.  

In particular, the package was designed to work within the [`Gridap`](https://github.com/gridap/Gridap.jl) ecosystem and specifically their embedded finite element extension [`GridapEmbedded`](https://github.com/gridap/GridapEmbedded.jl).

`EvolvingDomains` is a registered package ! You can install it via the Julia REPL

```julia
# Type ] to enter package mode
pkg> add EvolvingDomains 
```

For a full working example see `TestDumbellParabolic.jl` that recreates a test case from the 2025 paper of Olshanskii & Reusken. [arXiv:2504.14116](https://arxiv.org/pdf/2504.14116)

![Temperature evolution - Olshanskii & Reusken test case](TEMPERATURE_OLSHANSKII_REUSKEN.gif)

## Geometric tooling

The basic object of the package is the `EvolvingDiscreteGeometry`. It is an all-in-one object for basic routines regarding implicitly defined level-set geometry that evolves.

In particular it provides the following functionalities:

- **Reinitialization of the level set** to a signed distance function.

  After advection the level-set gradient `|∇φ|` drifts away from 1. Reinitialization restores the signed distance property by solving the Eikonal equation `|∇φ| = 1` on the background grid.
  The implementation uses the **Fast Sweeping Method** (Zhao 2005).

  ```julia
  reinitialize!(geom)   # restores |∇φ| ≈ 1 everywhere
  ```

- **A lazy cache system** for stencils and derived data.

  `EvolvingDiscreteGeometry` holds a `GeometryCache` that stores the GridapEmbedded cut geometry, the active node set, and the transfer and extension operators. Cache-aware mutating operations such as `set_levelset!` and `reinitialize!` invalidate derived entries, which are recomputed on first access. Direct mutation of `geom.levelset` is deliberately low-level and must be followed by `invalidate!(geom.cache)`.

  ```julia
  set_levelset!(geom, phi_new)        # updates φ, invalidates all cache entries
  indices = get_active_indices(geom)  # recomputes and caches active IN+CUT nodes
  ```
  **Active indices** are used to couple effectively with the CutFEM method provided by `GridapEmbedded`.

  Many moving domain problems involve fields living on the geometry. To handle this the package provides a dedicated structure `CartesianMeshField`. It wraps the flat nodal data array and provides clamped 2D indexing and a bilinear interpolant (via `get_interpolator`), which is used internally by both the WENO5 stencils and the transfer operators.

- **A robust explicit curvature handling** 

  Because curvature is an essential modeling asset, we provide a simple way to compute it from the evolving discrete geometry. Our goal is to provide a simple method for fast prototpying, However to fit the low level philosophy, it's not plug and play for a semi implicit approach.

- **A Topological filter**
  To remove subgrid artifcats we have a dedicated topological filter, together with reinitialization of the SDF property, this module is to restore good health of the level set function.

## Kinematic tooling

## Transfer

Transfer between different levels of discretization. The package provides a `GridMeshTransfer` operator that follows the [`TransferOperator.jl`](https://github.com/naoufelcre/TransferOperator.jl) protocol, exposing two directions:

- **`restrict` (Grid → Mesh):** evaluates the `CartesianMeshField` at FE mesh nodes via bilinear interpolation, then projects into the target `FESpace` using Gridap's `interpolate`.
- **`prolong` (Mesh → Grid):** maps DOF values back onto the background grid using direct index mapping when the mesh topology allows it, falling back to batch point evaluation otherwise.

```julia
transfer_op = setup_transfer(geom, V)      # V is the current AgFEM FESpace
u_mesh = grid_to_mesh(geom, u_grid)        # restrict: Cartesian field → FE function
u_grid = mesh_to_grid(geom, u_mesh)        # prolong:  FE function  → Cartesian field
```

### Level-set advection — WENO5 + SSP-RK3

We advance the level set by solving the *transport* equation on the domain
```math
∂φ/∂t + v·∇φ = 0`. 
```

The spatial discretization uses the **fifth-order WENO** scheme (Jiang & Shu 1996) with Jiang-Peng smoothness indicators (Jiang & Peng 2000): at each node the upwind-biased directional derivative is selected based on the sign of `v`, and non-linear weights suppress oscillations near discontinuities while recovering fifth-order accuracy on smooth regions. 

Time integration uses the **third-order Strong-Stability-Preserving Runge-Kutta** (SSP-RK3, Shu & Osher 1988).

It can be used with any well-defined velocity field, sampled onto the grid via `sample_velocity`:

```julia
using EvolvingDomains.Geometric: CartesianMeshField

vel = StaticFunctionVelocity(x -> VectorValue(-ω*(x[2]-0.5), ω*(x[1]-0.5)))
v_nodes = sample_velocity(vel, grid_info(grid), t)
v_field = CartesianMeshField(v_nodes, grid_info(grid))
advance!(geom, v_field, Δt)   # WENO5 + SSP-RK3 step on the level set
```

We refer to the classical text *Osher, S., & Fedkiw, R. (2002).* `Level Set Methods and Dynamic Implicit Surfaces`

### Fields advection: CCISL or CIP

For a scalar field `c` we expose two methods to solve either the *continuity* equation
```math
\partial_t c + \nabla\cdot(cv) = 0.
```
or the *transport* equation
```math
\partial_t c + v\cdot\nabla c = 0,
```

Differently to the previous method for the transport of the level set on the whole domain, those methods are aware of the implicitly defined geometry and they are intended to be coupled to the deforming geometry.

To advect a field
```julia
ensure_cut!(geom)                               # preserve Ωⁿ on the next update
advance!(geom, v_field, Δt)                     # construct Ωⁿ⁺¹
k_map = TransportMap(geom, v_field, Δt)         # validates velocity grid metadata
advect!(new_field, current_field, k_map; type=:conservative) #For the continuity equation
```
or
```julia
advect!(new_field, current_field, k_map; type=:intensive) #For the transport equation
```

This implementation uses a first-order x-then-y Lie split with Euler characteristics. The velocity is a frozen nodal vector field for the step, normally the same sampled field passed to `advance!`.

Bare velocity vectors are interpreted in the source field's grid ordering; wrapping the velocity in a `CartesianMeshField` additionally validates its grid metadata.

```julia
next = CartesianMeshField(similar(current.data), current.grid)
advect!(next, current, v_field, Δt)  # type=:intensive is the default
```

#### CCISL

We implements the **Conservative Cell-Integrated Semi-Lagrangian** method (Lentine, Grétarsson & Fedkiw 2011).

#### CIP

`WORK IN PROGRESS`

We implement a *Directionally split Constrained Interpolation Profile* method. 

Without a cache, the Hermite interpolation profile is reconstructed from `current` on every call. An optional cache carries the profile and reuses all scratch arrays:

```julia
cache = CIPCache(current)
for step in 1:nsteps
    advect!(next, current, v_field, Δt; cache=cache)
    current, next = next, current
end
```

The cache keeps a snapshot of its last output. If another operator changes the next source field, `advect!` detects the changed nodal values and rebuilds the profile.
Cached and uncached repeated transport are different discretizations: the cached form carries the CIP moments, while the uncached form reconstructs them each step.



## Examples and independent packages

The numerical examples run without a terminal renderer:

```sh
julia --project=. test/TestGeometryEvolution.jl
julia --project=. test/TestDumbellParabolic.jl
julia --project=. test/TestHeleShawST.jl
```

Terminal plotting is not part of the ED API. ED does not depend on TPlot, even as a weak dependency.
Production applications can pass `current_levelset(geom)` and `grid_info(geom.grid).dims` to a renderer through ordinary arrays.
The existing optional CairoMakie extension remains separate from terminal plotting.

In the scientific-stack workspace, `Stack/examples/terminal.jl` supplies visual entry points for these same numerical examples.
The visual entry points use the shared simulation code, not copies of the solvers.
Run `julia setup_stack.jl` from the stack root with Julia 1.11 or later to create that separate environment.
See the workspace `STACK.md` for commands and registration notes. No setup script is required for ED alone.

If you are interested in this work please feel free to contact me at: `naoufel.cresson@inria.fr`

my personal page: https://www.ljll.fr/~cresson/

# License
MIT License — see [LICENSE](LICENSE) for details.
