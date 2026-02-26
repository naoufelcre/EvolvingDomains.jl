# EvolvingDomains 0.2

**A Julia package for solving PDEs on moving domains.**

Mini package intended to provide a set of utilities to write 2D moving domain problems in the gridap ecosystem. It leverages `GridapEmbedded`.

The paradigm the package is built on is a decoupling between kinematics and dynamics. In particular the package provides tools to handle the kinematics side by dedicated geometric structures. The dynamics part is intended to be handled on the active mesh by FEM solving with `GridapEmbedded`.

For a full working example see `TestDumbellParabolic.jl` that recreates a test case from the 2025 paper of Olshanskii & Reusken. [arXiv:2504.14116](https://arxiv.org/pdf/2504.14116)

![Temperature evolution - Olshanskii & Reusken test case](TEMPERATURE_OLSHANSKII_REUSKEN.gif)

## Geometric

The basic object of the package is the `EvolvingDiscreteGeometry`. It is an all-in-one object for basic routines regarding implicitly defined level-set geometry that evolves.

In particular it provides the following functionalities:

- **Reinitialization of the level set** to a signed distance function.

  After advection the level-set gradient `|∇φ|` drifts away from 1. Reinitialization restores the signed distance property by solving the Eikonal equation `|∇φ| = 1` on the background grid.
  The implementation uses the **Fast Sweeping Method** (Zhao 2005): interface nodes are seeded first via the **Russo-Smereka** geometric subcell estimate (Russo & Smereka 2000) — a quadratic interpolation across sign-changing edges — then four alternating-direction sweeps propagate the distance outward.

  ```julia
  reinitialize!(geom)   # restores |∇φ| ≈ 1 everywhere
  ```

- **A lazy cache system** for stencils and derived data.

  `EvolvingDiscreteGeometry` holds a `GeometryCache` that stores the GridapEmbedded cut geometry, the active node set, and the transfer and extension operators. All entries are invalidated automatically whenever the level set changes (via `set_levelset!` or `reinitialize!`) and recomputed on first access.

  ```julia
  set_levelset!(geom, phi_new)        # updates φ, invalidates all cache entries
  indices = get_active_indices(geom)  # recomputes and caches active IN+CUT nodes
  ```
  **Active indices** are used to couple effectively with the CutFEM method provided by `GridapEmbedded`.


Many moving domain problems involve fields living on the geometry. To handle this the package provides a dedicated structure `CartesianMeshField`. It wraps the flat nodal data array and provides clamped 2D indexing and a bilinear interpolant (via `get_interpolator`), which is used internally by both the WENO5 stencils and the transfer operators.

## Kinematics

The package provides two distinct transport operators.

### Level-set advection — WENO5 + SSP-RK3

The first operator advances the level set by solving `∂φ/∂t + v·∇φ = 0`. 

The spatial discretization uses the **fifth-order WENO** scheme (Jiang & Shu 1996) with Jiang-Peng smoothness indicators (Jiang & Peng 2000): at each node the upwind-biased directional derivative is selected based on the sign of `v`, and non-linear weights suppress oscillations near discontinuities while recovering fifth-order accuracy on smooth regions. 

Time integration uses the **third-order Strong-Stability-Preserving Runge-Kutta** (SSP-RK3, Shu & Osher 1988).

It can be used with any well-defined velocity field, sampled onto the grid via `sample_velocity`:

```julia
vel = StaticFunctionVelocity(x -> VectorValue(-ω*(x[2]-0.5), ω*(x[1]-0.5)))
v_field = sample_velocity(vel, grid_info(grid), t)
advance!(geom, v_field, Δt)   # WENO5 + SSP-RK3 step on the level set
```

### Field advection — Conservative Semi-Lagrangian (CCISL)

The second operator advects fields coupled to the deforming geometry. It implements the **Conservative Cell-Integrated Semi-Lagrangian** method (Lentine, Grétarsson & Fedkiw 2011). 

```julia
k_map = TransportMap(geom, vel, Δt)           # build the transport map (geometry-aware)
advect!(new_data, current_field.data, k_map)  # apply it to any scalar field
```

The main drawback of this method is significant numerical diffusion. See the rotating checkerboard test `TestConservativeTransport.jl` for a illustration.

## Transfer

The most critical capability in hybrid workflows (also found in multigrid methods) is accurate transfer between different levels of discretization. The package provides a `GridMeshTransfer` operator that follows the `TransferOperator.jl` protocol, exposing two directions:

- **`restrict` (Grid → Mesh):** evaluates the `CartesianMeshField` at FE mesh nodes via bilinear interpolation, then projects into the target `FESpace` using Gridap's `interpolate`.
- **`prolong` (Mesh → Grid):** maps DOF values back onto the background grid using direct index mapping when the mesh topology allows it, falling back to batch point evaluation otherwise.

```julia
transfer_op = setup_transfer(geom, V)      # V is the current AgFEM FESpace
u_mesh = grid_to_mesh(geom, u_grid)        # restrict: Cartesian field → FE function
u_grid = mesh_to_grid(geom, u_mesh)        # prolong:  FE function  → Cartesian field
```

## License
MIT License — see [LICENSE](LICENSE) for details.
