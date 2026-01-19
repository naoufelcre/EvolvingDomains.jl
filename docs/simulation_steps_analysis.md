# Moving Domain Simulation: Step-by-Step Analysis

This document breaks down the lifecycle of a moving domain simulation, analyzing the contribution of `EvolvingDomains.jl` at each stage. We distinguish between **Logic** (workflow simplification) and **Performance** (optimizations).

---

## 1. Geometry Initialization
**The Problem**: Defining complex initial shapes and generating a computational mesh that respects these shapes.

### Logic (Abstraction)
- **CSG & SDFs**: Provides a high-level API for Constructive Solid Geometry (CSG) (Union, Intersection, Difference) to compose shapes intuitively.
- **Dual-Model Architecture**: Hides the complexity of managing two distinct meshes:
    1.  **Physics Mesh**: Unstructured, adapted to the interface.
    2.  **LSM Mesh**: Euclidean Cartesian grid for Level Set evolution.

### Performance (Optimization)
- **Quadtree Generation**: Generates adaptive background meshes rapidly ($O(N \log N)$), refining only where necessary (near the interface).
- **Memory Economy**: Reduces cell count by orders of magnitude (e.g., 48x reduction) compared to uniform grids, minimizing memory footprint for the subsequent physics solve.

---

## 2. Velocity Field Handling
**The Problem**: The domain moves according to a velocity field $v$, which may be analytical or computed from a physics solver (FSI).

### Logic (Abstraction)
- **Unified Interface**: Treats analytical functions (`t -> v(t)`) and Finite Element fields (`FEFunction`) interchangeably via `AbstractVelocity`.
- **FSI Coupling**: `FEVelocitySource` handles the tricky mapping of separate FE solutions onto the geometry interaction interface.

### Performance (Optimization)
- **Cached Extension**: When extending velocity from the fluid domain to the void (needed for level set transport), it uses cached geometric queries (kd-trees) to avoid recomputing nearest neighbors every time step.

---

## 3. Interface Evolution (Advection)
**The Problem**: Accurately transporting the interface $\Gamma(t)$ based on the velocity field.

### Logic (Abstraction)
- **`advance!(geo, dt)`**: A single call handles the entire advection step, properly coordinating the level set update on the background LSM mesh.

### Performance (Optimization)
- **Narrow Band (External)**: Relies on `LevelSetMethods.jl` optimized schemes (Runge-Kutta, WENO) which operate effectively on the supporting Cartesian grid.

---

## 4. Mesh Adaptation (Remeshing)
**The Problem**: As the domain moves, the mesh must adapt to the new configuration to maintain accuracy and prevent distortion.

### Logic (Abstraction)
- **Automatic Trigger**: The package monitors mesh quality/distortion and automatically triggers remeshing when limits are exceeded (though often done every step for cut-cell methods).
- **Template Paving**: Converts the non-conforming Quadtree (with hanging nodes) into a strictly conforming triangular mesh valid for standard FEM solvers.

### Performance (Optimization)
- **Local Updates**: (Planned/Partial) Ideally only updates regions near the moving interface. Currently, the Quadtree regeneration is fast enough that full regeneration is often viable.
- **2:1 Balancing**: Ensures gradual grading of element sizes, preventing numerical ill-conditioning in the physics solve.

---

## 5. Physical Solve (CutFEM/AgFEM)
**The Problem**: Solving the PDE on the complex, moving domain.

### Logic (Abstraction)
- **Gridap Integration**: Seamlessly exposes the adapted mesh as a standard `Gridap` discrete model.
- **AgFEM Support**: Provides necessary topology for Aggregated Finite Elements to stabilize small cut cells.

### Performance (Optimization)
- **DOF Reduction**: The adaptive Quadtree mesh places Degrees of Freedom (DoFs) only where physics demands them (the interface), leading to massive speedups in linear system assembly and solve times (e.g., ~25x faster solves).
