---
name: evolvingdomains
description: Guides writing and reviewing Julia simulations with EvolvingDomains.jl. Applies to moving-domain and level-set problems, external PDE-solver coupling, field transport, and audits of API correctness, snapshot ordering, storage aliasing, and numerical evidence.
---

# Verifiable evolving-domain simulations

Make the program's mathematics, data flow, mutations, and checks visible. Use
EvolvingDomains as shared data and numerical primitives, not as a hidden simulation engine.
These instructions guide implementation and review; they do not certify scientific accuracy.

## 1. Establish the problem before choosing methods

Identify the domain, dimension, units, level-set sign convention, physical state variables,
initial/boundary conditions, governing equations or prescribed motion, and final time.
Distinguish intensive quantities, densities, integrated quantities, and geometric fields.
For `cut`, the interior is the negative phase of the level set; convert other conventions
explicitly rather than silently reversing the domain.

Ask about missing information that changes the model. Do not invent a constitutive law,
boundary condition, external solver, or conservation claim to complete a runnable example.
For an underspecified coupling, provide an explicit interface and identify what remains
undefined. For prescribed motion, say that no coupled physical dynamics are being solved.

## 2. Check the installed public API

Before running, check that Julia and the project's dependencies are available. If execution
or native cutter loading is unavailable, report a blocked check rather than a pass.
Use the active Julia project's EvolvingDomains version. Inspect exports, method signatures,
and docstrings before writing calls; locate its source with `pkgdir(EvolvingDomains)` when
needed. `src/EvolvingDomains.jl` is the export declaration site relative to that package root.
Do not depend on private bindings, underscored helpers, or qualified internal kernels.

The README's opening simulation loop is schematic, not executable API. In particular:

- Use `Field(discretization, Dict(level_id => values))`, not a guessed `Field(grid, values)`.
- The public transport names are `eulerianTransport!` and `lagrangianTransport!`.
- `dynamics`, `velocityFromRates`, `evolve_state`, and `update` are application functions,
  not EvolvingDomains exports. Define any such function you use or mark it as unimplemented.
- The exported cut-cell accessor is `cutCellIds`, not `cut_cell_ids`.

Every helper in runnable code must be defined or imported; do not invent a helper name to
stand in for an API you have not checked. Do not change the library API to fit the schematic.
Use selective or qualified imports when combining toolkits: Ferrite and EvolvingDomains
both export `Grid`, for example. An import that succeeds does not prove an unqualified
constructor or callback is callable; probe the actual call and its return representation.

## 3. Keep ownership and representation explicit

- Build `Grid`, `Discretization`, and separate geometry/state `Field`s. Use named `Symbol`
  levels explicitly rather than dictionary iteration order or an accidental default.
- `Field` retains its input arrays. Distinct field containers do not guarantee distinct
  storage; use independent source/target arrays for transport and swap their roles afterwards.
- Reuse the same grid identity for fields, cuts, and maps that belong together. Equal
  coordinates in a newly constructed `Grid` do not give it the original mesh identity.
- Do not mutate grid axes underneath existing cuts, maps, or caches. An intentional grid
  change requires explicit data transfer and rebuilding dependent objects.
- Distinguish the general data containers from supported numerical methods: for the combined
  cut/evolution workflow, use finite dense nodal `Float64` arrays on uniform 2D or 3D grids.
- For `eulerianTransport!`, supply one `SVector{D,Float64}` per node in a vector matching
  `vec(values)` order. Facet data or normal speeds are not already that nodal vector field;
  any conversion/extension and its normal orientation must be explicit.

## 4. Write an explicit simulation loop

Keep initialization, parameters, application dynamics, stepping, and diagnostics separate.
A small named tuple or ordinary variables are sufficient; do not add an engine, callback
registry, or configuration framework merely to wrap the loop. Prefer plain functions for
solver adapters; do not introduce unused solver-wrapper types.

Expose these responsibilities in the caller:

1. Read the state and geometry at a declared time; evaluate the prescribed dynamics or
   external solver. Make any mutation by that solver explicit.
2. Obtain the velocity representation needed by the chosen transport method.
3. Choose a finite, positive stepping interval and clip the last step to the remaining
   time. The package does not choose a stable timestep for you.
4. Advance geometry and any physical state in the declared splitting/time-level order.
5. Perform the checks relevant to this step and advance the recorded time consistently.

Do not hide state updates in plotting, output, diagnostics, or global variables. Put the
loop in a function to avoid Julia global soft-scope surprises. Reuse valid work buffers and
caches outside the loop, without reusing a map whose geometry, velocity, or timestep changed.

## 5. Choose the state-transport contract deliberately

| Need | Public operation | What not to assume |
| --- | --- | --- |
| Advect a nodal level set | `eulerianTransport!(geometry, velocity, dt; level_id=...)` | It does not also update physical state or choose `dt`. |
| Transport an intensive nodal field using CIP | `lagrangianTransport!(target, source, velocity, dt; level_id=...)` | This is not the conservative cut-to-cut map method. |
| Redistribute a nodal state between cut supports | `TransportMap(previous, current, grid, velocity, dt)` then `lagrangianTransport!(target, source, map)` | Its conditional unweighted nodal-sum balance is not physical-volume-weighted mass conservation. |

For the cut-to-cut path, the required order is:

```julia
# geometry is the level-set Field; source holds the transported state at t.
previous = cut(geometry, level_id)
eulerianTransport!(geometry, velocity, dt; level_id)
current = cut(geometry, level_id)
transport = TransportMap(previous, current, grid, velocity, dt)
lagrangianTransport!(target, source, transport)
source, target = target, source
```

This fragment assumes initialized, compatible objects and independent transport buffers;
it is not a complete simulation. Both snapshots must refer to the same grid/level, but to
the intended before/after geometry. A `Cut` owns a snapshot: it does not update when the
level set changes, and mesh identity alone cannot prove that its time level is correct.
Before the transport, `source` must contain the state at the time of `previous`; after the
swap it contains the state at the time of `current`. Historical snapshots remain useful;
rebuild a cut whenever current geometry is required.

A physical source term or external PDE solve is a separate, visible update with its own
boundary conditions and time level. Do not silently substitute one of these transport
methods for that solve.

## 6. Build external-solver bridges in executable stages

Do not write an entire coupled application before testing the integration boundary.
If execution is delegated, hand off one minimal stage and await its observed result:

1. In a separate application environment, record the solver versions and run a minimal
   backend solve. Exercise callbacks with their real argument/return shapes; finding a
   function name in source is not execution evidence.
2. Test data exchange independently: state where backend DoFs live, how the solved field
   is sampled/projected, and how coordinates, components and `vec` order map to Field data.
   Use a known field on different, non-square grids, including boundary/interior probes.
   Do not reshape an unrelated FE/DG storage vector onto the EvolvingDomains grid.
3. Run one coupled step, then a short multi-step case. Carry evolved numerical state across
   steps rather than resetting it to analytic data. For time integrators, check the return
   code and attained horizon; for algebraic solves, check a residual. Make the time of each
   exchanged quantity explicit, including diagnostics labelled final.
4. Check both directions of coupling where the model requires them: geometry must actually
   affect the solver, and the computed solver output must actually drive geometry. An import
   test, a frozen coefficient, or an analytic substitute for the solver is not that check.

Keep authoritative solver-native state where needed; use explicit adapters rather than
advecting a pressure solve or restarting an evolved DG state from its analytic initial data.
Prefer documented interpolation APIs. Isolate any unavoidable version-specific storage
access in the adapter, assert its layout, and report the tested version and limitation.
If a stage fails, repair that stage and rerun it before composing more physics. Do not
replace the requested backend or relax a tolerance silently. Distinguish a background-box
PDE coupled to a cut from a PDE actually assembled on the moving/cut domain.

## 7. Keep geometric corrections opt-in

Do not insert `reinitialize!`, `filterUnresolvedFeatures!`, aggregation, or extension as
unexplained boilerplate. State why an operation is needed, where it belongs in the splitting,
and which snapshot it consumes. Rebuild subsequent geometry-dependent data when necessary.
Filtering requires approximately signed-distance input and does not reinitialize it. It can
remove genuine sub-resolution features and does not conserve area or volume; it is not an
accuracy fix. Do not catch numerical errors merely to continue with stale data.

## 8. Make the result checkable

Before an expensive run, execute a small case with assertions appropriate to the model:

- Check shapes, finite values, selected levels, independent transport storage, and progress
  of the time loop. Verify that historical cuts remain unchanged when fields evolve.
- Check a stationary or prescribed-motion case where the expected behavior is known.
- Define any balance using its actual discrete quantity, active support, sources, and
  boundary fluxes. The map conserves the unweighted nodal sum only for contributions that
  reach the target support; do not assert global physical mass conservation from that alone.
- For a no-loss case like the example, capture
  `expected = sum(vec(source.data[level_id])[transport.source_mask])` before transport, then
  compare with `sum(target.data[level_id])` before swapping buffers. Read `transport.source_mask`
  for this diagnostic; do not modify it or invent a support-mask helper. If contributions can
  leave the target support, establish the loss/flux budget instead of assuming this equality.
- If claiming accuracy, separately measure errors against a known solution and refine space
  and time. A finite output, conservation check, or passing API test is not accuracy evidence.

Start from [the runnable prescribed-motion example](examples/prescribed_motion.jl), not the
README pseudocode. From the root of an instantiated EvolvingDomains checkout, run:

```sh
julia --startup-file=no --project=. skills/evolvingdomains/examples/prescribed_motion.jl
julia --startup-file=no --project=. test/run_conformance.jl
```

Success prints `Prescribed-motion smoke checks passed`, followed by steps, time, and nodal sum.
From a separate application project, use its environment and the actual path to the example;
that project must declare EvolvingDomains and StaticArrays as dependencies. Published cutter
artifacts cover x86_64/aarch64 on Linux (glibc ≥ 2.28) or macOS ≥ 13; not Windows or musl Linux.
Do not override library paths or add a source build by default. The conformance command is
for the library checkout, not a requirement to copy its test suite into each application.

## 9. Deliver evidence, not assurances

Provide the runnable entry point and exact command; explain the physical inputs, numerical
choices, state/geometry update order, and diagnostics. Report checks actually executed,
observed results, and unresolved assumptions. Clearly label unimplemented coupling code,
untested paths, and unsupported accuracy or conservation claims. Never report a test as
passing because it appears likely to pass.
