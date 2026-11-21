# Original-package numerical regressions

**Verdict: baseline evidence**, separate from the redesign contract suite. This directory is
the porting pool for future cases.

```sh
julia --startup-file=no validation/archive/legacy_numerics/regressions.jl
```

This runner resolves the active environment's `EvolvingDomains` dev package and preserves
the established WENO, transport, geometry, aggregation, curvature, reinitialization and
transfer checks. Selected checks passed in the earlier local baseline. It is not a test of
the redesign API and must not be run in the redesigned environment and reported as such.

**Limits:** runtime and environment depend on the configured original package; large
plotting/scenario scripts are not part of the baseline.

## What each file is for

| File | Target once ported |
|---|---|
| `regressions.jl` | the ported entry point; the list below in order |
| `TestHeleShawST.jl`, `TestHeleShawDispersion.jl`, `examples/heleShaw*.jl` | `hele_shaw/` — curvature and surface tension |
| `examples/heatDumbell.jl`, `TestDumbellParabolic.jl` | `heat_dumbbell/` — parabolic diffusion on an evolving domain |
| `TestCurvature.jl`, `TestReinitialization.jl` | curvature/reinitialization case, or `test/` contracts |
| `TestWENO5.jl`, `TestConservativeTransport.jl`, `TestIntensiveTransport.jl` | covered by `test/` contracts; keep as baseline only |
| `TestBilinearInterpolation.jl`, `TestTransportWeightCache.jl`, `TestMultiComponentTransfer.jl`, `TestVisualTransfer.jl`, `TestAggregation.jl`, `TestGeometryEvolution.jl` | helper coverage; port into `test/` if still meaningful |
| `Transfer/` | grid/mesh transfer operators; port into `test/` if the new API exposes an equivalent |
| `Helpers/`, `examples/`, `PublicationFigures.jl` | old plotting/aggregation scaffolding, not part of the baseline; delete when the last porting target is done |

A file leaves this directory only when its replacement runs on the current API.
