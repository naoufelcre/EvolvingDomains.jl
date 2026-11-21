# Downstream StaticAssembler integration

**Verdict: accepted** for the real downstream quadrature/assembly path.

```sh
EVOLVINGDOMAINS_CUTCELLS_LIBRARY=/absolute/path/to/library \
STATICASSEMBLER_PATH=/path/to/StaticAssembler.jl/src \
julia --project=validation/env --startup-file=no validation/integration/static_assembler/run.jl
```

`STATICASSEMBLER_PATH` is the directory that contains `StaticAssembler.jl`, i.e. the
downstream package's `src/`.

The run exercises real quadrature/normals/area/assembly for cut, mixed cut/full/outside and
aligned cases. The native CutCells backend is the only producer; the GridapEmbedded producer
was removed with the local descriptor conversion replacing it. Recorded result: **52 checks
passed** (16, 16 and 20 across the three cases).

**Limits:** downstream version and environment are external; this is not a 3D assembly
contract. The 3D cut currently exposes unwelded volume/interface vertices.
