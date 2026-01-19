# EvolvingDomains V2

**A minimal Julia package for solving PDEs on moving domains.**

EvolvingDomains.jl simplifies workflows within the Gridap Ecosystem for problems where the domain evolves with respect to the PDE solution (FSI) or prescribed velocities. It leverages `GridapEmbedded` and `LevelSetMethods.jl` to handle the complexities of geometry handling, mesh adaptation, and velocity coupling.

## Key Features
- **External Mesh Support**: Integrate with external mesh generators via generic locator interface.
- **Dual-Model Architecture**: Seamlessly handles a Physics Mesh (unstructured) and a Level Set Mesh (Cartesian).
- **Velocity Coupling**: Unified interface for analytical and FE-based velocity fields.

## Documentation
- [Workflows](examples/WORKFLOW.md): Step-by-step guide to setting up a simulation.




## License
MIT License — see [LICENSE](LICENSE) for details.
