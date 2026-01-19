# Changelog

All notable changes to this project will be documented in this file.

## [Unreleased] - 2026-01-16

### Added
- **Benchmark**: New script `benchmark_poisson_quadtree.jl` comparing `EvolvingDomains` Quadtree meshing vs standard Uniform meshes for a Poisson problem (Crescent Moon geometry).
  - Demonstrates ~25-35x speedup in solver time for fine resolutions.
  - Includes VTK visualization output in `output/` directory.

### Fixed
- **Robustness**: Identified and documented a **Grid Alignment Singularity** at Level 6 ($L=6$) which caused `AgFEM` aggregation failure (`AssertionError: all_aggregated`).
  - **Cause**: Perfect alignment of the interface with background grid lines led to ambiguous topology and "sliver" cuts.
  - **Resolution**: A micro-perturbation of the domain origin ($\delta \approx 1e^{-5}$) resolves the singularity.
  - **Technical Detail**: See `docs/benchmark_report.md` for a mathematical explanation of the connectivity graph failure in degenerate alignment configurations.
