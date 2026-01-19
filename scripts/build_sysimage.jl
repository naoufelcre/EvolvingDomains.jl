# =============================================================================
# build_sysimage.jl — Create a precompiled sysimage for fast startup
# =============================================================================
#
# Usage:
#   julia --project=. scripts/build_sysimage.jl
#
# This will create 'evolving_domains.so' (Linux/Mac) or 'evolving_domains.dll' (Windows)
#
# To use the sysimage:
#   julia --sysimage=evolving_domains.so your_script.jl
#
# =============================================================================

using Pkg

# Ensure PackageCompiler is available
if !haskey(Pkg.project().dependencies, "PackageCompiler")
    @info "Installing PackageCompiler..."
    Pkg.add("PackageCompiler")
end

using PackageCompiler

# Configuration
const PACKAGES = [
    "EvolvingDomains",
    "Gridap",
    "GridapEmbedded",
    "LevelSetMethods",
    "StaticArrays",
]

const SYSIMAGE_NAME = Sys.iswindows() ? "evolving_domains.dll" : "evolving_domains.so"
const WARMUP_SCRIPT = joinpath(@__DIR__, "warmup_sysimage.jl")

# Check warmup script exists
if !isfile(WARMUP_SCRIPT)
    error("""
    Warmup script not found: $WARMUP_SCRIPT
    Please create scripts/warmup_sysimage.jl first.
    """)
end

@info "Building sysimage..." packages=PACKAGES output=SYSIMAGE_NAME

# Build the sysimage
create_sysimage(
    PACKAGES;
    sysimage_path = SYSIMAGE_NAME,
    precompile_execution_file = WARMUP_SCRIPT,
    cpu_target = "generic",  # For portability; use "native" for max perf on this machine
)

@info """
✅ Sysimage created successfully!

To use it:
    julia --sysimage=$(SYSIMAGE_NAME) your_simulation.jl

Or set as default in VS Code Julia extension:
    "julia.additionalArgs": ["--sysimage=/path/to/$(SYSIMAGE_NAME)"]
"""
