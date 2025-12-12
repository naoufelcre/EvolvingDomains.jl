using Documenter
using EvolvingDomains

# Generate documentation
makedocs(
    sitename = "EvolvingDomains.jl",
    authors = "Naoufel Cresson",
    modules = [EvolvingDomains],
    format = Documenter.HTML(
        prettyurls = get(ENV, "CI", nothing) == "true",
        canonical = "https://naoufelcre.github.io/EvolvingDomains.jl",
        assets = String[],
    ),
    pages = [
        "Home" => "index.md",
        "Getting Started" => "getting_started.md",
        "User Guide" => "user_guide.md",
        "API Reference" => [
            "Core Types" => "api/core.md",
            "Velocity Sources" => "api/velocity.md",
            "External Solver Interface" => "api/external.md",
        ],
        "Examples" => "examples.md",
    ],
    doctest = false,  # Disable doctests for now
    checkdocs = :none,  # Skip check for missing docstrings (TestVisualization is optional)
)

# Deploy to GitHub Pages (only in CI)
deploydocs(
    repo = "github.com/naoufelcre/EvolvingDomains.jl.git",
    devbranch = "main",
    push_preview = true,
)
