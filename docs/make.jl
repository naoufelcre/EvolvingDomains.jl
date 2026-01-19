using EvolvingDomains
using Documenter

DocMeta.setdocmeta!(EvolvingDomains, :DocTestSetup, :(using EvolvingDomains); recursive=true)

makedocs(;
    modules=[EvolvingDomains],
    authors="Naoufel Cresson",
    sitename="EvolvingDomains.jl",
    format=Documenter.HTML(;
        canonical="https://naoufelcre.github.io/EvolvingDomains.jl",
        edit_link="main",
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
        "API Reference" => [
            "Geometry" => "geometry.md",
            "Velocity" => "velocity.md",
            "Internals" => "internals.md",
        ],
    ],
)

deploydocs(;
    repo="github.com/naoufelcre/EvolvingDomains.jl",
    devbranch="main",
)
