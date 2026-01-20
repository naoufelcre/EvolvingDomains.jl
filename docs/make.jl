using EvolvingDomains
using Documenter

using Literate

# Generate examples
examples_dir = joinpath(@__DIR__, "..", "examples")
output_dir = joinpath(@__DIR__, "src", "examples")

# Create output directory if it doesn't exist
mkpath(output_dir)

# List of examples to process
example_files = [
    "colliding_balls_minimal.jl",
    "csg.jl",
    "independent_motion_minimal.jl",
    "stokes_driven_minimal.jl",
    "zalesak_minimal.jl",
]

generated_examples = []

for file in example_files
    input_path = joinpath(examples_dir, file)
    Literate.markdown(input_path, output_dir; documenter=true)
    push!(generated_examples, "examples/$(replace(file, ".jl" => ".md"))")
end

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
        "Examples" => generated_examples,
    ],
)

deploydocs(;
    repo="github.com/naoufelcre/EvolvingDomains.jl",
    devbranch="main",
)
