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
    # 1. Pre-generate locally to ensure artifacts exist
    # usage of @__DIR__ in scripts implies they save to their execution location.
    # We execute them in their source directory to generate the GIFs/PNGs there.
    println("Pre-generating artifact for $file ...")
    cmd = `julia --project=docs $(joinpath(examples_dir, file))`
    run(cmd)

    # 2. Copy artifacts to docs/src/examples
    # Identify generated file (basename of script + .gif or .png usually)
    base = splitext(file)[1]
    for ext in [".gif", ".png"]
        artifact_src = joinpath(examples_dir, base * ext)
        if isfile(artifact_src)
            cp(artifact_src, joinpath(output_dir, base * ext); force=true)
            println("Copied $base$ext to output directory.")
        end
    end

    # 3. Process with Literate
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
