using Test
using TOML
using SHA
using Pkg.Artifacts
using Base.BinaryPlatforms

@testset "Native artifact manifest" begin
    mktempdir() do directory
        labels = ("x86_64-linux-gnu", "aarch64-linux-gnu", "x86_64-apple-darwin",
                  "aarch64-apple-darwin", "x86_64-w64-mingw32")
        for label in labels
            tree = joinpath(directory, "unpacked", label)
            mkpath(tree)
            # These are manifest fixtures, not platform binaries.
            write(joinpath(tree, "fixture.txt"), label)
            run(`tar -czf $(joinpath(directory, "evolvingdomains-cutcells-$label.tar.gz")) -C $tree .`)
        end
        output = joinpath(directory, "Artifacts.toml")
        generator = joinpath(@__DIR__, "..", "..", "scripts", "generate_artifact_manifest.jl")
        command = `$(Base.julia_cmd()) --startup-file=no $generator $directory bindings-cutcells-v0.4.0-1 $output`
        run(addenv(command, "GITHUB_REPOSITORY" => "example/EvolvingDomains.jl"))
        entries = TOML.parsefile(output)["cutcells_bindings"]
        @test length(entries) == 5
        for entry in entries
            label = only(filter(l -> occursin(l, only(entry["download"])["url"]), labels))
            archive = joinpath(directory, "evolvingdomains-cutcells-$label.tar.gz")
            @test only(entry["download"])["sha256"] == bytes2hex(open(sha256, archive))
            @test get(entry, "lazy", false) == false
            @test startswith(only(entry["download"])["url"], "https://github.com/example/EvolvingDomains.jl/")
            @test isfile(joinpath(artifact_path(Base.SHA1(entry["git-tree-sha1"])), "fixture.txt"))
            if entry["os"] == "macos"
                @test VersionNumber(entry["os_version"]) == v"22"
            elseif entry["os"] == "linux"
                @test entry["libc"] == "glibc"
                @test !haskey(entry, "libc_version")
            end
        end
        @test artifact_hash("cutcells_bindings", output; platform=Platform("aarch64", "linux"; libc="musl")) === nothing
        @test artifact_hash("cutcells_bindings", output; platform=HostPlatform(Platform("x86_64", "macos"; os_version=v"21"))) === nothing
        @test artifact_hash("cutcells_bindings", output; platform=HostPlatform(Platform("x86_64", "macos"; os_version=v"22"))) !== nothing
        before = read(output)
        @test !success(pipeline(command; stdout=devnull, stderr=devnull))
        @test read(output) == before
    end
end
