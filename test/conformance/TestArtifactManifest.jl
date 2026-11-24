using Test
using TOML
using SHA
using Pkg.Artifacts
using Base.BinaryPlatforms

@testset "Native artifact manifest" begin
    mktempdir() do directory
        labels = ("x86_64-linux-gnu", "aarch64-linux-gnu", "x86_64-apple-darwin",
                  "aarch64-apple-darwin")
        for label in labels
            tree = joinpath(directory, "unpacked", label)
            mkpath(tree)
            # These are manifest fixtures, not platform binaries.
            write(joinpath(tree, "fixture.txt"), label)
            run(`tar -czf $(joinpath(directory, "evolvingdomains-cutcells-$label.tar.gz")) -C $tree .`)
            write(joinpath(tree, "not-in-archive.txt"), "staging directories are not authoritative")
        end
        output = joinpath(directory, "Artifacts.toml")
        generator = joinpath(@__DIR__, "..", "..", "scripts", "generate_artifact_manifest.jl")
        command = `$(Base.julia_cmd()) --startup-file=no $generator $directory cut-cells-binaries-1 $output`
        run(addenv(command, "GITHUB_REPOSITORY" => "example/EvolvingDomains.jl"))
        entries = TOML.parsefile(output)["cutcells_bindings"]
        @test length(entries) == 4
        for entry in entries
            label = only(filter(l -> occursin(l, only(entry["download"])["url"]), labels))
            archive = joinpath(directory, "evolvingdomains-cutcells-$label.tar.gz")
            @test only(entry["download"])["sha256"] == bytes2hex(open(sha256, archive))
            @test get(entry, "lazy", false) == false
            @test only(entry["download"])["url"] ==
                "https://github.com/example/EvolvingDomains.jl/releases/download/cut-cells-binaries-1/evolvingdomains-cutcells-$label.tar.gz"
            @test isfile(joinpath(artifact_path(Base.SHA1(entry["git-tree-sha1"])), "fixture.txt"))
            @test !isfile(joinpath(artifact_path(Base.SHA1(entry["git-tree-sha1"])), "not-in-archive.txt"))
            if entry["os"] == "macos"
                @test VersionNumber(entry["os_version"]) == v"22"
            elseif entry["os"] == "linux"
                @test entry["libc"] == "glibc"
                @test !haskey(entry, "libc_version")
            end
        end
        @test artifact_hash("cutcells_bindings", output; platform=Platform("x86_64", "windows")) === nothing
        @test artifact_hash("cutcells_bindings", output; platform=Platform("aarch64", "linux"; libc="musl")) === nothing
        @test artifact_hash("cutcells_bindings", output; platform=HostPlatform(Platform("x86_64", "macos"; os_version=v"21"))) === nothing
        @test artifact_hash("cutcells_bindings", output; platform=HostPlatform(Platform("x86_64", "macos"; os_version=v"22"))) !== nothing
        before = read(output)
        @test !success(pipeline(command; stdout=devnull, stderr=devnull))
        @test read(output) == before
        for tag in ("../unsafe", "release?query", "")
            rejected_output = joinpath(directory, "rejected.toml")
            rejected = `$(Base.julia_cmd()) --startup-file=no $generator $directory $tag $rejected_output`
            @test !success(pipeline(rejected; stdout=devnull, stderr=devnull))
            @test !ispath(rejected_output)
        end
    end
end

@testset "Published package artifact manifest" begin
    manifest = joinpath(@__DIR__, "..", "..", "Artifacts.toml")
    entries = TOML.parsefile(manifest)["cutcells_bindings"]
    @test length(entries) == 4
    @test Set((entry["arch"], entry["os"]) for entry in entries) ==
        Set((arch, os) for arch in ("aarch64", "x86_64"), os in ("linux", "macos"))
    for entry in entries
        platform = entry["os"] == "linux" ? Platform(entry["arch"], "linux"; libc="glibc") :
            HostPlatform(Platform(entry["arch"], "macos"; os_version=v"22"))
        @test string(artifact_hash("cutcells_bindings", manifest; platform)) == entry["git-tree-sha1"]
        label = entry["arch"] * (entry["os"] == "linux" ? "-linux-gnu" : "-apple-darwin")
        download = only(entry["download"])
        @test download["url"] ==
            "https://github.com/naoufelcre/EvolvingDomains.jl/releases/download/cut-cells-binaries-1/evolvingdomains-cutcells-$label.tar.gz"
        @test occursin(r"^[0-9a-f]{64}$", download["sha256"])
        @test !get(entry, "lazy", false)
    end
    @test artifact_hash("cutcells_bindings", manifest; platform=Platform("x86_64", "windows")) === nothing
    @test artifact_hash("cutcells_bindings", manifest; platform=Platform("aarch64", "linux"; libc="musl")) === nothing
    @test artifact_hash("cutcells_bindings", manifest; platform=HostPlatform(Platform("x86_64", "macos"; os_version=v"21"))) === nothing
end
