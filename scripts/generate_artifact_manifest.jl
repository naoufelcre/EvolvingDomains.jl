using Pkg.Artifacts
using SHA
using Base.BinaryPlatforms

length(ARGS) == 3 || error("usage: julia scripts/generate_artifact_manifest.jl ARCHIVE_DIR RELEASE_TAG OUTPUT_TOML")
archive_dir, release_tag, output_toml = ARGS
occursin(r"^bindings-cutcells-v0\.4\.0-[1-9][0-9]*$", release_tag) ||
    error("release tag must match bindings-cutcells-v0.4.0-N")

const TARGETS = (
    ("x86_64-linux-gnu", "x86_64-linux-gnu"),
    ("aarch64-linux-gnu", "aarch64-linux-gnu"),
    ("x86_64-apple-darwin", "x86_64-apple-darwin"),
    ("aarch64-apple-darwin", "aarch64-apple-darwin"),
    ("x86_64-w64-mingw32", "x86_64-w64-mingw32"),
)

ispath(output_toml) && error("refusing to replace existing artifact metadata: $output_toml")
for (label, _) in TARGETS
    isfile(joinpath(archive_dir, "evolvingdomains-cutcells-$(label).tar.gz")) ||
        error("missing platform archive: $label")
    isdir(joinpath(archive_dir, "unpacked", label)) || error("missing extracted artifact tree: $label")
end
mkpath(dirname(abspath(output_toml)))

for (label, triplet) in TARGETS
    archive = joinpath(archive_dir, "evolvingdomains-cutcells-$(label).tar.gz")
    unpacked = joinpath(archive_dir, "unpacked", label)
    isfile(archive) || error("missing platform archive: $archive")
    isdir(unpacked) || error("missing extracted artifact tree: $unpacked")

    tree_hash = create_artifact() do artifact_dir
        for (root, dirs, files) in walkdir(unpacked)
            relative = relpath(root, unpacked)
            destination = relative == "." ? artifact_dir : joinpath(artifact_dir, relative)
            mkpath(destination)
            for dir in dirs
                mkpath(joinpath(destination, dir))
            end
            for file in files
                cp(joinpath(root, file), joinpath(destination, file); force=true)
            end
        end
    end
    archive_sha256 = bytes2hex(open(sha256, archive))
    repository = get(ENV, "GITHUB_REPOSITORY", "naoufelcre/EvolvingDomains.jl")
    url = "https://github.com/$(repository)/releases/download/$(release_tag)/$(basename(archive))"
    platform = if endswith(triplet, "-linux-gnu")
        arch = startswith(triplet, "aarch64") ? "aarch64" : "x86_64"
        Platform(arch, "linux"; libc="glibc")
    elseif endswith(triplet, "-apple-darwin")
        arch = startswith(triplet, "aarch64") ? "aarch64" : "x86_64"
        # BinaryPlatforms records Darwin versions: macOS 13 corresponds to Darwin 22.
        Platform(arch, "macos"; os_version=v"22.0")
    else
        parse(Platform, triplet)
    end
    bind_artifact!(output_toml, "cutcells_bindings", tree_hash;
        platform,
        download_info=[(url, archive_sha256)],
        lazy=false, force=true)
end

println("Wrote platform-specific native artifact metadata to $output_toml")
