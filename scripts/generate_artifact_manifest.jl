using Pkg.Artifacts
using Pkg: PlatformEngines
using SHA
using Base.BinaryPlatforms

length(ARGS) == 3 || error("usage: julia scripts/generate_artifact_manifest.jl ARCHIVE_DIR RELEASE_TAG OUTPUT_TOML")
archive_dir, release_tag, output_toml = ARGS
occursin(r"^[A-Za-z0-9][A-Za-z0-9._-]*$", release_tag) ||
    error("release tag must start with a letter or digit and contain only letters, digits, dots, underscores, or hyphens")

const TARGETS = (
    ("x86_64-linux-gnu", "x86_64-linux-gnu"),
    ("aarch64-linux-gnu", "aarch64-linux-gnu"),
    ("x86_64-apple-darwin", "x86_64-apple-darwin"),
    ("aarch64-apple-darwin", "aarch64-apple-darwin"),
)

ispath(output_toml) && error("refusing to replace existing artifact metadata: $output_toml")
for (label, _) in TARGETS
    isfile(joinpath(archive_dir, "evolvingdomains-cutcells-$(label).tar.gz")) ||
        error("missing platform archive: $label")
end
mkpath(dirname(abspath(output_toml)))

for (label, triplet) in TARGETS
    archive = joinpath(archive_dir, "evolvingdomains-cutcells-$(label).tar.gz")
    tree_hash = create_artifact() do artifact_dir
        # Hash exactly the tree Pkg installs, not a separately extracted staging directory.
        PlatformEngines.unpack(archive, artifact_dir)
    end
    archive_sha256 = bytes2hex(open(sha256, archive))
    repository = get(ENV, "GITHUB_REPOSITORY", "naoufelcre/EvolvingDomains.jl")
    url = "https://github.com/$(repository)/releases/download/$(release_tag)/$(basename(archive))"
    platform = if endswith(triplet, "-linux-gnu")
        arch = startswith(triplet, "aarch64") ? "aarch64" : "x86_64"
        Platform(arch, "linux"; libc="glibc")
    else
        arch = startswith(triplet, "aarch64") ? "aarch64" : "x86_64"
        # BinaryPlatforms records Darwin versions: macOS 13 corresponds to Darwin 22.
        Platform(arch, "macos"; os_version=v"22.0")
    end
    bind_artifact!(output_toml, "cutcells_bindings", tree_hash;
        platform,
        download_info=[(url, archive_sha256)],
        lazy=false, force=true)
end

println("Wrote platform-specific native artifact metadata to $output_toml")
