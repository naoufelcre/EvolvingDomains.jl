using Libdl

@testset "Native library discovery" begin
    mktempdir() do root
        cp(joinpath(@__DIR__, "..", "..", "src"), joinpath(root, "src"))
        sandbox = Module(gensym(:NativeLibraryFixture))
        Base.include(sandbox, joinpath(root, "src", "EvolvingDomains.jl"))
        library_path = () -> Base.invokelatest(() -> sandbox.EvolvingDomains.CutCellsBinding._libraryPath())
        filename = (Sys.iswindows() ? "" : "lib") * "evolvingdomains_cut_cells." * Libdl.dlext
        local_library = joinpath(root, "deps", "usr", "lib", filename)

        withenv("EVOLVINGDOMAINS_CUTCELLS_LIBRARY" => nothing) do
            error = try
                library_path()
            catch err
                err
            end
            @test error isa ArgumentError
            @test occursin("build_cxx_bindings.sh", sprint(showerror, error))
            @test occursin("found automatically", sprint(showerror, error))

            mkpath(dirname(local_library))
            touch(local_library) # Discovery must not load native code in the contract suite.
            @test library_path() == local_library
            withenv("EVOLVINGDOMAINS_CUTCELLS_LIBRARY" => "") do
                @test library_path() == local_library
            end

            override = joinpath(root, filename)
            touch(override)
            withenv("EVOLVINGDOMAINS_CUTCELLS_LIBRARY" => override) do
                @test library_path() == override
            end
            withenv("EVOLVINGDOMAINS_CUTCELLS_LIBRARY" => joinpath(root, "missing")) do
                @test_throws ArgumentError library_path()
            end
            withenv("EVOLVINGDOMAINS_CUTCELLS_LIBRARY" => root) do
                @test_throws ArgumentError library_path()
            end
        end
    end
end
