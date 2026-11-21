using Test
using Libdl

const _native_path = get(ENV, "EVOLVINGDOMAINS_CUTCELLS_LIBRARY", "")
isempty(_native_path) && error("set EVOLVINGDOMAINS_CUTCELLS_LIBRARY to the built native library")

lib = Libdl.dlopen(_native_path)
const _errcap = 512

function _native_error(buffer)
    first_zero = findfirst(iszero, buffer)
    bytes = first_zero === nothing ? buffer : @view(buffer[1:first_zero-1])
    return String(copy(bytes))
end

function _create_quad(coordinates, values)
    result = Ref{Ptr{Cvoid}}(C_NULL)
    err = zeros(UInt8, _errcap)
    status = GC.@preserve coordinates values err ccall(
        Libdl.dlsym(lib, :edcut_quad_create), Cint,
        (Ptr{Cdouble}, UInt64, Ptr{Cdouble}, UInt64, Ref{Ptr{Cvoid}}, Ptr{UInt8}, UInt64),
        pointer(coordinates), length(coordinates), pointer(values), length(values),
        result, pointer(err), _errcap)
    return status, result[], _native_error(err)
end

function _readRegion(result, region)
    points_count, offsets_count = Ref{UInt64}(), Ref{UInt64}()
    connectivity_count, types_count = Ref{UInt64}(), Ref{UInt64}()
    err = zeros(UInt8, _errcap)
    status = GC.@preserve err ccall(
        Libdl.dlsym(lib, :edcut_region_sizes), Cint,
        (Ptr{Cvoid}, Cint, Ref{UInt64}, Ref{UInt64}, Ref{UInt64}, Ref{UInt64}, Ptr{UInt8}, UInt64),
        result, region, points_count, offsets_count, connectivity_count, types_count,
        pointer(err), _errcap)
    @test status == 0
    @test isempty(_native_error(err))

    points = zeros(Float64, points_count[])
    offsets = zeros(Int32, offsets_count[])
    connectivity = zeros(Int32, connectivity_count[])
    types = zeros(Int32, types_count[])
    err .= 0
    status = GC.@preserve points offsets connectivity types err ccall(
        Libdl.dlsym(lib, :edcut_region_copy), Cint,
        (Ptr{Cvoid}, Cint, Ptr{Cdouble}, UInt64, Ptr{Int32}, UInt64,
         Ptr{Int32}, UInt64, Ptr{Int32}, UInt64, Ptr{UInt8}, UInt64),
        result, region, pointer(points), length(points), pointer(offsets), length(offsets),
        pointer(connectivity), length(connectivity), pointer(types), length(types),
        pointer(err), _errcap)
    @test status == 0
    @test isempty(_native_error(err))
    return points, offsets, connectivity, types
end

function _region_area(data)
    points, offsets, connectivity, types = data
    area = 0.0
    for cell in eachindex(types)
        ids = (Int(offsets[cell]) + 1):Int(offsets[cell + 1])
        @test length(ids) == 3
        xy = [(points[2i - 1], points[2i]) for i in (Int.(connectivity[ids]) .+ 1)]
        area += abs(sum(xy[i][1] * xy[mod1(i + 1, 3)][2] -
                        xy[mod1(i + 1, 3)][1] * xy[i][2] for i in 1:3)) / 2
    end
    return area
end

try
    @test unsafe_string(ccall(Libdl.dlsym(lib, :edcut_upstream_version), Cstring, ())) == "0.4.0"

    coordinates = [0.0, 0.0, 1.0, 0.0, 1.0, 1.0, 0.0, 1.0]
    values = [-1.0, 1.0, 1.0, -1.0]
    result = Ref{Ptr{Cvoid}}(C_NULL)
    err = zeros(UInt8, _errcap)
    status = GC.@preserve coordinates values err ccall(
        Libdl.dlsym(lib, :edcut_quad_create), Cint,
        (Ptr{Cdouble}, UInt64, Ptr{Cdouble}, UInt64, Ref{Ptr{Cvoid}}, Ptr{UInt8}, UInt64),
        pointer(coordinates), length(coordinates), pointer(values), length(values),
        result, pointer(err), _errcap)
    @test status == 0
    @test result[] != C_NULL
    @test isempty(_native_error(err))

    try
        inside = _readRegion(result[], 0)
        interface = _readRegion(result[], 1)
        outside = _readRegion(result[], 2)
        @test all(inside[4] .== 2) && all(outside[4] .== 2)
        @test interface[4] == [1]
        @test all((0.0 .<= inside[1][1:2:end]) .& (inside[1][1:2:end] .<= 0.5))
        @test all((0.5 .<= outside[1][1:2:end]) .& (outside[1][1:2:end] .<= 1.0))
        @test all(isapprox.(interface[1][1:2:end], 0.5))
        @test interface[3] == [0, 1]
        @test interface[2] == [0, 2]
    finally
        ccall(Libdl.dlsym(lib, :edcut_result_free), Cvoid, (Ptr{Cvoid},), result[])
    end

    @testset "Cut geometry edge cases" begin
        for values in Iterators.product((-1.0, 0.0, 1.0), (-1.0, 0.0, 1.0),
                                        (-1.0, 0.0, 1.0), (-1.0, 0.0, 1.0))
            (any(<(0.0), values) && any(>(0.0), values)) || continue
            status, edge_result, message = _create_quad(coordinates, collect(values))
            @test status == 0
            @test edge_result != C_NULL
            @test isempty(message)
            negative_data = _readRegion(edge_result, 0)
            positive_data = _readRegion(edge_result, 2)
            @test all(negative_data[4] .== 2)
            @test all(positive_data[4] .== 2)
            @test isapprox(_region_area(negative_data) + _region_area(positive_data),
                           1.0; atol=1e-12)
            ccall(Libdl.dlsym(lib, :edcut_result_free), Cvoid, (Ptr{Cvoid},), edge_result)
        end

        for values in ([-1.0, 1.0, -1.0, 1.0], [1.0, -1.0, 1.0, -1.0])
            status, edge_result, _ = _create_quad(coordinates, values)
            @test status == 0
            inside = _readRegion(edge_result, 0)
            outside = _readRegion(edge_result, 2)
            @test isapprox(_region_area(inside), values[1] < 0 ? 0.25 : 0.75; atol=1e-12)
            @test isapprox(_region_area(outside), values[1] < 0 ? 0.75 : 0.25; atol=1e-12)
            @test isapprox(_region_area(inside) + _region_area(outside), 1.0; atol=1e-12)
            ccall(Libdl.dlsym(lib, :edcut_result_free), Cvoid, (Ptr{Cvoid},), edge_result)
        end

        for samples in ([-1.0, 1.0, 1.0, -1.0], [-1.0, 1.0, -1.0, 1.0]),
            scale in (1e-100, 1e100)
            status, reference_result, _ = _create_quad(coordinates, samples)
            @test status == 0
            status, scaled_result, _ = _create_quad(coordinates, scale .* samples)
            @test status == 0
            try
                for region in (0, 2)
                    @test isapprox(_region_area(_readRegion(reference_result, region)),
                                   _region_area(_readRegion(scaled_result, region)); atol=1e-12)
                end
            finally
                ccall(Libdl.dlsym(lib, :edcut_result_free), Cvoid, (Ptr{Cvoid},), reference_result)
                ccall(Libdl.dlsym(lib, :edcut_result_free), Cvoid, (Ptr{Cvoid},), scaled_result)
            end
        end

        status, empty_result, message = _create_quad(coordinates, zeros(4))
        @test status != 0
        @test empty_result == C_NULL
        @test !isempty(message)
    end

    @testset "ABI boundary errors" begin
        status, null_result, message = _create_quad(coordinates, [NaN, 1.0, 1.0, -1.0])
        @test status != 0
        @test null_result == C_NULL
        @test occursin("finite", message)

        zero_status, zero_result, zero_message = _create_quad(coordinates, zeros(4))
        @test zero_status != 0
        @test zero_result == C_NULL
        @test occursin("all level-set values are zero", zero_message)

        count_1, count_2, count_3, count_4 = (Ref{UInt64}() for _ in 1:4)
        err .= 0
        status = GC.@preserve err ccall(
            Libdl.dlsym(lib, :edcut_region_sizes), Cint,
            (Ptr{Cvoid}, Cint, Ref{UInt64}, Ref{UInt64}, Ref{UInt64}, Ref{UInt64}, Ptr{UInt8}, UInt64),
            C_NULL, 0, count_1, count_2, count_3, count_4, pointer(err), _errcap)
        @test status != 0
        @test occursin("invalid result", _native_error(err))

        status, valid_result, _ = _create_quad(coordinates, values)
        @test status == 0
        try
            err .= 0
            status = GC.@preserve err ccall(
                Libdl.dlsym(lib, :edcut_region_sizes), Cint,
                (Ptr{Cvoid}, Cint, Ref{UInt64}, Ref{UInt64}, Ref{UInt64}, Ref{UInt64}, Ptr{UInt8}, UInt64),
                valid_result, 3, count_1, count_2, count_3, count_4, pointer(err), _errcap)
            @test status != 0
            @test occursin("invalid result, region", _native_error(err))

            points_count, offsets_count, connectivity_count, types_count =
                [Ref{UInt64}() for _ in 1:4]
            err .= 0
            status = GC.@preserve err ccall(
                Libdl.dlsym(lib, :edcut_region_sizes), Cint,
                (Ptr{Cvoid}, Cint, Ref{UInt64}, Ref{UInt64}, Ref{UInt64}, Ref{UInt64}, Ptr{UInt8}, UInt64),
                valid_result, 0, points_count, offsets_count, connectivity_count, types_count,
                pointer(err), _errcap)
            @test status == 0
            points, offsets = fill(-7.0, points_count[]), fill(Int32(-7), offsets_count[])
            connectivity, types = fill(Int32(-7), connectivity_count[]), fill(Int32(-7), types_count[])
            before = (copy(points), copy(offsets), copy(connectivity), copy(types))
            err .= 0
            status = GC.@preserve points offsets connectivity types err ccall(
                Libdl.dlsym(lib, :edcut_region_copy), Cint,
                (Ptr{Cvoid}, Cint, Ptr{Cdouble}, UInt64, Ptr{Int32}, UInt64,
                 Ptr{Int32}, UInt64, Ptr{Int32}, UInt64, Ptr{UInt8}, UInt64),
                valid_result, 0, pointer(points), length(points) - 1, pointer(offsets), length(offsets),
                pointer(connectivity), length(connectivity), pointer(types), length(types), pointer(err), _errcap)
            @test status != 0
            @test occursin("capacity", _native_error(err))
            @test (points, offsets, connectivity, types) == before
        finally
            ccall(Libdl.dlsym(lib, :edcut_result_free), Cvoid, (Ptr{Cvoid},), valid_result)
        end
    end

    invalid_result = Ref{Ptr{Cvoid}}(C_NULL)
    err .= 0
    status = GC.@preserve coordinates values err ccall(
        Libdl.dlsym(lib, :edcut_quad_create), Cint,
        (Ptr{Cdouble}, UInt64, Ptr{Cdouble}, UInt64, Ref{Ptr{Cvoid}}, Ptr{UInt8}, UInt64),
        pointer(coordinates), length(coordinates) - 1, pointer(values), length(values),
        invalid_result, pointer(err), _errcap)
    @test status != 0
    @test invalid_result[] == C_NULL
    @test occursin("requires 8 coordinate scalars", _native_error(err))
finally
    Libdl.dlclose(lib)
end
