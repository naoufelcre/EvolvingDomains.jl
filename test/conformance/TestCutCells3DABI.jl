using Test
using Libdl
using LinearAlgebra

const _native_path = Main.EvolvingDomains.CutCellsBinding._libraryPath()

const _lib3d = Libdl.dlopen(_native_path)
const _errcap = 512

# VTK hexahedron order: (000,100,110,010,001,101,111,011).
const _unit_hex = Float64[0, 0, 0, 1, 0, 0, 1, 1, 0, 0, 1, 0, 0, 0, 1, 1, 0, 1, 1, 1, 1, 0, 1, 1]
const _mag = Float64[0.7, 1.3, 2.1, 0.9, 1.7, 0.35, 2.5, 1.1]

function _native_error(buffer)
    first_zero = findfirst(iszero, buffer)
    bytes = first_zero === nothing ? buffer : @view(buffer[1:first_zero-1])
    return String(copy(bytes))
end

function _create_hex(coordinates, values)
    result = Ref{Ptr{Cvoid}}(C_NULL)
    err = zeros(UInt8, _errcap)
    status = GC.@preserve coordinates values err ccall(
        Libdl.dlsym(_lib3d, :edcut_hex_create), Cint,
        (Ptr{Cdouble}, UInt64, Ptr{Cdouble}, UInt64, Ref{Ptr{Cvoid}}, Ptr{UInt8}, UInt64),
        pointer(coordinates), length(coordinates), pointer(values), length(values),
        result, pointer(err), _errcap)
    return status, result[], _native_error(err)
end

function _region_sizes(result, region)
    counts = [Ref{UInt64}() for _ in 1:4]
    err = zeros(UInt8, _errcap)
    status = GC.@preserve counts err ccall(
        Libdl.dlsym(_lib3d, :edcut_region_sizes), Cint,
        (Ptr{Cvoid}, Cint, Ref{UInt64}, Ref{UInt64}, Ref{UInt64}, Ref{UInt64}, Ptr{UInt8}, UInt64),
        result, region, counts[1], counts[2], counts[3], counts[4], pointer(err), _errcap)
    return status, UInt64[c[] for c in counts], _native_error(err)
end

function _readRegion(result, region)
    status, counts, message = _region_sizes(result, region)
    @test status == 0
    @test isempty(message)
    points = zeros(Float64, counts[1])
    offsets = zeros(Int32, counts[2])
    connectivity = zeros(Int32, counts[3])
    types = zeros(Int32, counts[4])
    err = zeros(UInt8, _errcap)
    status = GC.@preserve points offsets connectivity types err ccall(
        Libdl.dlsym(_lib3d, :edcut_region_copy), Cint,
        (Ptr{Cvoid}, Cint, Ptr{Cdouble}, UInt64, Ptr{Int32}, UInt64,
         Ptr{Int32}, UInt64, Ptr{Int32}, UInt64, Ptr{UInt8}, UInt64),
        result, region, pointer(points), length(points), pointer(offsets), length(offsets),
        pointer(connectivity), length(connectivity), pointer(types), length(types),
        pointer(err), _errcap)
    @test status == 0
    @test isempty(_native_error(err))
    return points, offsets, connectivity, types
end

_free(result) = ccall(Libdl.dlsym(_lib3d, :edcut_result_free), Cvoid, (Ptr{Cvoid},), result)

function _tet_volume(points, ids)
    i = Int.(ids) .- 1
    p = [(points[3k + 1], points[3k + 2], points[3k + 3]) for k in i]
    u = p[2] .- p[1]
    v = p[3] .- p[1]
    w = p[4] .- p[1]
    return abs(u[1] * (v[2] * w[3] - v[3] * w[2]) -
               u[2] * (v[1] * w[3] - v[3] * w[1]) +
               u[3] * (v[1] * w[2] - v[2] * w[1])) / 6
end

_cells(data) = ((Int(data[2][cell]) + 1):Int(data[2][cell + 1]) for cell in eachindex(data[4]))

function _volume(data)
    points, _, connectivity, types = data
    volume = 0.0
    for (cell, ids) in enumerate(_cells(data))
        @test types[cell] == 3
        @test length(ids) == 4
        vertex_ids = connectivity[ids] .+ 1
        volume += _tet_volume(points, vertex_ids)
    end
    return volume
end

function _triangle(points, connectivity, ids)
    return [[points[3 * i + 1 + d] for d in 0:2] for i in Int.(connectivity[ids])]
end

# Interface area and unit normals of its triangles.
function _interface_pieces(data)
    points, _, connectivity, types = data
    area = 0.0
    normals = Vector{Vector{Float64}}()
    centroids = Vector{Vector{Float64}}()
    facet_areas = Float64[]
    for (cell, ids) in enumerate(_cells(data))
        @test types[cell] == 2
        @test length(ids) == 3
        vertex = _triangle(points, connectivity, ids)
        u = vertex[2] .- vertex[1]
        v = vertex[3] .- vertex[1]
        n = [u[2] * v[3] - u[3] * v[2], u[3] * v[1] - u[1] * v[3], u[1] * v[2] - u[2] * v[1]]
        push!(facet_areas, norm(n) / 2)
        area += norm(n) / 2
        push!(normals, n ./ norm(n))
        push!(centroids, (vertex[1] .+ vertex[2] .+ vertex[3]) ./ 3)
    end
    return area, normals, centroids, facet_areas
end

_vertex(points, id) = [points[3id - 2], points[3id - 1], points[3id]]

# Facets that touch a sample lying exactly on the plane are built from a 1e-12 nudge,
# so their normals are only accurate to ~1e-4; ignore those when checking direction.
_main_normals(result) = [n for (n, a) in zip(result.normals, result.facet_areas) if a > 1e-9]

# Is p inside any of the region's tetrahedra (barycentric test)?
function _inside_region(data, p)
    points, _, connectivity, _ = data
    for ids in _cells(data)
        a, b, c, d = (_vertex(points, Int(connectivity[id]) + 1) for id in ids)
        m = hcat(b .- a, c .- a, d .- a)
        abs(det(m)) < 1e-14 && continue
        w = m \ (p .- a)
        if all(>(-1e-7), w) && sum(w) < 1 + 1e-7
            return true
        end
    end
    return false
end

# Every interface triangle must face phi>0: stepping 1e-5 along its normal leaves the
# phi<0 tetrahedra and the same step backwards leaves the phi>0 tetrahedra.
function _misoriented(result)
    bad = 0
    for (n, c) in zip(result.normals, result.centroids)
        if !_inside_region(result.positive, c .+ 1e-5 .* n)
            bad += 1
        elseif !_inside_region(result.negative, c .- 1e-5 .* n)
            bad += 1
        end
    end
    return bad
end

# Volume of {x+y+z <= s} and area of its cut for 0 <= s <= 3 in the unit cube.
function _simplex_volume(s)
    s <= 0 && return 0.0
    s <= 1 && return s^3 / 6
    s <= 2 && return (s^3 - 3 * (s - 1)^3) / 6
    s <= 3 && return 1 - (3 - s)^3 / 6
    return 1.0
end

_simplex_area(s) = (s <= 0 || s >= 3) ? 0.0 : sqrt(3.0) * (_simplex_volume(s + 1e-6) - _simplex_volume(s - 1e-6)) / 2e-6

function _check_cut(values; atol=1e-12)
    status, result, message = _create_hex(_unit_hex, values)
    @test status == 0
    @test result != C_NULL
    @test isempty(message)
    try
        negative = _readRegion(result, 0)
        interface = _readRegion(result, 1)
        positive = _readRegion(result, 2)
        negative_volume = _volume(negative)
        positive_volume = _volume(positive)
        area, normals, centroids, facet_areas = _interface_pieces(interface)
        for data in (negative, interface, positive)
            @test data[2][1] == 0
            @test issorted(data[2])
            @test data[2][end] == length(data[3])
            @test all(0 .<= data[3] .< length(data[1]) ÷ 3)
        end
        return (negative_volume=negative_volume, positive_volume=positive_volume,
                area=area, normals=normals, centroids=centroids, facet_areas=facet_areas,
                negative=negative, positive=positive)
    finally
        _free(result)
    end
end

try
    @test unsafe_string(ccall(Libdl.dlsym(_lib3d, :edcut_upstream_version), Cstring, ())) == "0.4.0"

    @testset "analytic planes" begin
        # phi = z - 0.5: half cube below, interface is the unit square.
        result = _check_cut([-0.5, -0.5, -0.5, -0.5, 0.5, 0.5, 0.5, 0.5])
        @test isapprox(result.negative_volume, 0.5; atol=1e-12)
        @test isapprox(result.positive_volume, 0.5; atol=1e-12)
        @test isapprox(result.area, 1.0; atol=1e-12)
        @test all(n[3] > 0.999 for n in result.normals)
        @test _misoriented(result) == 0

        # phi = x - 0.5 with scaled samples: same cut, same orientation (+x).
        result = _check_cut([-1e3, 1e3, 1e3, -1e3, -1e3, 1e3, 1e3, -1e3])
        @test isapprox(result.negative_volume, 0.5; atol=1e-12)
        @test isapprox(result.positive_volume, 0.5; atol=1e-12)
        @test isapprox(result.area, 1.0; atol=1e-12)
        @test all(n[1] > 0.999 for n in result.normals)

        # phi = x + y + z - s, normal (1,1,1)/sqrt(3); phi<0 is the far side.
        for s in (0.5, 1.0, 1.2, 1.5, 1.8, 2.5)
            values = [s - _unit_hex[3i - 2] - _unit_hex[3i - 1] - _unit_hex[3i] for i in 1:8]
            result = _check_cut(values)
            @test isapprox(result.negative_volume, 1 - _simplex_volume(s); atol=1e-10)
            @test isapprox(result.positive_volume, _simplex_volume(s); atol=1e-10)
            @test isapprox(result.area, _simplex_area(s); atol=1e-9)
            # phi>0 is the near side, so the outward interface normal is -(1,1,1)/sqrt(3).
            # s in {1,2} puts samples exactly on the plane, hence the loose tolerance.
            @test all(isapprox(n[1], -1 / sqrt(3); atol=1e-3) for n in _main_normals(result))
            @test _misoriented(result) == 0
        end

        # phi = x - y: the plane contains four vertices, so zeros exercise the same
        # path as the edge/vertex cases below.
        result = _check_cut([0.0, 1.0, 0.0, -1.0, 0.0, 1.0, 0.0, -1.0])
        @test isapprox(result.negative_volume, 0.5; atol=1e-9)
        @test isapprox(result.positive_volume, 0.5; atol=1e-9)
        @test isapprox(result.area, sqrt(2.0); atol=1e-9)
        @test all(isapprox(n[2], -1 / sqrt(2); atol=1e-3) for n in _main_normals(result))
    end

    @testset "pure full cells" begin
        for (sign, region) in ((-1.0, 0), (1.0, 2))
            status, result, message = _create_hex(_unit_hex, fill(sign, 8))
            @test status == 0
            @test isempty(message)
            try
                full = _readRegion(result, region)
                @test isapprox(_volume(full), 1.0; atol=1e-12)
                @test length(full[4]) == 6
                @test _volume(_readRegion(result, region == 0 ? 2 : 0)) == 0.0
                interface = _readRegion(result, 1)
                @test isempty(interface[4])
            finally
                _free(result)
            end
        end

        # A face lying in the phi=0 plane is not cut: the whole cell is positive.
        status, result, message = _create_hex(_unit_hex, [0.0, 0, 0, 0, 1.0, 1, 1, 1])
        @test status == 0
        @test isempty(message)
        try
            @test isapprox(_volume(_readRegion(result, 2)), 1.0; atol=1e-12)
            @test isempty(_readRegion(result, 1)[4])
        finally
            _free(result)
        end
    end

    @testset "zero, edge and vertex on the plane" begin
        # Well conditioned: a zero at one vertex, along one edge, or on the whole
        # diagonal, so both phases still have volume next to every interface facet.
        for values in (
            [0.0, 1.0, 1.0, 1.0, -1.0, -1.0, -1.0, -1.0],  # one vertex on the plane
            [0.0, 0.0, 1.0, 1.0, -1.0, -1.0, -1.0, -1.0],  # one edge on the plane
            [-1.0, 0.0, 1.0, 0.0, -1.0, 0.0, 1.0, 0.0],    # four vertices on the plane
        )
            result = _check_cut(values)
            @test isapprox(result.negative_volume + result.positive_volume, 1.0; atol=1e-9)
            @test all(isfinite, vcat(result.negative[1], result.positive[1], result.normals...))
            @test _misoriented(result) == 0
        end

        # Degenerate: three coplanar samples collapse one phase to a ~1e-12 sliver, so
        # only partition, finiteness and positive areas are meaningful.
        for values in (
            [0.0, 0.0, 0.0, -1.0, 1.0, 1.0, 1.0, 0.0],     # three coplanar vertices
            [0.0, 0.0, 0.0, 0.0, 1.0, 1.0, 1.0, 1.0],      # whole face on the plane
        )
            result = _check_cut(values)
            @test isapprox(result.negative_volume + result.positive_volume, 1.0; atol=1e-9)
            @test all(isfinite, vcat(result.negative[1], result.positive[1], result.normals...))
            @test all(>=(0), result.facet_areas)
            @test isapprox(result.area, sum(result.facet_areas); atol=1e-12)
        end
    end

    @testset "all mixed sign masks" begin
        mixed = [(mask, [isodd(mask >> bit) ? -_mag[bit + 1] : _mag[bit + 1] for bit in 0:7])
                 for mask in 1:254 if mask != 255 && mask != 0]
        # keep only masks with both signs present
        mixed = [(mask, values) for (mask, values) in mixed
                 if any(<(0), values) && any(>(0), values)]
        @test length(mixed) == 254
        for (mask, values) in mixed
            result = _check_cut(values; atol=1e-10)
            total = result.negative_volume + result.positive_volume
            @test result.negative_volume >= -1e-12
            @test result.positive_volume >= -1e-12
            # One property, one assertion: the old failure collector re-checked this and
            # only fed an @info log.
            @test isapprox(total, 1.0; atol=1e-10)
            @test result.area > 0
            @test _misoriented(result) == 0
            # Complement mask exchanges the phases and their volumes.
            complement = [isodd(mask >> bit) ? _mag[bit + 1] : -_mag[bit + 1] for bit in 0:7]
            other = _check_cut(complement; atol=1e-10)
            @test isapprox(other.negative_volume, result.positive_volume; atol=1e-10)
            @test isapprox(other.positive_volume, result.negative_volume; atol=1e-10)
            @test isapprox(other.area, result.area; atol=1e-10)
        end
    end

    @testset "ABI boundaries" begin
        for (coordinates, values, needle) in (
            (_unit_hex[1:end-3], fill(-1.0, 8), "24 coordinate scalars"),
            (_unit_hex, fill(-1.0, 7), "24 coordinate scalars"),
            (_unit_hex, fill(NaN, 8), "finite"),
            (vcat(_unit_hex[1:end-1], [Inf]), fill(-1.0, 8), "finite"),
            (_unit_hex, zeros(8), "all level-set values are zero"),
        )
            result = Ref{Ptr{Cvoid}}(Ptr{Cvoid}(1))
            err = zeros(UInt8, _errcap)
            status = GC.@preserve coordinates values err ccall(
                Libdl.dlsym(_lib3d, :edcut_hex_create), Cint,
                (Ptr{Cdouble}, UInt64, Ptr{Cdouble}, UInt64, Ref{Ptr{Cvoid}}, Ptr{UInt8}, UInt64),
                pointer(coordinates), length(coordinates), pointer(values), length(values),
                result, pointer(err), _errcap)
            @test status != 0
            @test result[] == C_NULL
            @test occursin(needle, _native_error(err))
        end

        err = zeros(UInt8, _errcap)
        status = GC.@preserve err ccall(
            Libdl.dlsym(_lib3d, :edcut_hex_create), Cint,
            (Ptr{Cdouble}, UInt64, Ptr{Cdouble}, UInt64, Ref{Ptr{Cvoid}}, Ptr{UInt8}, UInt64),
            C_NULL, 24, C_NULL, 8, Ref{Ptr{Cvoid}}(C_NULL), pointer(err), _errcap)
        @test status != 0
        @test occursin("24 coordinate scalars", _native_error(err))

        counts = [Ref{UInt64}() for _ in 1:4]
        err = zeros(UInt8, _errcap)
        status = GC.@preserve counts err ccall(
            Libdl.dlsym(_lib3d, :edcut_region_sizes), Cint,
            (Ptr{Cvoid}, Cint, Ref{UInt64}, Ref{UInt64}, Ref{UInt64}, Ref{UInt64}, Ptr{UInt8}, UInt64),
            C_NULL, 0, counts[1], counts[2], counts[3], counts[4], pointer(err), _errcap)
        @test status != 0
        @test occursin("invalid result", _native_error(err))

        status, result, _ = _create_hex(_unit_hex, [-0.5, -0.5, -0.5, -0.5, 0.5, 0.5, 0.5, 0.5])
        @test status == 0
        try
            err .= 0
            status = GC.@preserve counts err ccall(
                Libdl.dlsym(_lib3d, :edcut_region_sizes), Cint,
                (Ptr{Cvoid}, Cint, Ref{UInt64}, Ref{UInt64}, Ref{UInt64}, Ref{UInt64}, Ptr{UInt8}, UInt64),
                result, 3, counts[1], counts[2], counts[3], counts[4], pointer(err), _errcap)
            @test status != 0
            @test occursin("invalid result, region", _native_error(err))

            err .= 0
            status = GC.@preserve counts err ccall(
                Libdl.dlsym(_lib3d, :edcut_region_sizes), Cint,
                (Ptr{Cvoid}, Cint, Ptr{UInt64}, Ref{UInt64}, Ref{UInt64}, Ref{UInt64}, Ptr{UInt8}, UInt64),
                result, 2, C_NULL, counts[2], counts[3], counts[4], pointer(err), _errcap)
            @test status != 0
            @test occursin("output counts must not be null", _native_error(err))

            _, full_counts, _ = _region_sizes(result, 2)
            points = fill(-7.0, full_counts[1])
            offsets = fill(Int32(-7), full_counts[2])
            connectivity = fill(Int32(-7), full_counts[3])
            types = fill(Int32(-7), full_counts[4])
            before = (copy(points), copy(offsets), copy(connectivity), copy(types))
            err .= 0
            status = GC.@preserve points offsets connectivity types err ccall(
                Libdl.dlsym(_lib3d, :edcut_region_copy), Cint,
                (Ptr{Cvoid}, Cint, Ptr{Cdouble}, UInt64, Ptr{Int32}, UInt64,
                 Ptr{Int32}, UInt64, Ptr{Int32}, UInt64, Ptr{UInt8}, UInt64),
                result, 2, pointer(points), length(points) - 1, pointer(offsets), length(offsets),
                pointer(connectivity), length(connectivity), pointer(types), length(types),
                pointer(err), _errcap)
            @test status != 0
            @test occursin("capacity", _native_error(err))
            @test (points, offsets, connectivity, types) == before
        finally
            _free(result)
        end
    end

    @testset "ownership" begin
        status, first_result, _ = _create_hex(_unit_hex, [-0.5, -0.5, -0.5, -0.5, 0.5, 0.5, 0.5, 0.5])
        @test status == 0
        status, second_result, _ = _create_hex(_unit_hex, fill(-1.0, 8))
        @test status == 0
        @test first_result != second_result
        _free(first_result)
        # The surviving handle is still usable after the sibling is freed.
        @test isapprox(_volume(_readRegion(second_result, 0)), 1.0; atol=1e-12)
        _free(second_result)
        # Empty regions can be queried and copied with null storage.
        status, full_result, _ = _create_hex(_unit_hex, fill(1.0, 8))
        @test status == 0
        try
            # An empty region still reports one offset (cells + 1), as in the 2-D ABI.
            _, empty_counts, _ = _region_sizes(full_result, 0)
            @test empty_counts == [0, 1, 0, 0]
            offsets = zeros(Int32, 1)
            err = zeros(UInt8, _errcap)
            status = GC.@preserve offsets err ccall(
                Libdl.dlsym(_lib3d, :edcut_region_copy), Cint,
                (Ptr{Cvoid}, Cint, Ptr{Cdouble}, UInt64, Ptr{Int32}, UInt64,
                 Ptr{Int32}, UInt64, Ptr{Int32}, UInt64, Ptr{UInt8}, UInt64),
                full_result, 0, C_NULL, 0, pointer(offsets), 1, C_NULL, 0, C_NULL, 0, pointer(err), _errcap)
            @test status == 0
            @test isempty(_native_error(err))
            @test offsets == [0]
        finally
            _free(full_result)
        end
    end
finally
    Libdl.dlclose(_lib3d)
end
