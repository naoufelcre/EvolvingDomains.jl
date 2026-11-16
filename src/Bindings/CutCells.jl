# Private native CutCells backend. Users call `cut(field, level)`; this namespace
# holds the C ABI, library loading, and conversion helpers.
module CutCellsBinding

using Artifacts
using Libdl
using ..EvolvingDomains: Cut, Field, _cutFromNative, _levelgrid, _levelvalues,
    _UNIFORM_RTOL, _coordinateTolerance

const _library_handle = Ref{Ptr{Cvoid}}(C_NULL)
const _library_lock = ReentrantLock()
const _native_error_capacity = 512
const _native_regions = (0, 1, 2)
const _native_interval = Int32(1)
const _native_triangle = Int32(2)
const _native_tetrahedron = Int32(3)
const _package_uuid = Base.UUID("f9bd8acb-b7a8-42b3-8348-bb9bbfd76943")

function _checkGlibcVersion()
    Sys.islinux() || return nothing
    version = try
        Libdl.dlopen("libc.so.6") do handle
            VersionNumber(unsafe_string(ccall(Libdl.dlsym(handle, :gnu_get_libc_version), Cstring, ())))
        end
    catch err
        throw(ArgumentError("CutCells Linux artifacts require glibc 2.28 or newer: $(sprint(showerror, err))"))
    end
    version >= v"2.28" || throw(ArgumentError(
        "CutCells Linux artifacts require glibc 2.28 or newer; found $version"))
    return nothing
end

function _libraryPath()
    _checkGlibcVersion()
    override = get(ENV, "EVOLVINGDOMAINS_CUTCELLS_LIBRARY", "")
    if !isempty(override)
        path = abspath(override)
        isfile(path) || throw(ArgumentError("CutCells shared library does not exist: $path"))
        return path
    end
    artifact_toml = normpath(joinpath(@__DIR__, "..", "..", "Artifacts.toml"))
    isfile(artifact_toml) || throw(ArgumentError(_noArtifactMessage(
        "CutCells platform artifacts have not been published in Artifacts.toml")))
    hash = Artifacts.artifact_hash("cutcells_bindings", artifact_toml; pkg_uuid=_package_uuid)
    hash === nothing && throw(ArgumentError(_noArtifactMessage(
        "no CutCells artifact is registered for $(Base.BinaryPlatforms.triplet(Base.BinaryPlatforms.HostPlatform()))")))
    filename = (Sys.iswindows() ? "" : "lib") * "evolvingdomains_cut_cells." * Libdl.dlext
    path = joinpath(Artifacts.artifact_path(hash), "lib", filename)
    isfile(path) || throw(ArgumentError("CutCells artifact is missing its native library: $path"))
    return path
end

# A missing artifact is a dead end unless the two local ways forward are named.
function _noArtifactMessage(problem::AbstractString)
    root = normpath(joinpath(@__DIR__, "..", ".."))
    return string(problem, ". Either point EVOLVINGDOMAINS_CUTCELLS_LIBRARY at an existing ",
        "library, or build one locally with \`", joinpath(root, "scripts", "build_cxx_bindings.sh"),
        "\` and set that variable to its lib/ output.")
end

function _library()
    handle = _library_handle[]
    handle != C_NULL && return handle
    lock(_library_lock) do
        handle = _library_handle[]
        if handle == C_NULL
            path = _libraryPath()
            handle = try
                Libdl.dlopen(path)
            catch err
                throw(ArgumentError("could not load CutCells native library at $path: $(sprint(showerror, err))"))
            end
            _library_handle[] = handle
        end
    end
    return _library_handle[]
end

_nativeSymbol(name::Symbol) = Libdl.dlsym(_library(), name)

function _errorMessage(buffer)
    first_zero = findfirst(iszero, buffer)
    bytes = first_zero === nothing ? buffer : @view(buffer[1:first_zero-1])
    return String(copy(bytes))
end

function _checkNative(status, buffer, operation)
    status == 0 && return nothing
    message = _errorMessage(buffer)
    isempty(message) && (message = "native status $status")
    throw(ErrorException("CutCells $operation failed: $message"))
end

function _readRegion(result, region::Int, ndim::Int=2)
    points_count, offsets_count = Ref{UInt64}(), Ref{UInt64}()
    connectivity_count, types_count = Ref{UInt64}(), Ref{UInt64}()
    err = zeros(UInt8, _native_error_capacity)
    status = GC.@preserve err ccall(
        _nativeSymbol(:edcut_region_sizes), Cint,
        (Ptr{Cvoid}, Cint, Ref{UInt64}, Ref{UInt64}, Ref{UInt64}, Ref{UInt64}, Ptr{UInt8}, UInt64),
        result, region, points_count, offsets_count, connectivity_count, types_count,
        pointer(err), _native_error_capacity)
    _checkNative(status, err, "region size query")

    points = zeros(Float64, points_count[])
    offsets = zeros(Int32, offsets_count[])
    connectivity = zeros(Int32, connectivity_count[])
    types = zeros(Int32, types_count[])
    fill!(err, 0)
    status = GC.@preserve points offsets connectivity types err ccall(
        _nativeSymbol(:edcut_region_copy), Cint,
        (Ptr{Cvoid}, Cint, Ptr{Cdouble}, UInt64, Ptr{Int32}, UInt64,
         Ptr{Int32}, UInt64, Ptr{Int32}, UInt64, Ptr{UInt8}, UInt64),
        result, region, pointer(points), length(points), pointer(offsets), length(offsets),
        pointer(connectivity), length(connectivity), pointer(types), length(types),
        pointer(err), _native_error_capacity)
    _checkNative(status, err, "region copy")
    length(points) % ndim == 0 && all(isfinite, points) ||
        throw(ErrorException("CutCells returned invalid coordinates"))
    !isempty(offsets) && first(offsets) == 0 && last(offsets) == length(connectivity) &&
        length(types) + 1 == length(offsets) && issorted(offsets) ||
        throw(ErrorException("CutCells returned invalid CSR topology"))
    all(i -> 0 <= i < length(points) ÷ ndim, connectivity) ||
        throw(ErrorException("CutCells returned out-of-range point connectivity"))
    return points, offsets, connectivity, types
end

function _cutQuad(values::NTuple{4,Float64})
    coordinates = Float64[0, 0, 1, 0, 1, 1, 0, 1]
    samples = collect(values)
    result = Ref{Ptr{Cvoid}}(C_NULL)
    err = zeros(UInt8, _native_error_capacity)
    status = GC.@preserve coordinates samples err ccall(
        _nativeSymbol(:edcut_quad_create), Cint,
        (Ptr{Cdouble}, UInt64, Ptr{Cdouble}, UInt64, Ref{Ptr{Cvoid}}, Ptr{UInt8}, UInt64),
        coordinates, length(coordinates), samples, length(samples), result,
        pointer(err), _native_error_capacity)
    _checkNative(status, err, "quadrilateral cut")
    result[] == C_NULL && throw(ErrorException("CutCells returned a null result"))
    try
        return ntuple(i -> _readRegion(result[], _native_regions[i]), 3)
    finally
        ccall(_nativeSymbol(:edcut_result_free), Cvoid, (Ptr{Cvoid},), result[])
    end
end

function _cutHex(values::NTuple{8,Float64})
    coordinates = Float64[0, 0, 0, 1, 0, 0, 1, 1, 0, 0, 1, 0,
                          0, 0, 1, 1, 0, 1, 1, 1, 1, 0, 1, 1]
    samples = collect(values)
    result = Ref{Ptr{Cvoid}}(C_NULL)
    err = zeros(UInt8, _native_error_capacity)
    status = GC.@preserve coordinates samples err ccall(
        _nativeSymbol(:edcut_hex_create), Cint,
        (Ptr{Cdouble}, UInt64, Ptr{Cdouble}, UInt64, Ref{Ptr{Cvoid}}, Ptr{UInt8}, UInt64),
        coordinates, length(coordinates), samples, length(samples), result,
        pointer(err), _native_error_capacity)
    _checkNative(status, err, "hexahedral cut")
    result[] == C_NULL && throw(ErrorException("CutCells returned a null result"))
    try
        return ntuple(i -> _readRegion(result[], _native_regions[i], 3), 3)
    finally
        ccall(_nativeSymbol(:edcut_result_free), Cvoid, (Ptr{Cvoid},), result[])
    end
end

# Same tolerance and the same atol rule as `Geometric._gridMetrics`: the producer must
# accept exactly the axes that every numerical operator accepts.
function _uniformSpacing(axis)
    n = length(axis)
    x0, x1 = Float64(first(axis)), Float64(last(axis))
    spacing = (x1 - x0) / (n - 1)
    tolerance = _coordinateTolerance(x0, x1, spacing)
    all(i -> abs(Float64(axis[i]) - (x0 + (i - 1) * spacing)) <= tolerance, eachindex(axis)) ||
        throw(ArgumentError("CutCells adapter requires uniformly spaced grid axes"))
    return spacing
end

@inline function _facetState(a, b)
    a < 0 && b < 0 && return Int8(-1)
    a > 0 && b > 0 && return Int8(1)
    return Int8(0)
end

@inline function _cellState(phi)
    all(<(0.0), phi) && return Int8(-1)
    all(>(0.0), phi) && return Int8(1)
    return Int8(0)
end

function _appendPoints!(physical, reference, points, x0, y0, hx, hy)
    base = length(physical)
    for k in 1:2:length(points)
        ξ, η = points[k], points[k + 1]
        push!(reference, (ξ, η))
        push!(physical, (x0 + hx * ξ, y0 + hy * η))
    end
    return base
end

function _appendVolumeRegion!(physical, reference, states, parents, ids, pointers,
                                data, state, parent, x0, y0, hx, hy)
    points, offsets, connectivity, types = data
    base = _appendPoints!(physical, reference, points, x0, y0, hx, hy)
    ncells = max(length(offsets) - 1, 0)
    length(types) == ncells || throw(ErrorException("CutCells returned mismatched cell offsets and types"))
    for cell in 1:ncells
        types[cell] == _native_triangle || throw(ErrorException("CutCells returned a non-triangle volume subcell"))
        first_index, last_index = Int(offsets[cell]) + 1, Int(offsets[cell + 1])
        last_index - first_index + 1 == 3 || throw(ErrorException("CutCells returned a non-triangular connectivity row"))
        vertices = Int[connectivity[index] for index in first_index:last_index]
        all(id -> 0 <= id < length(points) ÷ 2, vertices) ||
            throw(ErrorException("CutCells returned out-of-range point connectivity"))
        x1, y1 = points[2vertices[1] + 1], points[2vertices[1] + 2]
        x2, y2 = points[2vertices[2] + 1], points[2vertices[2] + 2]
        x3, y3 = points[2vertices[3] + 1], points[2vertices[3] + 2]
        signed_area = (x2 - x1) * (y3 - y1) - (y2 - y1) * (x3 - x1)
        iszero(signed_area) && continue
        if signed_area < 0
            vertices[2], vertices[3] = vertices[3], vertices[2]
        end
        append!(ids, Int32[base + id + 1 for id in vertices])
        push!(pointers, Int32(length(ids) + 1))
        push!(states, state)
        push!(parents, Int32(parent))
    end
    return nothing
end

function _appendBoundaryRegion!(physical, reference, parents, ids, pointers, normals,
                                  data, parent, values, x0, y0, hx, hy)
    points, offsets, connectivity, types = data
    base = _appendPoints!(physical, reference, points, x0, y0, hx, hy)
    ncells = max(length(offsets) - 1, 0)
    length(types) == ncells || throw(ErrorException("CutCells returned mismatched interface offsets and types"))
    for cell in 1:ncells
        types[cell] == _native_interval || throw(ErrorException("CutCells returned a non-segment interface cell"))
        first_index, last_index = Int(offsets[cell]) + 1, Int(offsets[cell + 1])
        last_index - first_index + 1 == 2 || throw(ErrorException("CutCells returned a non-segment connectivity row"))
        pids = Int[]
        for index in first_index:last_index
            0 <= connectivity[index] < length(points) ÷ 2 ||
                throw(ErrorException("CutCells returned out-of-range interface connectivity"))
            id = Int(connectivity[index])
            push!(pids, id)
        end
        ξ1, η1 = points[2pids[1] + 1], points[2pids[1] + 2]
        ξ2, η2 = points[2pids[2] + 1], points[2pids[2] + 2]
        # Point-only contacts are not represented as boundary facets.
        ξ1 == ξ2 && η1 == η2 && continue
        append!(ids, Int32[base + id + 1 for id in pids])
        tx, ty = hx * (ξ2 - ξ1), hy * (η2 - η1)
        length_tangent = hypot(tx, ty)
        length_tangent > 0 || throw(ErrorException("CutCells returned a zero-length interface segment"))
        nx, ny = ty / length_tangent, -tx / length_tangent

        ξ, η = (ξ1 + ξ2) / 2, (η1 + η2) / 2
        gx = ((1 - η) * (values[2] - values[1]) + η * (values[3] - values[4])) / hx
        gy = ((1 - ξ) * (values[4] - values[1]) + ξ * (values[3] - values[2])) / hy
        orientation = nx * gx + ny * gy
        if orientation < 0
            nx, ny = -nx, -ny
        elseif orientation == 0
            δ = min(1e-6, max(0.0, min(ξ, 1 - ξ, η, 1 - η)) / 2)
            δ > 0 || throw(ErrorException("cannot orient a degenerate CutCells interface normal"))
            plus = _bilinearValue(values, ξ + δ * nx / hx, η + δ * ny / hy)
            minus = _bilinearValue(values, ξ - δ * nx / hx, η - δ * ny / hy)
            plus == minus && throw(ErrorException("cannot orient a CutCells interface normal from level-set values"))
            if plus < minus
                nx, ny = -nx, -ny
            end
        end
        push!(pointers, Int32(length(ids) + 1))
        push!(parents, Int32(parent))
        push!(normals, (nx, ny))
    end
    return nothing
end

@inline function _bilinearValue(v, ξ, η)
    return (1 - ξ) * (1 - η) * v[1] + ξ * (1 - η) * v[2] + ξ * η * v[3] + (1 - ξ) * η * v[4]
end

function _appendFullCell!(physical, reference, states, parents, ids, pointers, state,
                            parent, x0, y0, hx, hy)
    base = _appendPoints!(physical, reference,
        (0.0, 0.0, 1.0, 0.0, 1.0, 1.0, 0.0, 1.0), x0, y0, hx, hy)
    for vertex in (0, 1, 2, 0, 2, 3)
        push!(ids, Int32(base + vertex + 1))
    end
    append!(pointers, (Int32(length(ids) - 2), Int32(length(ids) + 1)))
    append!(states, (state, state))
    append!(parents, (Int32(parent), Int32(parent)))
    return nothing
end

@inline _facetState(values...) = all(<(0.0), values) ? Int8(-1) :
    all(>(0.0), values) ? Int8(1) : Int8(0)

function _append3DPoints!(physical, reference, points, axes, indices)
    base = length(physical)
    x0, y0, z0 = axes[1][indices[1]], axes[2][indices[2]], axes[3][indices[3]]
    hx, hy, hz = axes[1][indices[1] + 1] - x0,
                 axes[2][indices[2] + 1] - y0,
                 axes[3][indices[3] + 1] - z0
    for k in 1:3:length(points)
        ξ, η, ζ = points[k], points[k + 1], points[k + 2]
        push!(reference, (ξ, η, ζ))
        push!(physical, (x0 + hx * ξ, y0 + hy * η, z0 + hz * ζ))
    end
    return base
end

@inline function _det3(a, b, c)
    return a[1] * (b[2] * c[3] - b[3] * c[2]) -
           a[2] * (b[1] * c[3] - b[3] * c[1]) +
           a[3] * (b[1] * c[2] - b[2] * c[1])
end

function _append3DVolume!(physical, reference, states, parents, ids, pointers,
                             data, state, parent, axes, indices)
    points, offsets, connectivity, types = data
    base = _append3DPoints!(physical, reference, points, axes, indices)
    ncells = length(offsets) - 1
    length(types) == ncells || throw(ErrorException("CutCells returned mismatched tetrahedron offsets and types"))
    for cell in 1:ncells
        types[cell] == _native_tetrahedron || throw(ErrorException("CutCells returned a non-tetrahedron volume subcell"))
        first_index, last_index = Int(offsets[cell]) + 1, Int(offsets[cell + 1])
        last_index - first_index + 1 == 4 || throw(ErrorException("CutCells returned non-tetrahedral connectivity"))
        vertices = Int[connectivity[q] for q in first_index:last_index]
        all(id -> 0 <= id < length(points) ÷ 3, vertices) ||
            throw(ErrorException("CutCells returned out-of-range point connectivity"))
        p = ntuple(d -> points[3vertices[1] + d], 3)
        q = ntuple(d -> points[3vertices[2] + d] - p[d], 3)
        r = ntuple(d -> points[3vertices[3] + d] - p[d], 3)
        s = ntuple(d -> points[3vertices[4] + d] - p[d], 3)
        volume = _det3(q, r, s)
        iszero(volume) && continue
        if volume < 0
            vertices[3], vertices[4] = vertices[4], vertices[3]
        end
        append!(ids, Int32[base + id + 1 for id in vertices])
        push!(pointers, Int32(length(ids) + 1))
        push!(states, state)
        push!(parents, Int32(parent))
    end
    return nothing
end

function _trilinearGradient(v, ξ, η, ζ, axes, indices)
    hx, hy, hz = axes[1][indices[1] + 1] - axes[1][indices[1]],
                 axes[2][indices[2] + 1] - axes[2][indices[2]],
                 axes[3][indices[3] + 1] - axes[3][indices[3]]
    bx = ((1-η)*(1-ζ)*(v[2]-v[1]) + η*(1-ζ)*(v[3]-v[4]) +
          η*ζ*(v[7]-v[8]) + (1-η)*ζ*(v[6]-v[5])) / hx
    by = ((1-ξ)*(1-ζ)*(v[4]-v[1]) + ξ*(1-ζ)*(v[3]-v[2]) +
          ξ*ζ*(v[7]-v[6]) + (1-ξ)*ζ*(v[8]-v[5])) / hy
    bz = ((1-ξ)*(1-η)*(v[5]-v[1]) + ξ*(1-η)*(v[6]-v[2]) +
          ξ*η*(v[7]-v[3]) + (1-ξ)*η*(v[8]-v[4])) / hz
    return (bx, by, bz)
end

function _append3DBoundary!(physical, reference, parents, ids, pointers, normals,
                              data, parent, values, axes, indices; orient_by_levelset::Bool=false)
    points, offsets, connectivity, types = data
    base = _append3DPoints!(physical, reference, points, axes, indices)
    ncells = length(offsets) - 1
    length(types) == ncells || throw(ErrorException("CutCells returned mismatched interface offsets and types"))
    for cell in 1:ncells
        types[cell] == _native_triangle || throw(ErrorException("CutCells returned a non-triangle interface cell"))
        first_index, last_index = Int(offsets[cell]) + 1, Int(offsets[cell + 1])
        last_index - first_index + 1 == 3 || throw(ErrorException("CutCells returned non-triangular interface connectivity"))
        vertices = Int[connectivity[q] for q in first_index:last_index]
        all(id -> 0 <= id < length(points) ÷ 3, vertices) ||
            throw(ErrorException("CutCells returned out-of-range interface connectivity"))
        # Geometry from the emitted physical points, not the reference ones: a cell is the
        # unit cube stretched by the spacings, so a reference-space normal is only correct
        # when those are equal, and a triangle that collapses during the mapping still looks
        # non-degenerate in reference coordinates.
        p1 = physical[base + vertices[1] + 1]
        p2 = physical[base + vertices[2] + 1]
        p3 = physical[base + vertices[3] + 1]
        q = ntuple(d -> p2[d] - p1[d], 3)
        r = ntuple(d -> p3[d] - p1[d], 3)
        cx, cy, cz = q[2]*r[3]-q[3]*r[2], q[3]*r[1]-q[1]*r[3], q[1]*r[2]-q[2]*r[1]
        twice_area = hypot(cx, cy, cz)
        # Tangent to the simplex (a point or a segment of intersection): no facet to emit.
        twice_area > 0 || continue
        if orient_by_levelset
            ξ = sum(points[3id + 1] for id in vertices) / 3
            η = sum(points[3id + 2] for id in vertices) / 3
            ζ = sum(points[3id + 3] for id in vertices) / 3
            gradient = _trilinearGradient(values, ξ, η, ζ, axes, indices)
            hypot(gradient...) > 0 ||
                throw(ErrorException("cannot orient a degenerate CutCells interface normal"))
            if (cx * gradient[1] + cy * gradient[2] + cz * gradient[3]) < 0
                vertices[2], vertices[3] = vertices[3], vertices[2]
                cx, cy, cz = -cx, -cy, -cz
            end
        end
        append!(ids, Int32[base + id + 1 for id in vertices])
        push!(pointers, Int32(length(ids) + 1))
        push!(parents, Int32(parent))
        # Native triangles are already oriented toward phi>0 on the six-tet model; store
        # the plane normal, which need not equal the trilinear gradient.
        push!(normals, (cx / twice_area, cy / twice_area, cz / twice_area))
    end
    return nothing
end

function _append3DFullCell!(physical, reference, states, parents, ids, pointers,
                               state, parent, axes, indices)
    points = (0.,0.,0., 1.,0.,0., 1.,1.,0., 0.,1.,0.,
              0.,0.,1., 1.,0.,1., 1.,1.,1., 0.,1.,1.)
    base = _append3DPoints!(physical, reference, points, axes, indices)
    for tet in ((1,2,3,7),(1,3,4,7),(1,4,8,7),(1,8,5,7),(1,5,6,7),(1,6,2,7))
        append!(ids, Int32[base + i for i in tet])
        push!(pointers, Int32(length(ids) + 1))
        push!(states, state)
        push!(parents, Int32(parent))
    end
    return nothing
end

function _fromField3D(field, level_id)
    d = field.discretization
    _levelgrid(d, level_id)
    _levelvalues(field, level_id)
    g = d.levels[level_id]
    nx, ny, nz = ntuple(a -> length(g.axes[a]) - 1, 3)
    raw = field.data[level_id]
    size(raw) == (nx + 1, ny + 1, nz + 1) || throw(DimensionMismatch(
        "Field has size $(size(raw)); grid requires $((nx + 1, ny + 1, nz + 1))"))
    eltype(raw) <: Real || throw(ArgumentError("CutCells adapter requires real-valued Field data"))
    all(isfinite, raw) || throw(ArgumentError("Field level set must contain only finite values"))
    phi = Float64.(raw)
    all(isfinite, phi) || throw(ArgumentError("Field level set values must be representable as finite Float64 numbers"))
    all(a -> all(isfinite, a) && all(diff(a) .> 0), g.axes) ||
        throw(ArgumentError("CutCells adapter requires finite, increasing grid axes"))
    # Validation only: 3D maps reference to physical per axis and ignores the spacing.
    foreach(_uniformSpacing, g.axes)

    ncells = nx * ny * nz
    bg = Vector{Int8}(undef, ncells)
    facets = Vector{NTuple{6,Int8}}(undef, ncells)
    states, parents, ids, pointers = Int8[], Int32[], Int32[], Int32[1]
    physical, reference = NTuple{3,Float64}[], NTuple{3,Float64}[]
    bparents, bids, bptrs = Int32[], Int32[], Int32[1]
    bphysical, breference, normals = NTuple{3,Float64}[], NTuple{3,Float64}[], NTuple{3,Float64}[]
    corners = ((0,0,0),(1,0,0),(1,1,0),(0,1,0),(0,0,1),(1,0,1),(1,1,1),(0,1,1))
    for k in 1:nz, j in 1:ny, i in 1:nx
        parent = i + (j-1)*nx + (k-1)*nx*ny
        indices = (i,j,k)
        cv = ntuple(a -> phi[i+corners[a][1],j+corners[a][2],k+corners[a][3]], 8)
        has_negative, has_positive = any(<(0.0), cv), any(>(0.0), cv)
        has_negative || has_positive || throw(ArgumentError("level set is identically zero on background cell $parent; a whole-cell interface is unsupported"))
        state = _cellState(cv)
        bg[parent] = state
        faces = ((1,4,8,5),(2,3,7,6),(1,2,6,5),(4,3,7,8),(1,2,3,4),(5,6,7,8))
        facets[parent] = ntuple(q -> _facetState(cv[faces[q][1]],cv[faces[q][2]],cv[faces[q][3]],cv[faces[q][4]]), 6)
        if !has_negative || !has_positive
            volstate = has_negative ? Int8(-1) : Int8(1)
            !has_negative && (bg[parent] = Int8(1))
            _append3DFullCell!(physical, reference, states, parents, ids, pointers,
                                  volstate, parent, g.axes, indices)
            if has_negative
                aligned = Float64[]; offsets = Int32[0]; conn = Int32[]; types = Int32[]
                for (face, verts) in enumerate(faces)
                    all(a -> iszero(cv[a]), verts) || continue
                    for a in verts
                        append!(aligned, (Float64(corners[a][1]),Float64(corners[a][2]),Float64(corners[a][3])))
                    end
                    baseid = length(aligned) ÷ 3 - 4
                    append!(conn, Int32[baseid,baseid+1,baseid+2,baseid,baseid+2,baseid+3])
                    append!(offsets, (Int32(length(conn)-3), Int32(length(conn))))
                    append!(types, (Int32(2),Int32(2)))
                end
                _append3DBoundary!(bphysical, breference, bparents, bids, bptrs, normals,
                    (aligned,offsets,conn,types), parent, cv, g.axes, indices; orient_by_levelset=true)
            end
        else
            # Loaded here, outside the guard, so a missing library is reported as itself and
            # not as a failure on whichever background cell hit the first ccall. Only cells
            # that straddle the interface need the backend at all.
            _library()
            outputs = try _cutHex(cv) catch err
                throw(ErrorException("CutCells failed on background cell $parent: $(sprint(showerror, err))"))
            end
            _append3DVolume!(physical, reference, states, parents, ids, pointers,
                               outputs[1], Int8(-1), parent, g.axes, indices)
            _append3DVolume!(physical, reference, states, parents, ids, pointers,
                               outputs[3], Int8(1), parent, g.axes, indices)
            _append3DBoundary!(bphysical, breference, bparents, bids, bptrs, normals,
                outputs[2], parent, cv, g.axes, indices)
        end
    end
    return _cutFromNative(g, level_id; bg_states=bg, subcell_states=states,
        cell_to_bgcell=parents, point_ids=ids, point_pointers=pointers,
        physical_points=physical, reference_points=reference, cell_facet_states=facets,
        boundary_facet_to_bgcell=bparents, boundary_point_ids=bids,
        boundary_point_pointers=bptrs, boundary_physical_points=bphysical,
        boundary_reference_points=breference, boundary_normals=normals)
end

"""    _fromField(field, level_id)

Produce the arrays of a `Cut` from a level-set Field, cell by cell, delegating the straddling
cells to the CutCells library.
"""
function _fromField(field::Field, level_id::Symbol)
    d = field.discretization
    _levelgrid(d, level_id)
    _levelvalues(field, level_id)
    package_grid = d.levels[level_id]
    length(package_grid.axes) == 3 && return _fromField3D(field, level_id)
    length(package_grid.axes) == 2 || throw(ArgumentError("CutCells adapter supports only two- or three-dimensional grids"))
    nx, ny = length(package_grid.axes[1]) - 1, length(package_grid.axes[2]) - 1
    hx, hy = _uniformSpacing(package_grid.axes[1]), _uniformSpacing(package_grid.axes[2])
    (isfinite(hx) && isfinite(hy) && hx > 0 && hy > 0) ||
        throw(ArgumentError("CutCells adapter requires finite, increasing grid axes"))

    raw = field.data[level_id]
    size(raw) == (nx + 1, ny + 1) || throw(DimensionMismatch(
        "Field has size $(size(raw)); grid requires $((nx + 1, ny + 1))"))
    eltype(raw) <: Real || throw(ArgumentError("CutCells adapter requires real-valued Field data"))
    all(isfinite, raw) || throw(ArgumentError("Field level set must contain only finite values"))
    phi = Float64.(raw)
    all(isfinite, phi) || throw(ArgumentError("Field level set values must be representable as finite Float64 numbers"))

    ncells = nx * ny
    bg_states = Vector{Int8}(undef, ncells)
    facet_states = Vector{NTuple{4,Int8}}(undef, ncells)
    subcell_states, cell_to_bgcell = Int8[], Int32[]
    point_ids, point_pointers = Int32[], Int32[1]
    physical_points, reference_points = NTuple{2,Float64}[], NTuple{2,Float64}[]
    boundary_parents, boundary_ids, boundary_pointers = Int32[], Int32[], Int32[1]
    boundary_physical, boundary_reference = NTuple{2,Float64}[], NTuple{2,Float64}[]
    boundary_normals = NTuple{2,Float64}[]

    for j in 1:ny, i in 1:nx
        parent = i + (j - 1) * nx
        cell_values = (phi[i, j], phi[i + 1, j], phi[i + 1, j + 1], phi[i, j + 1])
        state = _cellState(cell_values)
        bg_states[parent] = state
        facet_states[parent] = (_facetState(cell_values[1], cell_values[4]),
                                _facetState(cell_values[2], cell_values[3]),
                                _facetState(cell_values[1], cell_values[2]),
                                _facetState(cell_values[4], cell_values[3]))
        x0, y0 = package_grid.axes[1][i], package_grid.axes[2][j]
        if state == -1 || state == 1
            _appendFullCell!(physical_points, reference_points, subcell_states,
                cell_to_bgcell, point_ids, point_pointers, state, parent, x0, y0, hx, hy)
        else
            has_negative, has_positive = any(<(0.0), cell_values), any(>(0.0), cell_values)
            has_negative || has_positive || throw(ArgumentError(
                "level set is identically zero on background cell $parent; a whole-cell interface is unsupported"))
            if !(has_negative && has_positive)
                volume_state = has_negative ? Int8(-1) : Int8(1)
                if !has_negative
                    bg_states[parent] = Int8(1)
                end
                _appendFullCell!(physical_points, reference_points, subcell_states,
                    cell_to_bgcell, point_ids, point_pointers, volume_state, parent, x0, y0, hx, hy)
                if has_negative
                    ref_vertices = ((0.0, 0.0), (1.0, 0.0), (1.0, 1.0), (0.0, 1.0))
                    edge_vertices = ((1, 2), (2, 3), (3, 4), (4, 1))
                    points, offsets, connectivity, types = Float64[], Int32[0], Int32[], Int32[]
                    for (a, b) in edge_vertices
                        if iszero(cell_values[a]) && iszero(cell_values[b])
                            append!(points, ref_vertices[a])
                            append!(points, ref_vertices[b])
                            append!(connectivity, (length(points) ÷ 2 - 2, length(points) ÷ 2 - 1))
                            push!(offsets, Int32(length(connectivity)))
                            push!(types, Int32(1))
                        end
                    end
                    aligned_interface = (points, offsets, connectivity, types)
                    _appendBoundaryRegion!(boundary_physical, boundary_reference,
                        boundary_parents, boundary_ids, boundary_pointers, boundary_normals,
                        aligned_interface, parent, cell_values, x0, y0, hx, hy)
                end
            else
                _library()
                outputs = try
                    _cutQuad(cell_values)
                catch err
                    throw(ErrorException("CutCells failed on background cell $parent: $(sprint(showerror, err))"))
                end
                _appendVolumeRegion!(physical_points, reference_points, subcell_states,
                    cell_to_bgcell, point_ids, point_pointers, outputs[1], Int8(-1), parent, x0, y0, hx, hy)
                _appendVolumeRegion!(physical_points, reference_points, subcell_states,
                    cell_to_bgcell, point_ids, point_pointers, outputs[3], Int8(1), parent, x0, y0, hx, hy)
                _appendBoundaryRegion!(boundary_physical, boundary_reference,
                    boundary_parents, boundary_ids, boundary_pointers, boundary_normals,
                    outputs[2], parent, cell_values, x0, y0, hx, hy)
            end
        end
    end

    return _cutFromNative(package_grid, level_id; bg_states, subcell_states, cell_to_bgcell,
        cell_facet_states=facet_states,
        point_ids, point_pointers, physical_points, reference_points,
        boundary_facet_to_bgcell=boundary_parents,
        boundary_point_ids=boundary_ids, boundary_point_pointers=boundary_pointers,
        boundary_physical_points=boundary_physical,
        boundary_reference_points=boundary_reference,
        boundary_normals)
end

end # module CutCellsBinding
