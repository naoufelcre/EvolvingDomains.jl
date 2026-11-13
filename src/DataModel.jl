# Core Cartesian data model. Include this file inside the EvolvingDomains module.

# Monotonic identity token. A `UInt64` survives `deepcopy` and serialization unchanged and
# keeps `Grid` a concrete (type-stable) struct, unlike the previous `Ref`-in-an-`Any` field
# which minted a fresh identity on every copy.
const _MESH_COUNTER = Threads.Atomic{UInt64}(0)
_nextMeshId() = Threads.atomic_add!(_MESH_COUNTER, UInt64(1)) + UInt64(1)

"""A Cartesian-product grid whose axes contain strictly increasing coordinates."""
struct Grid{D,T}
    axes::NTuple{D,Vector{T}}
    mesh_id::UInt64

    function Grid{D,T}(axes::NTuple{D,Vector{T}}) where {D,T}
        D > 0 || throw(ArgumentError("a grid must have at least one axis"))
        for (i, axis) in enumerate(axes)
            length(axis) >= 2 || throw(ArgumentError("grid axis $i must contain at least two coordinates"))
            all(isfinite, axis) || throw(ArgumentError("grid axis $i must contain finite coordinates"))
            all(diff(axis) .> zero(T)) || throw(ArgumentError("grid axis $i must be strictly increasing"))
        end
        new{D,T}(axes, _nextMeshId())
    end
end

"""
    Grid(axes::Tuple)

Create a structured Cartesian grid from one or more coordinate vectors. Each grid
gets a fresh `mesh_id` token; equal coordinates do not imply mesh identity.
"""
function Grid(axes::Tuple)
    isempty(axes) && throw(ArgumentError("a grid must have at least one axis"))
    all(axis -> axis isa AbstractVector, axes) ||
        throw(ArgumentError("each grid axis must be a coordinate vector"))
    all(axis -> eltype(axis) <: Real, axes) ||
        throw(ArgumentError("grid coordinates must be real numbers"))
    T = promote_type(map(eltype, axes)...)
    T <: Real || throw(ArgumentError("grid coordinates must have a common real type"))
    converted = map(axis -> Vector{T}(axis), axes)
    Grid{length(converted),T}(converted)
end

Grid(axes::AbstractVector...) = Grid(axes)

# Identity, not coordinates: two grids over the same points are different meshes by
# construction, which is what `Cut`, `TransportMap` and aggregation rely on. All three
# methods must agree, or `==`, `Set` and `Dict` contradict each other (a `deepcopy` keeps
# the mesh identity, so it must compare equal — it is the same mesh).
Base.:(==)(a::Grid, b::Grid) = a.mesh_id == b.mesh_id
Base.isequal(a::Grid, b::Grid) = isequal(a.mesh_id, b.mesh_id)
Base.hash(g::Grid, h::UInt) = hash(g.mesh_id, h)

"""
A named collection of compatible-dimensional grids over one physical box.

# Known limitation: `levels` is a plain mutable `Dict`

The shared-domain invariant is enforced in the constructor, so every `Discretization` that
is *built* has it. `levels` itself stays a mutable `Dict` because `Field` holds its
discretization **by value**: an immutable topology could not propagate a refinement to
Fields that already exist, and every alternative (a shared handle inside `Field`, manual
rebinding, or a validated choke point that users can still bypass) trades that ergonomics
for enforcement. We chose ergonomics and documented the hole.

Consequence: assigning `d.levels[name] = foreign_grid` directly bypasses the check. Inside
the package only `refine!` and `coarsen!` add levels, and both preserve the box by
construction, so the hole is reachable only by hand. A hand-grafted level is still a
self-consistent single-level problem, so `cut` and `eulerianTransport!` will run on it; the
mismatch surfaces only where two levels must agree — `transfer!` refuses it.
"""
struct Discretization{D,T}
    levels::Dict{Symbol,Grid{D,T}}

    # The invariant lives in the inner constructor so that no outer path -- including the
    # implicit `Discretization(levels::Dict{Symbol,Grid{D,T}})` Julia would otherwise
    # generate -- can build a discretization whose levels span different boxes.
    function Discretization{D,T}(levels::Dict{Symbol,Grid{D,T}}) where {D,T}
        isempty(levels) &&
            throw(ArgumentError("a discretization must contain at least one level"))
        reference_id = first(keys(levels))
        reference_ends = _axisEnds(levels[reference_id])
        for (id, grid) in levels
            _axisEnds(grid) == reference_ends || throw(ArgumentError(
                "level $id spans $(_axisEnds(grid)) but level $reference_id spans " *
                "$reference_ends; all levels of a Discretization must share one physical domain"))
        end
        new{D,T}(levels)
    end
end

"""
    Discretization(axes::AbstractVector...)

Single-level shorthand: `Discretization(range(0,1,9), range(0,1,9))` is
`Discretization(Dict(:init => Grid((range(0,1,9), range(0,1,9)))))`. Further levels are
added with `refine!`/`coarsen!`.
"""
Discretization(axes::AbstractVector...) = Discretization(Dict(:init => Grid(axes)))

"""A mutable collection of dense arrays, one per attached grid level."""
struct Field{T,N,D}
    discretization::D
    data::Dict{Symbol,Array{T,N}}

    Field{T,N,D}(discretization::D, data::Dict{Symbol,Array{T,N}}) where {T,N,D} =
        new{T,N,D}(discretization, data)
end

"""
    Field(d::Discretization, data)

Attach dense arrays to named levels. The arrays must share an element type and
rank, but their dimensions and rank need not match grid axes: layout conventions
belong to numerical methods. Arrays are retained (not copied), so a raw
`field.data[id] = array` keeps your array but converts its element type and does not check
its shape; `attach!` is the validated path. Retained means
`field.data[id]` is mutable shared state.
"""
function Field(d::Discretization, data::AbstractDict)
    isempty(data) && throw(ArgumentError("a field must contain values for at least one level"))
    all(id -> id isa Symbol, keys(data)) || throw(ArgumentError("field level IDs must be Symbols"))
    first_id = first(keys(data))
    haskey(d.levels, first_id) || throw(ArgumentError("field references unknown level $first_id"))
    first_array = data[first_id]
    first_array isa Array || throw(ArgumentError("field values must be dense N-dimensional Arrays"))
    T, N = eltype(first_array), ndims(first_array)
    typed_data = Dict{Symbol,Array{T,N}}()
    for (id, array) in data
        haskey(d.levels, id) || throw(ArgumentError("field references unknown level $id"))
        array isa Array{T,N} || throw(ArgumentError("all field values must be Arrays with the same element type and dimension"))
        typed_data[id] = array
    end
    Field{T,N,typeof(d)}(d, typed_data)
end


# Pending by design: refine!/coarsest/finest and nonuniform-grid numerical methods.

# =============================================================================
# Level topology
#
# A `Discretization` is an ordered-by-resolution collection of grids over one physical
# box. Refinement and coarsening create *levels* (grids); they never touch `Field` data,
# which is what `attach!` and `transfer!` are for. The three responsibilities are kept
# apart on purpose: topology, data wiring, algebra.

# The single tolerance for "these axes are the same uniform nodes". Uniformity is the
# backward-error bound of the node formula `x_a(i) = x_1 + (i-1)h`, so it scales with the axis
# magnitude, not with its extent: a box far from the origin has coarser absolute resolution.
# `_gridMetrics`, the cut producer and `_nested` all call `_coordinateTolerance`, so a grid is
# accepted by every operator or by none. It lives here, with the data model, because the
# standalone native runner loads this file without the geometric operators.
const _UNIFORM_RTOL = 32eps(Float64)

@inline _coordinateTolerance(first_, last_, h) =
    _UNIFORM_RTOL * max(abs(Float64(first_)), abs(Float64(last_)), abs(Float64(h)))

_nodeDims(g::Grid) = ntuple(k -> length(g.axes[k]), length(g.axes))

# Factors coarsen! would accept on an axis of n nodes: divisors of (n-1) that are at least 2.
_validFactors(n::Integer) = [f for f in 2:max(2, n - 1) if (n - 1) % f == 0]

"""
    _levelgrid(d, level) -> Grid

The single entry point for *using* a level's topology: return its `Grid` or throw. Every
operator goes through this pair, so the strict half of the contract lives in one place.

The permissive half is the container itself: `d.levels[name] = grid` is not intercepted.
Mutating `levels` by hand is allowed; it is only ever validated here, when a level is used.
"""
function _levelgrid(d::Discretization, level::Symbol)
    grid = get(d.levels, level, nothing)
    grid === nothing &&
        throw(ArgumentError("level $level is not present in the Discretization"))
    return grid
end

"""
    _levelvalues(field, level) -> Array

The single entry point for *using* a level's data: return the attached array or throw.
"""
function _levelvalues(field::Field, level::Symbol)
    values = get(getproperty(field, :data), level, nothing)
    values === nothing &&
        throw(ArgumentError("level $level has no attached values; use attach! first"))
    return values
end
_resolution(g::Grid) = prod(_nodeDims(g))

# The invariant that makes "level" mean something: one physical box for every level.
# Resolution stays free, so unmatched levels (17² beside 33²) remain legal.
_axisEnds(g::Grid) = ntuple(k -> (first(g.axes[k]), last(g.axes[k])), length(g.axes))

"""
    Discretization(levels)

Create a stable Symbol-keyed collection of grids that share one physical domain: every
level must start and end on the same coordinate along each axis. Resolution is free, so
levels of different refinement coexist; a genuinely different box is rejected because
`transfer!` and any multilevel algorithm assume one domain.
"""
function Discretization(levels::AbstractDict)
    all(id -> id isa Symbol, keys(levels)) || throw(ArgumentError("level IDs must be Symbols"))
    first_grid = first(values(levels))
    first_grid isa Grid || throw(ArgumentError("each level must contain a Grid"))
    D, T = length(first_grid.axes), eltype(first_grid.axes[1])
    typed_levels = Dict{Symbol,Grid{D,T}}()
    for (id, g) in levels
        g isa Grid{D,T} || throw(ArgumentError(
            "all grids must have the same dimension and coordinate type"))
        typed_levels[id] = g
    end
    Discretization{D,T}(typed_levels)   # inner constructor enforces the shared domain
end

"""Refine one axis by inserting `factor - 1` nodes in each interval, endpoints kept."""
function _refinedAxis(axis::Vector{T}, factor::Integer) where {T}
    refined = Vector{T}(undef, (length(axis) - 1) * factor + 1)
    index = 1
    for i in 1:(length(axis) - 1)
        a, b = axis[i], axis[i + 1]
        for k in 0:(factor - 1)
            refined[index] = T(a + (b - a) * (k // factor))
            index += 1
        end
    end
    refined[index] = axis[end]
    return refined
end

function _freshLevelName(d::Discretization, base::Symbol, tag::Symbol, factor::Integer,
                         name)
    if name !== nothing
        name isa Symbol || throw(ArgumentError("a level name must be a Symbol"))
        haskey(d.levels, name) && throw(ArgumentError(
            "level $name already exists; pass another name"))
        return name
    end
    candidate = Symbol(base, :_, tag, factor)
    suffix = 1
    while haskey(d.levels, candidate)
        suffix += 1
        candidate = Symbol(base, :_, tag, factor, :_, suffix)
    end
    return candidate
end

"""
    refine!(d, level; factor = 2, name = nothing) -> Symbol
    refine!(d, level, name::Symbol; factor = 2) -> Symbol

Add a refined level to `d`, interpolating node positions inside every interval so the
physical domain is preserved. Returns the new level's name. The new level carries no data;
attach it with `attach!` and fill it with `transfer!`. Omitting `name` derives one
(`:init_r2`) with a numeric suffix on collision.

Requires floating-point coordinates: a midpoint is not representable on a whole-number axis,
so such a grid is refused rather than silently rounded. `coarsen!` needs no fractions.
"""
function refine!(d::Discretization, level::Symbol; factor::Integer=2, name=nothing)
    g = _levelgrid(d, level)
    factor >= 2 || throw(ArgumentError("refine! requires factor >= 2; got $factor"))
    # A whole-number mesh is a legitimate container, but refinement asks for midpoints,
    # which the integers do not have. Refuse here, where the fraction is actually needed.
    T = eltype(g.axes[1])
    # Integer and rational axes are refused even when a particular midpoint happens to be
    # representable: the refinement rule is only specified for floating-point coordinates.
    T <: AbstractFloat || throw(ArgumentError(
        "refine! requires a floating-point axis but level $level stores $T; build the grid " *
        "with floating-point coordinates (for example 0.0:1.0)"))
    factor <= typemax(Int) ÷ max(1, length(g.axes[1]) - 1) || throw(ArgumentError(
        "refine! factor $factor would overflow the axis length"))
    new = _freshLevelName(d, level, :r, factor, name)
    d.levels[new] = Grid(ntuple(k -> _refinedAxis(g.axes[k], factor), length(g.axes)))
    return new
end
refine!(d::Discretization, level::Symbol, name::Symbol; factor::Integer=2) =
    refine!(d, level; factor=factor, name=name)

"""
    coarsen!(d, level; factor = 2, name = nothing) -> Symbol
    coarsen!(d, level, name::Symbol; factor = 2) -> Symbol

Add a level that keeps every `factor`-th node of `level`. `(n - 1)` must be divisible by
`factor` on every axis, otherwise the subsample could not keep both endpoints and the
shared-domain invariant would break. Returns the new level's name.
"""
function coarsen!(d::Discretization, level::Symbol; factor::Integer=2, name=nothing)
    g = _levelgrid(d, level)
    factor >= 2 || throw(ArgumentError("coarsen! requires factor >= 2; got $factor"))
    for (k, axis) in enumerate(g.axes)
        (length(axis) - 1) % factor == 0 || throw(ArgumentError(
            "coarsen! cannot keep both endpoints of axis $k: $(length(axis) - 1) intervals " *
            "is not divisible by $factor" *
            (isempty(_validFactors(length(axis))) ? "; this axis is too short to coarsen" :
             "; divisors of $(length(axis) - 1) that are at least 2: " *
             join(_validFactors(length(axis)), ", "))))
    end
    new = _freshLevelName(d, level, :c, factor, name)
    d.levels[new] = Grid(ntuple(k -> g.axes[k][1:factor:end], length(g.axes)))
    return new
end
coarsen!(d::Discretization, level::Symbol, name::Symbol; factor::Integer=2) =
    coarsen!(d, level; factor=factor, name=name)

"""
    attach!(field, level, data) -> field

Attach `data` to an existing level of `field`'s discretization. The array is **retained**,
not copied, matching `Field`'s ownership model; element type, rank and shape are checked
first so that a level cannot hold an array its grid does not describe.
"""
function attach!(field::Field{T,N}, level::Symbol, data) where {T,N}
    g = _levelgrid(getproperty(field, :discretization), level)
    length(g.axes) == N || throw(ArgumentError(
        "level $level is $(length(g.axes))-dimensional but the Field holds rank-$N arrays"))
    data isa Array{T,N} || throw(ArgumentError(
        "level $level expects an Array{$T,$N}; got $(typeof(data))"))
    size(data) == _nodeDims(g) || throw(DimensionMismatch(
        "level $level expects size $(_nodeDims(g)); got $(size(data))"))
    getproperty(field, :data)[level] = data
    return field
end

"""    detach!(field, level) -> field

Remove the array attached to `level`, leaving the level in the discretization."""
function detach!(field::Field, level::Symbol)
    data = getproperty(field, :data)
    haskey(data, level) || throw(ArgumentError("level $level has no attached values"))
    delete!(data, level)
    return field
end

# Total order on levels: fewest nodes first, ties broken by name. Ties are ordinary --
# coarsening back onto an existing resolution is a normal workflow -- so they are ordered
# deterministically rather than rejected, and navigation always returns a neighbour.
function _byResolution(d::Discretization)
    sort!(collect(keys(d.levels));
          by=id -> (_resolution(d.levels[id]), string(id)))
end

function _neighbour(d::Discretization, level::Symbol, step::Int)
    _levelgrid(d, level)
    order = _byResolution(d)
    i = findfirst(==(level), order)
    j = i + step
    1 <= j <= length(order) || return nothing
    return order[j]
end

"""    immediateCoarser(d, level) -> Symbol or nothing

The next level down in node count, or `nothing` at the coarsest level. Levels of equal
resolution are ordered by name, so the neighbour can be equally resolved; whether data can
actually move between the two levels is `transfer!`'s question, not this one."""
immediateCoarser(d::Discretization, level::Symbol) = _neighbour(d, level, -1)

"""    immediateFiner(d, level) -> Symbol or nothing

The next level up in node count, or `nothing` at the finest level. Levels of equal
resolution are ordered by name, so the neighbour can be equally resolved."""
immediateFiner(d::Discretization, level::Symbol) = _neighbour(d, level, +1)

"""    coarsest(d) -> Symbol

Name of the level with the fewest nodes."""
coarsest(d::Discretization) = first(_byResolution(d))

"""    finest(d) -> Symbol

Name of the level with the most nodes."""
finest(d::Discretization) = last(_byResolution(d))
