# Core Cartesian data model. Include this file inside the EvolvingDomains module.

# Monotonic identity token. A `UInt64` survives `deepcopy` and serialization unchanged and
# keeps `Grid` a concrete (type-stable) struct, unlike the previous `Ref`-in-an-`Any` field
# which minted a fresh identity on every copy.
const _MESH_COUNTER = Threads.Atomic{UInt64}(0)
_nextMeshId() = Threads.atomic_add!(_MESH_COUNTER, UInt64(1)) + UInt64(1)

"""
    Grid{D,T}(axes::NTuple{D,Vector{T}})

Inner constructor: create a `D`-dimensional grid from one coordinate vector per dimension.
"""
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

Outer constructor: create a grid from a tuple of axes, one coordinate vector per dimension.
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
    Discretization{D,T}

A dictionary of named grids.

Discretizations are mutable: you can add grids directly, including grids that do not
share the same domain.
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

Initialize a `Discretization` from a single grid. The input grid is indexed as `:init`
in the dictionary.
"""
Discretization(axes::AbstractVector...) = Discretization(Dict(:init => Grid(axes)))

"""
    Field{T,N,D}

Stores `N`-dimensional nodal values at the levels of a `Discretization`. A `Field` is meant
to represent the same field over the different levels. The struct is a container; accessing
a level that is not stored requires an interpolator.
"""
struct Field{T,N,D}
    discretization::D
    data::Dict{Symbol,Array{T,N}}

    Field{T,N,D}(discretization::D, data::Dict{Symbol,Array{T,N}}) where {T,N,D} =
        new{T,N,D}(discretization, data)
end

"""
    Field(d::Discretization, data::AbstractDict)

Initialize a `Field` from a `Discretization` and a dictionary of nodal-value arrays sharing
type and rank, keyed by keys that already exist in the discretization.
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

"""    _levelgrid(d::Discretization, level::Symbol)

Return the grid for `level`, or throw if the level is absent."""
function _levelgrid(d::Discretization, level::Symbol)
    grid = get(d.levels, level, nothing)
    grid === nothing &&
        throw(ArgumentError("level $level is not present in the Discretization"))
    return grid
end

"""    _levelvalues(field::Field, level::Symbol)

Return the array attached to `level`, or throw if no array is attached."""
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
    Discretization(levels::AbstractDict)

Initialize a `Discretization` from an already existing dictionary of grids of common
dimension and axis endpoints.
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
    refine!(d::Discretization, level::Symbol; factor=2, name=nothing)

Provided a discretization level, generate another level with a given refinement ratio.
Optionally, the refined level can be named. The return value is the new level id.
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
    coarsen!(d::Discretization, level::Symbol; factor=2, name=nothing)

Provided a discretization level, generate another level with a given coarsening ratio.
Optionally, the coarsened level can be named. The return value is the new level id.
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

"""    attach!(field::Field{T,N}, level::Symbol, data) -> field

Attach rank-`N` nodal values to a discretization level of `field`. Overwrites any data already
attached to that level."""
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

"""    detach!(field::Field, level::Symbol) -> field

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

"""    immediateCoarser(d::Discretization, level::Symbol) -> Symbol or nothing

The next level down in node count, or `nothing` at the coarsest level. Equal resolutions are
ordered by name, so the neighbour may be equally resolved; whether data can move between the
two is `transfer!`'s question."""
immediateCoarser(d::Discretization, level::Symbol) = _neighbour(d, level, -1)

"""    immediateFiner(d::Discretization, level::Symbol) -> Symbol or nothing

The next level up in node count, or `nothing` at the finest level; equal resolutions are
ordered by name."""
immediateFiner(d::Discretization, level::Symbol) = _neighbour(d, level, +1)

"""    coarsest(d::Discretization) -> Symbol

Name of the level with the fewest nodes, breaking ties by name."""
coarsest(d::Discretization) = first(_byResolution(d))

"""    finest(d::Discretization) -> Symbol

Name of the level with the most nodes, breaking ties by name."""
finest(d::Discretization) = last(_byResolution(d))
