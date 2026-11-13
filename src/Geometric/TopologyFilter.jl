module TopologyFilter

using ..Geometric: _checklevel, _gridMetrics
using ....EvolvingDomains: Field, _levelgrid, _levelvalues

"""
    filterSmallPhaseIslands!(phi::AbstractArray; phase = :negative, min_nodes = 5,
                             magnitude = one(eltype(phi)))
    filterSmallPhaseIslands!(field::Field, level::Symbol; phase = :negative,
                             min_nodes = 5, magnitude = nothing)

Remove every connected component of the selected phase that has fewer than `min_nodes`
nodes, in place, by overwriting its nodes with a value of the opposite sign:
`+magnitude` for `phase = :negative`, `-abs(magnitude)` for `phase = :positive`. Component
connectivity is face adjacency on the node graph — four neighbours in 2D, six in 3D.
Returns `(n_components, n_removed_components, n_flipped_nodes)`.

**A field-hygiene operation, not an accuracy improvement.** It deletes geometry the
discretization cannot resolve; it does not make transport better. Two consequences:

- \\(\\phi\\) is no longer a signed-distance approximation after a flip, whatever it was
  before. Signs and surviving zero crossings are what the filter preserves.
- **Build a fresh `cut` afterwards.** Any snapshot taken before the call describes the old
  geometry. `reinitialize!(field, cut(field, level))` then restores a usable level set.

Node count is the threshold, so it is resolution-dependent: five nodes is a different
physical size on every grid. Only *isolated* components are touched — a thin channel or
bridge between two large regions belongs to a large component and survives. Filtering by
feature width is a different operation, deliberately not folded in here.

`phase = :negative` covers `phi <= 0` and `phase = :positive` covers `phi > 0`, so a node
on the zero level belongs to the negative phase only and is never counted twice.

The array method works in any dimension. The `Field` method reads the level's array and
defaults `magnitude` to the finest node spacing, so it accepts only what `_gridMetrics`
accepts: a uniformly spaced 2D or 3D grid.
"""
function filterSmallPhaseIslands!(phi::AbstractArray{T,N}; phase::Symbol=:negative,
                                  min_nodes::Integer=5,
                                  magnitude::Real=one(T)) where {T<:AbstractFloat,N}
    phase in (:negative, :positive) ||
        throw(ArgumentError("phase must be :negative or :positive; got $(repr(phase))"))
    min_nodes > 0 ||
        throw(ArgumentError("min_nodes must be positive; got $min_nodes"))
    delta = T(magnitude)
    isfinite(delta) && delta > 0 ||
        throw(ArgumentError("magnitude must be positive and finite; got $magnitude"))

    # One closure, one branch: a node is in exactly one phase, zero included.
    inphase = phase === :negative ? x -> x <= 0 : x -> x > 0
    replacement = phase === :negative ? delta : -delta

    visited = falses(size(phi))
    stack = CartesianIndex{N}[]
    # Only the first `min_nodes` nodes of a component decide whether it is small, so a kept
    # component costs `O(min_nodes)` memory instead of `O(nodes)`.
    small = CartesianIndex{N}[]
    directions = ntuple(d -> CartesianIndex(ntuple(k -> Int(k == d), N)), N)

    n_components = 0
    n_removed_components = 0
    n_flipped_nodes = 0

    for start in CartesianIndices(phi)
        (visited[start] || !inphase(phi[start])) && continue

        n_components += 1
        empty!(stack)
        empty!(small)
        push!(stack, start)
        visited[start] = true
        nodes = 0

        while !isempty(stack)
            I = pop!(stack)
            nodes += 1
            length(small) < min_nodes && push!(small, I)
            for e in directions, J in (I + e, I - e)
                checkbounds(Bool, phi, J) || continue
                (visited[J] || !inphase(phi[J])) && continue
                visited[J] = true
                push!(stack, J)
            end
        end

        if nodes < min_nodes
            for I in small
                phi[I] = replacement
            end
            n_removed_components += 1
            n_flipped_nodes += nodes
        end
    end

    return (n_components=n_components, n_removed_components=n_removed_components,
            n_flipped_nodes=n_flipped_nodes)
end

function filterSmallPhaseIslands!(field::Field, level::Symbol; phase::Symbol=:negative,
                                  min_nodes::Integer=5, magnitude=nothing)
    info = _gridMetrics(_levelgrid(getproperty(field, :discretization), level))
    phi = _checklevel(_levelvalues(field, level), nothing, info, "Level-set Field values")
    delta = magnitude === nothing ? minimum(info.spacing) : magnitude
    return filterSmallPhaseIslands!(phi; phase=phase, min_nodes=min_nodes,
                                    magnitude=delta)
end

# =============================================================================
# Sub-resolution feature filtering
#
# The structuring element is the digital ball B_k = {v in Z^N : |v|^2 <= k^2}. Both
# morphological passes reduce to one primitive, the exact Euclidean distance transform:
# erosion keeps a node iff the nearest counter-phase node is farther than k, dilation adds
# a node iff the nearest in-phase node is within k. That makes the criterion isotropic up
# to lattice quantization instead of the sqrt(2) (2D) / sqrt(3) (3D) anisotropy a box
# element imposes at every radius.
#
# shortcone: the ball is spherical in NODE units, so on a grid with unequal spacings it is
# an ellipsoid physically; upgrade to a physical-space transform if anisotropic grids need
# the same criterion in every direction.

const _EDT_FAR = 1e20

# 1-D squared distance transform: the lower envelope of the parabolas f[q] + (i - q)^2,
# scanned once (Felzenszwalb & Huttenlocher, 2012). `f` holds 0 at a source node and
# `_EDT_FAR` elsewhere; a large finite value instead of `Inf` keeps the envelope arithmetic
# free of NaN when two unreachable sources are compared.
function _edtLine!(f::Vector{Float64}, d::Vector{Float64}, v::Vector{Int},
                   z::Vector{Float64}, n::Int)
    k = 0
    v[1] = 1
    z[1] = -Inf
    z[2] = Inf
    @inbounds for q in 2:n
        s = ((f[q] + q * q) - (f[v[k+1]] + v[k+1] * v[k+1])) / (2q - 2v[k+1])
        while s <= z[k+1]
            k -= 1
            s = ((f[q] + q * q) - (f[v[k+1]] + v[k+1] * v[k+1])) / (2q - 2v[k+1])
        end
        k += 1
        v[k+1] = q
        z[k+1] = s
        z[k+2] = Inf
    end
    k = 0
    @inbounds for q in 1:n
        while z[k+2] < q
            k += 1
        end
        d[q] = (q - v[k+1])^2 + f[v[k+1]]
    end
    return d
end

"""
    _squaredDistanceTransform(source) -> Array{Float64}

Squared Euclidean distance, in node units, from every node to the nearest `true` node of
`source`. One separability pass per axis, `O(D * nnodes)`, dimension-agnostic.
"""
function _squaredDistanceTransform(source::AbstractArray{Bool,N}) where {N}
    dims = size(source)
    d2 = Array{Float64}(undef, dims)
    @inbounds for I in CartesianIndices(d2)
        d2[I] = source[I] ? 0.0 : _EDT_FAR
    end
    nmax = maximum(dims)
    f = Vector{Float64}(undef, nmax)
    d = Vector{Float64}(undef, nmax)
    v = Vector{Int}(undef, nmax + 1)
    z = Vector{Float64}(undef, nmax + 2)
    for axis in 1:N
        n = dims[axis]
        for line in CartesianIndices(ntuple(k -> k == axis ? 1 : dims[k], N))
            @inbounds for i in 1:n
                f[i] = d2[CartesianIndex(ntuple(k -> k == axis ? i : line[k], N))]
            end
            _edtLine!(f, d, v, z, n)
            @inbounds for i in 1:n
                d2[CartesianIndex(ntuple(k -> k == axis ? i : line[k], N))] = d[i]
            end
        end
    end
    return d2
end

# Both passes see the same extended mask: the domain boundary is replicated k nodes out, so
# "outside the grid" means "same as the edge node" for erosion and dilation alike. Padding
# both is what keeps the opening anti-extensive -- eroding against a padded complement alone
# would let a positive edge node survive erosion and be dilated back as in-phase, so a node
# outside the phase would become phase, which no opening may do.
function _paddedMask(mask::AbstractArray{Bool,N}, k::Int) where {N}
    dims = size(mask)
    padded = Array{Bool}(undef, ntuple(d -> dims[d] + 2k, N))
    @inbounds for I in CartesianIndices(padded)
        padded[I] = mask[CartesianIndex(ntuple(d -> clamp(I[d] - k, 1, dims[d]), N))]
    end
    return padded
end

_croppedMask(padded::AbstractArray{Bool,N}, k::Int, dims) where {N} =
    padded[ntuple(d -> (k + 1):(k + dims[d]), N)...]

# Erosion by B_k: the ball of radius k around a kept node must be entirely in phase, so
# the nearest counter-phase node has to be strictly farther than k.
function _ballErode(mask::AbstractArray{Bool,N}, k::Int) where {N}
    dims = size(mask)
    padded = _paddedMask(mask, k)
    return _croppedMask(map(>(k * k), _squaredDistanceTransform(map(!, padded))), k, dims)
end

# Dilation by B_k: some in-phase node lies within distance k.
function _ballDilate(mask::AbstractArray{Bool,N}, k::Int) where {N}
    dims = size(mask)
    padded = _paddedMask(mask, k)
    return _croppedMask(map(<=(k * k), _squaredDistanceTransform(padded)), k, dims)
end

# Number of connected components of a mask under face adjacency.
function _componentCount(mask::AbstractArray{Bool,N}) where {N}
    visited = falses(size(mask))
    stack = CartesianIndex{N}[]
    directions = ntuple(d -> CartesianIndex(ntuple(k -> Int(k == d), N)), N)
    components = 0
    for start in CartesianIndices(mask)
        (visited[start] || !mask[start]) && continue
        components += 1
        push!(stack, start)
        visited[start] = true
        while !isempty(stack)
            I = pop!(stack)
            for e in directions, J in (I + e, I - e)
                checkbounds(Bool, mask, J) || continue
                (visited[J] || !mask[J]) && continue
                visited[J] = true
                push!(stack, J)
            end
        end
    end
    return components
end

# Background cells that straddle the interface, i.e. whose corners are not single-signed.
# The cut's cell states are a function of exactly this pattern, so it is the cut-level
# statistic that needs no snapshot.
function _cutCellCount(mask::AbstractArray{Bool,N}) where {N}
    dims = size(mask)
    corners = CartesianIndices(ntuple(_ -> 0:1, N))
    cells = 0
    for cell in CartesianIndices(ntuple(d -> dims[d] - 1, N))
        reference = mask[CartesianIndex(ntuple(d -> cell[d], N))]
        mixed = false
        for corner in corners
            mixed |= mask[CartesianIndex(ntuple(d -> cell[d] + corner[d], N))] != reference
        end
        cells += mixed
    end
    return cells
end

"""
    filterSubResolutionFeatures!(phi::AbstractArray; phase = :negative, mode = :opening,
                                 radius = 1, magnitude = one(eltype(phi)))
    filterSubResolutionFeatures!(field::Field, level::Symbol; phase = :negative,
                                 mode = :opening, radius = 1, magnitude = nothing)

Remove the geometry of one phase that is below the chosen resolution, in place, by a single
morphological pass with the digital ball of radius `radius` nodes:

- `mode = :opening` — γ_B(A) = (A ⊖ B) ⊕ B: removes islands, filaments, protrusions and
  thin bridges of the phase. Can split a region.
- `mode = :closing` — φ_B(A) = (A ⊕ B) ⊖ B: fills cavities, slits and narrow channels of
  the opposite phase. Can merge two regions.

`phase = :negative` filters `A = {φ <= 0}` and `phase = :positive` filters `A = {φ > 0}`,
so a node on the zero level belongs to the negative phase only. Changed nodes are written as
`+magnitude` / `-magnitude` (default `h`, the node spacing, for the `Field` method;
`one(T)` for arrays); every other node keeps its value bit-for-bit.

**The criterion, precisely.** `γ_B(A) = ⋃{x + B : x + B ⊆ A}`: what survives is the union of
all ball translates that fit inside the phase. It is not "delete every locally thin
section", and the difference is measurable. A thin neck is eroded away but the dilation
restores the `radius` nodes nearest each root, so a neck of nodal length `<= 2k` survives and
one longer than `2k` is severed (2D, measured). A convex corner is never coverable by a ball,
so an axis-aligned block loses its corner nodes — exactly four at `radius = 1`, whatever its
size. `B_k` is the lattice ball `{v ∈ Z^N : |v|₂ <= k}`, the same element in 2D and 3D.

**Radius 1 in 3D is rough.** The six-neighbour cross is a poor digital ball: a 5x5x5 block
loses 44 of its 125 nodes and a 1-node-thick slab disappears entirely. The anisotropy ratio
falls 1.70 -> 1.11 -> 1.06 at `k = 1 -> 3 -> 8`, so use `radius >= 2` in 3D whenever
well-resolved corners matter.

**Boundary.** Both passes see the same mask replicated `radius` nodes out of the grid, so the
domain boundary behaves as its edge node does and the opening stays anti-extensive (no node
ever enters the phase). A sub-resolution channel that reaches the boundary keeps one-node
notches at its mouths rather than being sealed, because a mouth corner cannot be covered.

**What the threshold means.** An axis-aligned band of `t` nodes is removed iff `t <= 2k`
and kept iff `t >= 2k+1`; with physical width `w = (t-1)h` that is `k = 1` removing `w <= h`,
`k = 2` removing `w <= 3h`. Obliquely the threshold is lattice-dependent: the required
perpendicular width is `2kh` in the continuum, and on a lattice the measured anisotropy ratio
(max over directions / min over directions) is 1.41 / 1.12 / 1.05 at `k = 1/3/8` in 2D and
1.70 / 1.11 / 1.06 in 3D. So `radius = 1` removes features thinner than roughly `2h` in any
orientation, up to that lattice spread.

**Consequences.** This is geometric regularization, not an accuracy improvement: it decides
which geometry the discretization is asked to resolve, and it may change connectivity by
design. It breaks the signed-distance property wherever it flips a sign, it conserves
neither area nor volume, and it cannot tell numerical junk from a genuine sub-resolution
feature — that judgement is the caller's, which is why the call is explicit and why the
statistics are returned. Build a fresh `cut` afterwards, and reinitialize if `φ` is used as
a distance:

```julia
stats = filterSubResolutionFeatures!(field, :init; mode=:closing, radius=1)
stats.n_changed_nodes > 0 && reinitialize!(field, cut(field, :init))
```

Returns `(n_flipped_to_positive, n_flipped_to_negative, n_changed_nodes,
 n_cut_cells_before, n_cut_cells_after, n_components_before, n_components_after, radius)`,
with component counts on the face-adjacency node graph of the filtered phase.

Only uniform 2D and 3D grids are accepted by the `Field` method (`_gridMetrics`); the array
method works in any dimension and ignores geometry, working in node units.
"""
function filterSubResolutionFeatures!(phi::AbstractArray{T,N}; phase::Symbol=:negative,
                                      mode::Symbol=:opening, radius::Integer=1,
                                      magnitude::Real=one(T)) where {T<:AbstractFloat,N}
    phase in (:negative, :positive) ||
        throw(ArgumentError("phase must be :negative or :positive; got $(repr(phase))"))
    mode in (:opening, :closing) ||
        throw(ArgumentError("mode must be :opening or :closing; got $(repr(mode))"))
    radius >= 1 ||
        throw(ArgumentError("radius must be at least one node; got $radius"))
    delta = T(magnitude)
    isfinite(delta) && delta > 0 ||
        throw(ArgumentError("magnitude must be positive and finite; got $magnitude"))

    negative = map(<=(zero(T)), phi)
    material = phase === :negative ? negative : map(!, negative)
    filtered = mode === :opening ?
        _ballDilate(_ballErode(material, radius), radius) :
        _ballErode(_ballDilate(material, radius), radius)
    target = phase === :negative ? filtered : map(!, filtered)

    cut_cells_before = _cutCellCount(negative)
    components_before = _componentCount(negative)

    n_to_positive = 0
    n_to_negative = 0
    for I in CartesianIndices(phi)
        negative[I] == target[I] && continue
        if target[I]
            phi[I] = -delta
            n_to_negative += 1
        else
            phi[I] = delta
            n_to_positive += 1
        end
    end

    return (n_flipped_to_positive=n_to_positive, n_flipped_to_negative=n_to_negative,
            n_changed_nodes=n_to_positive + n_to_negative,
            n_cut_cells_before=cut_cells_before, n_cut_cells_after=_cutCellCount(target),
            n_components_before=components_before,
            n_components_after=_componentCount(target), radius=Int(radius))
end

function filterSubResolutionFeatures!(field::Field, level::Symbol; phase::Symbol=:negative,
                                      mode::Symbol=:opening, radius::Integer=1,
                                      magnitude=nothing)
    info = _gridMetrics(_levelgrid(getproperty(field, :discretization), level))
    phi = _checklevel(_levelvalues(field, level), nothing, info, "Level-set Field values")
    delta = magnitude === nothing ? minimum(info.spacing) : magnitude
    return filterSubResolutionFeatures!(phi; phase=phase, mode=mode, radius=radius,
                                        magnitude=delta)
end

end # module
