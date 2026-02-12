module Geometric

using StaticArrays: SVector
using GridapEmbedded
using Gridap: CartesianDiscreteModel
using Gridap.Geometry: get_node_coordinates, get_grid_topology, get_faces, num_nodes
using GridapEmbedded: cut
using Gridap.TensorValues: VectorValue

include("GridInfo.jl")
include("CartesianField.jl")
include("Stencils.jl")
include("GeometryDesign.jl")

using .CartesianField
using .Stencils

export CartesianMeshField, meshsize, get_interpolator
export D⁺, D⁻, D⁰, D²
export weno5⁻, weno5⁺

export GeometryCache, EvolvingDiscreteGeometry
export set_levelset!, current_levelset, current_cut, get_active_indices
export CartesianGridInfo, grid_info
export AbstractGeometry, Circle, Rectangle, Translate, signed_distance


mutable struct GeometryCache
    transfer_op::Any
    extension_op::Any
    cut::Union{Nothing,GridapEmbedded.Interfaces.EmbeddedDiscretization}
    prev_cut::Union{Nothing,GridapEmbedded.Interfaces.EmbeddedDiscretization}

    function GeometryCache()
        new(nothing, nothing, nothing, nothing)
    end
end

mutable struct EvolvingDiscreteGeometry
    levelset::Vector{Float64}
    grid::CartesianDiscreteModel
    cache::GeometryCache
end

"""
    EvolvingDiscreteGeometry(levelset::Vector{Float64}, grid::CartesianDiscreteModel)

Create a geometry object from a levelset vector and a grid.
The cache is initialized automatically.
"""
function EvolvingDiscreteGeometry(levelset::Vector{Float64}, grid::CartesianDiscreteModel)
    return EvolvingDiscreteGeometry(levelset, grid, GeometryCache())
end

"""
    EvolvingDiscreteGeometry(grid::CartesianDiscreteModel)

Initialize an empty geometry (all zeros) for a given grid.
Useful for cases where set_levelset! will be called immediately after.
"""
function EvolvingDiscreteGeometry(grid::CartesianDiscreteModel)
    n_nodes = num_nodes(grid)
    phi = zeros(n_nodes)
    return EvolvingDiscreteGeometry(phi, grid)
end

"""
    current_levelset(geom::EvolvingDiscreteGeometry) -> Vector{Float64}

Return the current level set field.
"""
current_levelset(geom::EvolvingDiscreteGeometry) = geom.levelset

"""
    current_cut(geom::EvolvingDiscreteGeometry)

Return the current GridapEmbedded cut geometry.
"""
current_cut(geom::EvolvingDiscreteGeometry) = geom.cache.cut

"""
    set_levelset!(geom::EvolvingDiscreteGeometry, ϕ::Vector{Float64})

Update the level set function and invalidate the geometry cache.
Previous geometry cut is moved to `prev_cut` for history access.
"""
function set_levelset!(geom::EvolvingDiscreteGeometry, ϕ::Vector{Float64})
    @assert length(ϕ) == length(geom.levelset) "Dimension mismatch"
    copyto!(geom.levelset, ϕ)

    # Shift current cut to history
    geom.cache.prev_cut = geom.cache.cut
    
    # Invalidate current cache
    geom.cache.cut = nothing
    geom.cache.transfer_op = nothing
    geom.cache.extension_op = nothing
    return geom
end

"""
    get_active_indices(geom::EvolvingDiscreteGeometry, state::Symbol=:current) -> Vector{Int}

Return the linear indices of the active **nodes** (nodes belonging to IN or CUT cells) for the specified state.
State can be `:current` or `:prev`.
"""
function get_active_indices(geom::EvolvingDiscreteGeometry, state::Symbol=:current)
    cut_geo = nothing
    if state == :current
        if isnothing(geom.cache.cut)
            # We need to wrap the levelset vector into a geometry object that 'cut' accepts.
            flat_coords = vec(collect(get_node_coordinates(geom.grid)))
            geo = GridapEmbedded.LevelSetCutters.DiscreteGeometry(geom.levelset, flat_coords)
            geom.cache.cut = cut(geom.grid, geo)
        end
        cut_geo = geom.cache.cut
    elseif state == :prev
        cut_geo = geom.cache.prev_cut
        if isnothing(cut_geo)
            return Int[]
        end
    else
        error("Unknown state: $state. Use :current or :prev.")
    end

    # 1. Access Cell Classification
    const_IN = -1
    const_CUT = 0
    
    raw_status = cut_geo.ls_to_bgcell_to_inoutcut
    status = if eltype(raw_status) <: AbstractVector
        raw_status[1]
    else
        raw_status
    end
    
    active_cell_indices = findall(x -> x == const_IN || x == const_CUT, status)
    
    if isempty(active_cell_indices)
        return Int[]
    end

    # 2. Map Cells -> Nodes (Topology)
    grid = geom.grid
    topo = get_grid_topology(grid)
    cell_to_nodes = get_faces(topo, 2, 0)
    
    active_nodes = Int[]
    sizehint!(active_nodes, length(active_cell_indices) * 4)
    for cell_idx in active_cell_indices
        nodes = cell_to_nodes[cell_idx]
        append!(active_nodes, nodes)
    end
    sort!(unique!(active_nodes))
    return active_nodes
end

include("Reinitialization.jl")
using .Reinitialization
export reinitialize!

end # module