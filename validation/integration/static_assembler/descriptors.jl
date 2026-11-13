# StaticAssembler-specific records. Intentionally outside the package API.
struct MeshDescriptor{T}
    mesh_id::Any
    nx::Int
    ny::Int
    origin_x::T
    origin_y::T
    hx::T
    hy::T
end

struct BoundaryDescriptor{T}
    facet_to_bgcell::Vector{Int32}
    point_ids::Vector{Int32}
    point_pointers::Vector{Int32}
    normals::Vector{NTuple{2,T}}
    physical_points::Vector{NTuple{2,T}}
    reference_points::Vector{NTuple{2,T}}
end

struct CutDescriptor{T}
    mesh_id::Any
    bg_states::Vector{Int8}
    subcell_states::Vector{Int8}
    cell_to_bgcell::Vector{Int32}
    point_ids::Vector{Int32}
    point_pointers::Vector{Int32}
    physical_points::Vector{NTuple{2,T}}
    reference_points::Vector{NTuple{2,T}}
    boundary::BoundaryDescriptor{T}
end

function to_descriptors(grid, cut)
    cut.mesh_id === grid.mesh_id || throw(ArgumentError("Cut and Grid mesh identities differ"))
    length(grid.axes) == 2 || throw(ArgumentError("descriptor conversion requires a 2D Grid"))
    nx, ny = length(grid.axes[1]) - 1, length(grid.axes[2]) - 1
    hx = grid.axes[1][2] - grid.axes[1][1]
    hy = grid.axes[2][2] - grid.axes[2][1]
    geometry_type(points) = isempty(points) ? eltype(grid.axes[1]) : typeof(first(points)[1])
    T = promote_type(eltype(grid.axes[1]), eltype(grid.axes[2]),
        geometry_type(cut.physical_points), geometry_type(cut.boundary.physical_points),
        geometry_type(cut.boundary.normals))
    mesh = MeshDescriptor{T}(grid.mesh_id, nx, ny, T(first(grid.axes[1])),
        T(first(grid.axes[2])), T(hx), T(hy))
    b = cut.boundary
    boundary = BoundaryDescriptor{T}(copy(b.facet_to_bgcell), copy(b.point_ids),
        copy(b.point_pointers), copy(b.normals), copy(b.physical_points),
        copy(b.reference_points))
    cells = findall(k -> cut.bg_states[cut.cell_to_bgcell[k]] == 0,
        eachindex(cut.subcell_states))
    ids, pointers = Int32[], Int32[1]
    for k in cells
        append!(ids, @view cut.point_ids[cut.point_pointers[k]:cut.point_pointers[k + 1] - 1])
        push!(pointers, length(ids) + 1)
    end
    cutdata = CutDescriptor{T}(cut.mesh_id, copy(cut.bg_states), cut.subcell_states[cells],
        cut.cell_to_bgcell[cells], ids, pointers, copy(cut.physical_points),
        copy(cut.reference_points), boundary)
    return mesh, cutdata
end
