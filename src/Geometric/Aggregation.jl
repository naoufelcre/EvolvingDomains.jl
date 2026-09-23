import GridapEmbedded.AgFEM: _find_best_neighbor, _touch_aggregated_cells!
using Gridap.Arrays: array_cache
using GridapEmbedded.Interfaces: EmbeddedDiscretization, IN, CUT
using Gridap.Geometry: get_cell_coordinates, get_grid_topology, get_faces,
                       get_triangulation, num_cell_dims

"""
    aggregate_cut_cells(cutgeo)

Aggregate every cut cell to a valid root. Use fully interior roots when
available; for an all-cut component, use its largest physical cell.
"""
function aggregate_cut_cells(cutgeo::EmbeddedDiscretization)
    raw_states = cutgeo.ls_to_bgcell_to_inoutcut
    states = eltype(raw_states) <: AbstractVector ? only(raw_states) : raw_states
    facet_states = GridapEmbedded.compute_bgfacet_to_inoutcut(cutgeo.bgmodel, cutgeo.geo)
    triangulation = get_triangulation(cutgeo.bgmodel)
    topology = get_grid_topology(cutgeo.bgmodel)
    dimension = num_cell_dims(cutgeo.bgmodel)
    cell_to_faces = get_faces(topology, dimension, dimension - 1)
    face_to_cells = get_faces(topology, dimension - 1, dimension)
    _aggregate_all_cut_cells(
        cutgeo, states, facet_states, get_cell_coordinates(triangulation),
        cell_to_faces, face_to_cells)
end

function _physical_cell_areas(cutgeo, n_cells)
    areas = zeros(Float64, n_cells)
    subcells = cutgeo.subcells
    connectivity = subcells.cell_to_points
    subcell_states = only(cutgeo.ls_to_subcell_to_inout)
    for subcell in eachindex(subcells.cell_to_bgcell)
        subcell_states[subcell] == IN || continue
        first_point = Int(connectivity.ptrs[subcell])
        next_point = Int(connectivity.ptrs[subcell + 1])
        next_point == first_point + 3 ||
            error("cut-cell aggregation requires triangular cut subcells")
        p1 = subcells.point_to_coords[connectivity.data[first_point]]
        p2 = subcells.point_to_coords[connectivity.data[first_point + 1]]
        p3 = subcells.point_to_coords[connectivity.data[first_point + 2]]
        areas[subcells.cell_to_bgcell[subcell]] += abs(
            (p2[1] - p1[1]) * (p3[2] - p1[2]) -
            (p3[1] - p1[1]) * (p2[2] - p1[2])) / 2
    end
    areas
end

function _aggregate_all_cut_cells(
    cutgeo, states, facet_states, cell_coordinates, cell_to_faces, face_to_cells)

    n_cells = length(states)
    roots = zeros(Int32, n_cells)
    touched = falses(n_cells)
    for cell in 1:n_cells
        if states[cell] == IN
            roots[cell] = cell
            touched[cell] = true
        end
    end

    face_cache = array_cache(cell_to_faces)
    neighbor_cache = array_cache(face_to_cells)
    cell_cache = array_cache(cell_coordinates)
    root_cache = array_cache(cell_coordinates)
    fallback_roots = Int[]
    areas = nothing

    for _ in 1:n_cells
        unaggregated = 0
        made_progress = false
        for cell in 1:n_cells
            if !touched[cell] && states[cell] == CUT
                neighbor = _find_best_neighbor(
                    face_cache, neighbor_cache, cell_cache, root_cache, cell,
                    cell_to_faces, face_to_cells, cell_coordinates, touched,
                    roots, facet_states, IN)
                if neighbor > 0
                    roots[cell] = roots[neighbor]
                    made_progress = true
                else
                    unaggregated += 1
                end
            end
        end
        unaggregated == 0 && break
        _touch_aggregated_cells!(touched, roots)

        if !made_progress
            isnothing(areas) && (areas = _physical_cell_areas(cutgeo, n_cells))
            root = 0
            best_area = 0.0
            for cell in 1:n_cells
                if !touched[cell] && states[cell] == CUT && areas[cell] > best_area
                    root = cell
                    best_area = areas[cell]
                end
            end
            root > 0 || break
            roots[root] = root
            touched[root] = true
            push!(fallback_roots, root)
        end
    end

    unresolved = findall(i -> roots[i] == 0 && states[i] == CUT, 1:n_cells)
    isempty(unresolved) || error(
        "cut-cell aggregation failed for $(length(unresolved)) cells: " *
        join(unresolved[1:min(10, length(unresolved))], ", "))
    isempty(fallback_roots) || @warn(
        "AgFEM component has no interior cell; using its largest cut cell as root",
        roots=fallback_roots, maxlog=3)
    roots
end
