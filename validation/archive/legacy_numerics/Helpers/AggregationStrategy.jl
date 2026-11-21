# Override of GridapEmbedded.AgFEM._aggregate_by_threshold_barrier
#
# Purpose:
#   On ultra-fine grids (or near pinch-off / neck singularities), the
#   aggregation chain from a CUT cell to the nearest root cell can exceed
#   the upstream hardcoded limit of 20 hops. This causes a bare
#   `@assert all_aggregated` to fire with no diagnostic information.
#
# Changes vs upstream (GridapEmbedded/src/AgFEM/CellAggregation.jl):
#   1. max_iters: 20 -> 100
#   2. @assert all_aggregated -> descriptive error(...) with cell counts
#
# Usage:
#   include("Aggregation/AggregationStrategy.jl")  # after `using GridapEmbedded`
#
# The existing call site:
#   strategy = AggregateCutCellsByThreshold(1.0)
#   aggregates = aggregate(strategy, cutgeo)
# dispatches through _aggregate_by_threshold -> _aggregate_by_threshold_barrier,
# so this override is picked up automatically.

import GridapEmbedded.AgFEM: _aggregate_by_threshold_barrier,
                             _find_best_neighbor,
                             _touch_aggregated_cells!
using Gridap.Arrays: array_cache, getindex!
using GridapEmbedded.Interfaces: CUT

function _aggregate_by_threshold_barrier(
  threshold, cell_to_unit_cut_meas, facet_to_inoutcut, cell_to_inoutcut,
  loc, cell_to_coords, cell_to_faces, face_to_cells)

  n_cells = length(cell_to_unit_cut_meas)
  cell_to_cellin = zeros(Int32, n_cells)
  cell_to_touched = fill(false, n_cells)

  for cell in 1:n_cells
    if cell_to_unit_cut_meas[cell] >= threshold
      cell_to_cellin[cell] = cell
      cell_to_touched[cell] = true
    end
  end

  c1 = array_cache(cell_to_faces)
  c2 = array_cache(face_to_cells)
  c3 = array_cache(cell_to_coords)
  c4 = array_cache(cell_to_coords)

  max_iters = 400  # upstream: 20 — raised for ultra-fine grids

  all_aggregated = false
  for iter in 1:max_iters
    all_aggregated = true
    for cell in 1:n_cells
      if !cell_to_touched[cell] && cell_to_inoutcut[cell] == CUT
        neigh_cell = _find_best_neighbor(
          c1, c2, c3, c4, cell,
          cell_to_faces,
          face_to_cells,
          cell_to_coords,
          cell_to_touched,
          cell_to_cellin,
          facet_to_inoutcut,
          loc)
        if neigh_cell > 0
          cellin = cell_to_cellin[neigh_cell]
          cell_to_cellin[cell] = cellin
        else
          all_aggregated = false
        end
      end
    end
    if all_aggregated
      break
    end
    _touch_aggregated_cells!(cell_to_touched, cell_to_cellin)
  end

  if !all_aggregated
    n_cut = count(i -> cell_to_inoutcut[i] == CUT, 1:n_cells)
    n_unagg = count(
      i -> cell_to_cellin[i] == 0 && cell_to_inoutcut[i] == CUT,
      1:n_cells)

    # Save diagnostic info about unaggregated cells
    unagg_cells = findall(
      i -> cell_to_cellin[i] == 0 && cell_to_inoutcut[i] == CUT,
      1:n_cells)

    error(
      "Aggregation failed: $n_unagg / $n_cut CUT cells remain unaggregated " *
      "after $max_iters iterations (n_cells=$n_cells). " *
      "Unaggregated cell IDs: $(join(unagg_cells[1:min(10, length(unagg_cells))], ", "))" *
      (length(unagg_cells) > 10 ? "..." : "") *
      ". Consider increasing max_iters in Aggregation/AggregationStrategy.jl " *
      "or inspecting the geometry near pinch-off regions."
    )
  end

  cell_to_cellin
end
