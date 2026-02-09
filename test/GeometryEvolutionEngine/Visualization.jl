module Visualization

using CairoMakie
using Gridap
using EvolvingDomains
using EvolvingDomains.Geometric: grid_info

using EvolvingDomains.Kinematic: RayMaps
using EvolvingDomains.Geometric: get_active_indices
using Gridap.Geometry: get_node_coordinates
using Gridap.TensorValues: VectorValue

export plot_cell_classification!, plot_raymaps!, get_raymap_bounding_box

function plot_cell_classification!(ax, geom::EvolvingDiscreteGeometry; state::Symbol=:current)
    info = grid_info(geom.grid)
    nx, ny = info.dims
    ox, oy = info.origin
    dx, dy = info.spacing

    # --- 1. PREPARE COORDINATES ---
    # Cell Centers (for Heatmap) - Size: (nx-1) x (ny-1)
    xs_cells = range(ox + dx / 2, step=dx, length=nx - 1)
    ys_cells = range(oy + dy / 2, step=dy, length=ny - 1)
    
    # Node Coordinates (for Contour & Grid Lines) - Size: nx x ny
    x_nodes = range(ox, step=dx, length=nx)
    y_nodes = range(oy, step=dy, length=ny)

    # --- 2. GET STATUS GRID ---
    status_grid = zeros(Int8, nx - 1, ny - 1)
    has_phi = false
    
    if state == :current
        # Compute from current levelset
        phi = geom.levelset
        has_phi = true
        
        for j in 1:ny-1
            for i in 1:nx-1
                # Indices of the 4 corners
                n1 = i + (j - 1) * nx
                n2 = (i + 1) + (j - 1) * nx
                n3 = i + j * nx
                n4 = (i + 1) + j * nx
    
                vals = (phi[n1], phi[n2], phi[n3], phi[n4])
                n_pos = count(v -> v > 0, vals)
                n_neg = count(v -> v < 0, vals)
    
                if n_pos == 4
                    status_grid[i, j] = 0 # OUT
                elseif n_neg == 4
                    status_grid[i, j] = 1 # IN
                else
                    status_grid[i, j] = 2 # CUT
                end
            end
        end
        
    elseif state == :prev
        # Retrieve from cache
        prev_cut = geom.cache.prev_cut
        
        if isnothing(prev_cut)
            @warn "No previous geometry state available (prev_cut is nothing). Ensure get_active_indices(geom) was called BEFORE the last advance! step."
            return ax
        end
        
        # Raw status: -1 (IN), 0 (CUT), 1 (OUT)
        # We MUST collect to ensure it's a standard array before reshape
        raw_status = prev_cut.ls_to_bgcell_to_inoutcut
        flat_status = if eltype(raw_status) <: AbstractVector
            collect(raw_status[1])
        else
            collect(raw_status)
        end
        
        # Reshape to grid (nx-1, ny-1)
        try
            reshaped_raw = reshape(flat_status, nx-1, ny-1)
             for j in 1:ny-1
                for i in 1:nx-1
                    val = reshaped_raw[i, j]
                    if val == -1 # IN
                        status_grid[i, j] = 1
                    elseif val == 0 # CUT
                        status_grid[i, j] = 2
                    else # OUT (1)
                        status_grid[i, j] = 0
                    end
                end
            end
        catch e
             @warn "Failed to reshape previous status: $e. Flat size: $(length(flat_status)), expected: $((nx-1)*(ny-1))"
             return ax
        end
    end

    # --- 3. PLOT HEATMAP ---
    cmap = cgrad([:white, :lightblue, :gray50], 3, categorical=true)
    heatmap!(ax, xs_cells, ys_cells, status_grid, colormap=cmap, colorrange=(0, 2))

    # --- 4. PLOT CONTOUR (Only if phi available) ---
    if has_phi
        phi_matrix = reshape(geom.levelset, nx, ny)
        contour!(ax, x_nodes, y_nodes, phi_matrix, levels=[0.0], color=:black, linewidth=1.5)
    end

    # --- 5. PLOT GRID LINES ---
    vlines!(ax, x_nodes, color=:black, linewidth=0.5, alpha=0.3)
    hlines!(ax, y_nodes, color=:black, linewidth=0.5, alpha=0.3)

    return ax
end

function get_raymap_bounding_box(geom::EvolvingDiscreteGeometry, maps::RayMaps; 
                                 stride=1, type=:backward, margin=0.1)
    # 1. Select Map and Indices based on Type
    if type == :backward
        active_indices = get_active_indices(geom, :current)
        target_points = maps.backward
    elseif type == :forward
        active_indices = get_active_indices(geom, :prev)
        target_points = maps.forward
    else
        error("Unknown map type: $type. Use :backward or :forward.")
    end

    if isempty(active_indices)
        return (0.0, 1.0, 0.0, 1.0) # Default fallback
    end

    coords = get_node_coordinates(geom.grid) 
    indices = 1:stride:length(active_indices)
    
    min_x, max_x = Inf, -Inf
    min_y, max_y = Inf, -Inf
    
    found_points = false

    for i in indices
        found_points = true
        node_idx = active_indices[i]
        pt_start = coords[node_idx]
        pt_end = target_points[i]
        
        min_x = min(min_x, pt_start[1], pt_end[1])
        max_x = max(max_x, pt_start[1], pt_end[1])
        min_y = min(min_y, pt_start[2], pt_end[2])
        max_y = max(max_y, pt_start[2], pt_end[2])
    end
    
    if !found_points
         return (0.0, 1.0, 0.0, 1.0)
    end

    dx = max_x - min_x
    dy = max_y - min_y
    
    mx = dx * margin
    my = dy * margin
    
    if dx ≈ 0 dx = 0.1 end
    if dy ≈ 0 dy = 0.1 end

    return (min_x - mx, max_x + mx, min_y - my, max_y + my)
end

function plot_raymaps!(ax, geom::EvolvingDiscreteGeometry, maps::RayMaps; 
                       stride=1, color=:blue, arrowsize=8, linewidth=1, 
                       type=:backward, fit_to_data=false, margin=0.1)
    
    # 1. Select Map and Indices based on Type
    if type == :backward
        # Backward Map: Defined on Current Grid (t_{n+1}) looking back to t_n.
        # We visualize flow t_n -> t_{n+1}.
        # Tail: Departure Point (t_n) -> maps.backward
        # Head: Grid Node (t_{n+1})     -> coords
        active_indices = get_active_indices(geom, :current)
        points_depart = maps.backward
        points_arrive = nothing # Will get from coords
    elseif type == :forward
        # Forward Map: Defined on Previous Grid (t_n) looking forward to t_{n+1}.
        # We visualize flow t_n -> t_{n+1}.
        # Tail: Grid Node (t_n)         -> coords
        # Head: Arrival Point (t_{n+1}) -> maps.forward
        active_indices = get_active_indices(geom, :prev)
        points_depart = nothing # Will get from coords
        points_arrive = maps.forward
    else
        error("Unknown map type: $type. Use :backward or :forward.")
    end
    
    # Check consistency
    mapped_points = isnothing(points_depart) ? points_arrive : points_depart
    if length(active_indices) != length(mapped_points)
        @warn "Map size mismatch: $(length(mapped_points)) points vs $(length(active_indices)) active nodes."
        return ax
    end
    
    if isempty(active_indices)
        return ax
    end

    coords = get_node_coordinates(geom.grid) 
    
    # Subsample
    indices = 1:stride:length(active_indices)
    
    xs = Float64[]
    ys = Float64[]
    us = Float64[]
    vs = Float64[]
    
    sizehint!(xs, length(indices))
    
    # For bounding box calculation
    min_x, max_x = Inf, -Inf
    min_y, max_y = Inf, -Inf

    for i in indices
        node_idx = active_indices[i]
        
        # Determine Start (Tail) and End (Head) of the arrow
        if type == :backward
            pt_start = coords[node_idx]       # Tail is on the Grid (Current)
            pt_end   = points_depart[i]       # Head is the Departure Point
        else # type == :forward
            pt_start = coords[node_idx]       # Tail is on the Grid (Previous)
            pt_end   = points_arrive[i]       # Head is the Arrival Point
        end
        
        # Vector: Always Start -> End (Flow Direction)
        vec = pt_end - pt_start
        
        push!(xs, pt_start[1])
        push!(ys, pt_start[2])
        push!(us, vec[1])
        push!(vs, vec[2])

        if fit_to_data
            min_x = min(min_x, pt_start[1], pt_end[1])
            max_x = max(max_x, pt_start[1], pt_end[1])
            min_y = min(min_y, pt_start[2], pt_end[2])
            max_y = max(max_y, pt_start[2], pt_end[2])
        end
    end
    
    arrows!(ax, xs, ys, us, vs; 
        arrowsize=arrowsize, 
        linewidth=linewidth,
        color=color,
        lengthscale=1.0 # Use actual length
    )
    
    if fit_to_data && !isempty(xs)
        dx = max_x - min_x
        dy = max_y - min_y
        
        # Add margin
        mx = dx * margin
        my = dy * margin
        
        # Ensure non-zero size if single point or flat
        if dx ≈ 0 dx = 0.1 end
        if dy ≈ 0 dy = 0.1 end

        limits!(ax, min_x - mx, max_x + mx, min_y - my, max_y + my)
    end

    return ax
end

end # module
