"""
    filter_small_phase_islands!(levelset, info; phase=:negative, min_component_size=5,
                                connectivity=:edge, eps=nothing)

Remove small connected components from the selected phase by flipping level-set signs.

- `phase=:negative` targets nodes where `phi <= 0`.
- `phase=:positive` targets nodes where `phi >= 0`.
- `min_component_size` is measured in node count.
- `connectivity=:edge` uses 4-neighbor connectivity on the Cartesian node graph.

Returns a named tuple with component and flip statistics.
"""
function filter_small_phase_islands!(levelset::AbstractVector{<:Real}, info;
    phase::Symbol=:negative,
    min_component_size::Int=5,
    connectivity::Symbol=:edge,
    eps=nothing,
)
    if connectivity != :edge
        error("Unsupported connectivity=$connectivity. Use :edge for mesh-edge connectivity.")
    end
    if min_component_size <= 0
        return (n_components=0, n_removed_components=0, n_flipped_nodes=0)
    end

    nx, ny = info.dims
    n_nodes = nx * ny
    if length(levelset) != n_nodes
        error("Level-set size mismatch: got $(length(levelset)), expected $n_nodes from grid dims ($nx, $ny).")
    end

    in_phase = BitVector(undef, n_nodes)
    if phase == :negative
        @inbounds for idx in 1:n_nodes
            in_phase[idx] = levelset[idx] <= 0.0
        end
    elseif phase == :positive
        @inbounds for idx in 1:n_nodes
            in_phase[idx] = levelset[idx] >= 0.0
        end
    else
        error("Unsupported phase=$phase. Use :negative or :positive.")
    end

    visited = falses(n_nodes)
    stack = Int[]
    component_nodes = Int[]

    n_components = 0
    n_removed_components = 0
    n_flipped_nodes = 0

    @inbounds for start_idx in 1:n_nodes
        if !in_phase[start_idx] || visited[start_idx]
            continue
        end

        n_components += 1
        empty!(stack)
        empty!(component_nodes)
        push!(stack, start_idx)
        visited[start_idx] = true

        while !isempty(stack)
            idx = pop!(stack)
            push!(component_nodes, idx)

            i = ((idx - 1) % nx) + 1
            j = ((idx - 1) ÷ nx) + 1

            if i > 1
                nidx = idx - 1
                if in_phase[nidx] && !visited[nidx]
                    visited[nidx] = true
                    push!(stack, nidx)
                end
            end
            if i < nx
                nidx = idx + 1
                if in_phase[nidx] && !visited[nidx]
                    visited[nidx] = true
                    push!(stack, nidx)
                end
            end
            if j > 1
                nidx = idx - nx
                if in_phase[nidx] && !visited[nidx]
                    visited[nidx] = true
                    push!(stack, nidx)
                end
            end
            if j < ny
                nidx = idx + nx
                if in_phase[nidx] && !visited[nidx]
                    visited[nidx] = true
                    push!(stack, nidx)
                end
            end
        end

        if length(component_nodes) < min_component_size
            n_removed_components += 1
            local_eps = isnothing(eps) ? min(info.spacing...) : eps

            if phase == :negative
                for idx in component_nodes
                    levelset[idx] = local_eps
                end
            else
                for idx in component_nodes
                    levelset[idx] = -abs(local_eps)
                end
            end

            n_flipped_nodes += length(component_nodes)
        end
    end

    return (
        n_components=n_components,
        n_removed_components=n_removed_components,
        n_flipped_nodes=n_flipped_nodes,
    )
end
