module Reinitialization

using LinearAlgebra
using Gridap.Geometry: get_node_coordinates, get_grid_topology, get_faces
using GridapEmbedded: cut
using GridapEmbedded.LevelSetCutters: DiscreteGeometry

using ..CartesianField: CartesianMeshField
using ..Geometric: EvolvingDiscreteGeometry, current_cut, get_active_indices, grid_info, CartesianGridInfo
using ..Stencils: weno5⁻, weno5⁺

export reinitialize!

"""
    reinitialize!(geom::EvolvingDiscreteGeometry; eps_grad=1e-12)

Reinitialize the level set to a signed distance function using Fast Sweeping Method (FSM).
Uses WENO5 for high-precision interface anchoring.
"""
function reinitialize!(geom::EvolvingDiscreteGeometry; eps_grad=1e-12)
    phi_orig = geom.levelset
    info = grid_info(geom.grid)
    nx, ny = info.dims
    hx, hy = info.spacing
    
    # 1. Initialize distance field with large values
    t = fill(1e10, nx * ny)
    frozen = fill(false, nx * ny)
    
    # 2. Anchoring Phase (Seeding)
    phi_wrapper = CartesianMeshField(phi_orig, info)
    
    cut_nodes = _get_cut_nodes(geom)
    for idx in cut_nodes
        i = (idx - 1) % nx + 1
        j = (idx - 1) ÷ nx + 1
        I = CartesianIndex(i, j)
        
        # Compute |grad phi| using WENO5
        # We use a central-like average of weno stencils
        gx = 0.5 * (weno5⁻(phi_wrapper, I, 1) + weno5⁺(phi_wrapper, I, 1))
        gy = 0.5 * (weno5⁻(phi_wrapper, I, 2) + weno5⁺(phi_wrapper, I, 2))
        
        gnorm = sqrt(gx^2 + gy^2)
        
        # Dirichlet seed: d = |phi| / |grad phi|
        t[idx] = abs(phi_orig[idx]) / (gnorm + eps_grad)
        frozen[idx] = true
    end
    
    # 3. Fast Sweeping Phase
    t_mesh = CartesianMeshField(t, info)
    _fast_sweeping!(t_mesh, frozen)
    
    # 4. Finalize: Restore Sign
    @inbounds for i in eachindex(phi_orig)
        geom.levelset[i] = sign(phi_orig[i]) * t[i]
    end
    
    # Invalidate cache
    geom.cache.cut = nothing
    geom.cache.transfer_op = nothing
    geom.cache.extension_op = nothing
    
    return geom
end

function _get_cut_nodes(geom::EvolvingDiscreteGeometry)
    cut_geo = current_cut(geom)
    if isnothing(cut_geo)
        # Force cut if not available
        coords = vec(collect(get_node_coordinates(geom.grid)))
        geo = DiscreteGeometry(geom.levelset, coords)
        cut_geo = cut(geom.grid, geo)
        geom.cache.cut = cut_geo
    end
    
    raw_status = cut_geo.ls_to_bgcell_to_inoutcut
    status = if eltype(raw_status) <: AbstractVector
        raw_status[1]
    else
        raw_status
    end
    
    const_CUT = 0
    cut_cell_indices = findall(x -> x == const_CUT, status)
    
    topo = get_grid_topology(geom.grid)
    cell_to_nodes = get_faces(topo, 2, 0)
    
    cut_nodes = Int[]
    sizehint!(cut_nodes, length(cut_cell_indices) * 4)
    for cell_idx in cut_cell_indices
        nodes = cell_to_nodes[cell_idx]
        append!(cut_nodes, nodes)
    end
    unique!(cut_nodes)
    return cut_nodes
end

function _fast_sweeping!(t_mesh::CartesianMeshField, frozen)
    nx, ny = t_mesh.grid.dims
    hx, hy = t_mesh.grid.spacing
    
    # 4 Sweeps (2D Gray Code order)
    # (1,1), (1,-1), (-1,-1), (-1,1)
    for sy in (1, -1), sx in (1, -1)
        range_x = sx > 0 ? (1:nx) : (nx:-1:1)
        range_y = sy > 0 ? (1:ny) : (ny:-1:1)
        
        for j in range_y, i in range_x
            idx = i + (j-1)*nx
            if frozen[idx]
                continue
            end
            
            # Neighbors in the "backwards" direction of the sweep
            im = i - sx
            jm = j - sy
            
            tx = (1 <= im <= nx) ? t_mesh.data[im + (j-1)*nx] : 1e10
            ty = (1 <= jm <= ny) ? t_mesh.data[i + (jm-1)*nx] : 1e10
            
            # Solve Eikonal equation locally
            t_new = _solve_eikonal_local(tx, ty, hx, hy)
            
            if t_new < t_mesh.data[idx]
                t_mesh.data[idx] = t_new
            end
        end
    end
end

@inline function _solve_eikonal_local(tx, ty, hx, hy)
    # Quadratic: (t-tx)^2/hx^2 + (t-ty)^2/hy^2 = 1
    inv_hx2 = 1.0 / (hx*hx)
    inv_hy2 = 1.0 / (hy*hy)
    
    A = inv_hx2 + inv_hy2
    B = -2.0 * (tx * inv_hx2 + ty * inv_hy2)
    C = tx^2 * inv_hx2 + ty^2 * inv_hy2 - 1.0
    
    Δ = B^2 - 4*A*C
    
    if Δ < 0
        return min(tx + hx, ty + hy)
    end
    
    sol = (-B + sqrt(Δ)) / (2*A)
    
    if sol > tx && sol > ty
        return sol
    else
        return min(tx + hx, ty + hy)
    end
end

end # module
