module Reinitialization

using LinearAlgebra
using Gridap.Geometry: get_node_coordinates, get_grid_topology, get_faces
using GridapEmbedded: cut  # also available from parent Geometric module scope
using GridapEmbedded.LevelSetCutters: DiscreteGeometry

using ..CartesianField: CartesianMeshField
using ..Geometric: EvolvingDiscreteGeometry, current_cut, get_active_indices, grid_info, CartesianGridInfo
export reinitialize!

# Geometric Subcell Helpers (Russo-Smereka with Geometric Clamping)

# Maximum allowed curvature relative to local level set magnitude
# C = 2.0 allows representation of the sharpest possible grid feature (sawtooth pattern)
const CURVATURE_CLAMP_FACTOR = 2.0

@inline function _linear_subcell_dist(ϕᵢ, ϕⱼ, h)
    h * abs(ϕᵢ) / abs(ϕᵢ - ϕⱼ)
end

function _quadratic_subcell_dist(phi::Vector{Float64}, info::CartesianGridInfo,
    idx::Int, dim::Int, direction::Int)
    nx, ny = info.dims
    hx, hy = info.spacing
    h = (dim == 1) ? hx : hy

    i = (idx - 1) % nx + 1
    j = (idx - 1) ÷ nx + 1

    # Get neighbor index with bounds checking
    neighbor_idx = if dim == 1
        idx + direction
    else
        idx + direction * nx
    end

    # Early exit if neighbor is out of bounds - use first-order approximation
    if neighbor_idx < 1 || neighbor_idx > nx * ny
        return abs(phi[idx])
    end

    # Check if neighbor is in valid row (x-direction) - ensures we're not wrapping
    if dim == 1
        ni = (neighbor_idx - 1) % nx + 1
        if ni < 1 || ni > nx
            return _linear_subcell_dist(phi[idx], phi[idx+direction], h)
        end
    end

    ϕᵢ = phi[idx]
    ϕⱼ = phi[neighbor_idx]

    # Compute second derivative using central difference (second-order accurate for SDFs)
    # Standard central difference: ϕ[i-1] - 2ϕ[i] + ϕ[i+1]
    ϕₓₓ = 0.0

    if dim == 1  # x-direction
        if i > 1 && i < nx
            im = idx - 1
            ip = idx + 1
            ϕₓₓ = phi[im] - 2ϕᵢ + phi[ip]
        end
    else  # y-direction
        if j > 1 && j < ny
            jm = idx - nx
            jp = idx + nx
            ϕₓₓ = phi[jm] - 2ϕᵢ + phi[jp]
        end
    end

    # Geometric Monotonicity Limiter: Clamp curvature based on local level set magnitude
    # This prevents the parabola from swinging wildly outside the data range
    # C = 2.0 allows representation of the sharpest possible grid feature (sawtooth pattern)
    S = abs(ϕᵢ) + abs(ϕⱼ)
    if S > 0
        ϕₓₓ = clamp(ϕₓₓ, -CURVATURE_CLAMP_FACTOR * S, CURVATURE_CLAMP_FACTOR * S)
    end

    # Choose between quadratic and linear based on curvature
    ε = 1e-10
    if abs(ϕₓₓ) > ε
        # Quadratic interpolation: solve ϕᵢ + s*(ϕⱼ-ϕᵢ) + s²*ϕₓₓ/2 = 0
        a = ϕₓₓ / 2
        b = ϕⱼ - ϕᵢ
        c = ϕᵢ

        discriminant = b^2 - 4 * a * c

        if discriminant >= 0
            sqrt_disc = sqrt(discriminant)
            s = (-b + sign(ϕᵢ) * sqrt_disc) / (2 * a)
            # Ensure s is in valid range (0, 1)
            if 0 < s < 1
                return abs(s) * h
            end
        end
    end

    # Fallback to linear (robust and guaranteed solution exists via IVT)
    return _linear_subcell_dist(ϕᵢ, ϕⱼ, h)
end

@inline function _has_sign_change(ϕᵢ, ϕⱼ)
    ϕᵢ * ϕⱼ < 0
end

function _compute_geometric_distance(phi::Vector{Float64}, info::CartesianGridInfo,
    idx::Int, i::Int, j::Int)
    nx, ny = info.dims
    hx, hy = info.spacing

    d_min = Inf
    ϕᵢ = phi[idx]

    # Check 4 cardinal neighbors for sign changes
    # Right
    if i < nx && _has_sign_change(ϕᵢ, phi[idx+1])
        d = _quadratic_subcell_dist(phi, info, idx, 1, 1)
        d_min = min(d_min, d)
    end

    # Left
    if i > 1 && _has_sign_change(ϕᵢ, phi[idx-1])
        d = _quadratic_subcell_dist(phi, info, idx, 1, -1)
        d_min = min(d_min, d)
    end

    # Up (positive y direction)
    if j < ny && _has_sign_change(ϕᵢ, phi[idx+nx])
        d = _quadratic_subcell_dist(phi, info, idx, 2, 1)
        d_min = min(d_min, d)
    end

    # Down (negative y direction)
    if j > 1 && _has_sign_change(ϕᵢ, phi[idx-nx])
        d = _quadratic_subcell_dist(phi, info, idx, 2, -1)
        d_min = min(d_min, d)
    end

    # If no sign changes found, fall back to absolute value
    return isfinite(d_min) ? d_min : abs(ϕᵢ)
end

"""
    reinitialize!(geom::EvolvingDiscreteGeometry)

Reinitialize the level set to a signed distance function using Fast Sweeping Method (FSM).
Uses Russo-Smereka geometric subcell anchoring for robust interface detection.
"""
function reinitialize!(geom::EvolvingDiscreteGeometry)
    phi_orig = geom.levelset
    info = grid_info(geom.grid)
    nx, ny = info.dims
    hx, hy = info.spacing

    # 1. Initialize distance field with large values
    t = fill(1e10, nx * ny)
    frozen = fill(false, nx * ny)

    # 2. Anchoring Phase (Seeding) - Russo-Smereka Geometric Subcell
    cut_nodes = _get_cut_nodes(geom)
    for idx in cut_nodes
        i = (idx - 1) % nx + 1
        j = (idx - 1) ÷ nx + 1

        # Compute geometric distance to interface (no gradients!)
        t[idx] = _compute_geometric_distance(phi_orig, info, idx, i, j)
        frozen[idx] = true
    end

    # 3. Fast Sweeping Phase
    t_mesh = CartesianMeshField(t, info)
    _fast_sweeping!(t_mesh, frozen)

    # 4. Finalize: Restore Sign
    @inbounds for i in eachindex(phi_orig)
        geom.levelset[i] = sign(phi_orig[i]) * t[i]
    end

    # Invalidate cache (including active_nodes which depends on the cut)
    geom.cache.cut = nothing
    geom.cache.active_nodes = nothing
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

    # 4 alternating-direction sweeps (standard FSM for 2D Eikonal).
    # Actual iteration order: (sy=+1,sx=+1), (sy=+1,sx=-1), (sy=-1,sx=+1), (sy=-1,sx=-1).
    for sy in (1, -1), sx in (1, -1)
        range_x = sx > 0 ? (1:nx) : (nx:-1:1)
        range_y = sy > 0 ? (1:ny) : (ny:-1:1)

        for j in range_y, i in range_x
            idx = i + (j - 1) * nx
            if frozen[idx]
                continue
            end

            # Neighbors in the "backwards" direction of the sweep
            im = i - sx
            jm = j - sy

            tx = (1 <= im <= nx) ? t_mesh.data[im+(j-1)*nx] : 1e10
            ty = (1 <= jm <= ny) ? t_mesh.data[i+(jm-1)*nx] : 1e10

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
    inv_hx2 = 1.0 / (hx * hx)
    inv_hy2 = 1.0 / (hy * hy)

    A = inv_hx2 + inv_hy2
    B = -2.0 * (tx * inv_hx2 + ty * inv_hy2)
    C = tx^2 * inv_hx2 + ty^2 * inv_hy2 - 1.0

    Δ = B^2 - 4 * A * C

    if Δ < 0
        return min(tx + hx, ty + hy)
    end

    sol = (-B + sqrt(Δ)) / (2 * A)

    if sol > tx && sol > ty
        return sol
    else
        return min(tx + hx, ty + hy)
    end
end

end # module
