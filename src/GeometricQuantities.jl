# =============================================================================
# GeometricQuantities — Narrow-Band Curvature/Normal for CutFEM
# =============================================================================
# Computes curvature and normals restricted to GridapEmbedded's active mesh,
# avoiding wasteful full-grid calculations.

using StaticArrays

# =============================================================================
# Cell Classification Utilities
# =============================================================================

"""
    get_cut_cells(cut_geo) -> Vector{Int}

Extract indices of cells that are cut by the interface.
Uses GridapEmbedded's internal cell classification.
"""
function get_cut_cells(cut_geo::GridapEmbedded.Interfaces.EmbeddedDiscretization)
    cell_class = cut_geo.ls_to_bgcell_to_inoutcut[1]
    return findall(==(GridapEmbedded.Interfaces.CUT), cell_class)
end

"""
    expand_to_band(cut_cells, n_layers, partition) -> Vector{Int}

Expand a set of cell indices by `n_layers` to form a narrow band.
Uses 4-neighbor connectivity for Cartesian grids.

# Arguments
- `cut_cells`: Initial cell indices (typically CUT cells)
- `n_layers`: Number of layers to expand (default: 2 for curvature stencil)
- `partition`: Grid cell counts (nx_cells, ny_cells)
"""
function expand_to_band(cut_cells::Vector{Int}, n_layers::Int, partition::NTuple{2,Int})
    nx, ny = partition
    
    # Convert linear index to (i,j)
    function cell_to_ij(cell)
        j = div(cell - 1, nx) + 1
        i = mod(cell - 1, nx) + 1
        return (i, j)
    end
    
    # Convert (i,j) to linear index
    function ij_to_cell(i, j)
        return i + (j - 1) * nx
    end
    
    band = Set(cut_cells)
    
    for _ in 1:n_layers
        new_cells = Set{Int}()
        for cell in band
            i, j = cell_to_ij(cell)
            # Add 4-neighbors
            for (di, dj) in ((-1, 0), (1, 0), (0, -1), (0, 1))
                ni, nj = i + di, j + dj
                if 1 <= ni <= nx && 1 <= nj <= ny
                    push!(new_cells, ij_to_cell(ni, nj))
                end
            end
        end
        union!(band, new_cells)
    end
    
    return collect(band)
end

"""
    cells_to_nodes(cells, partition) -> Set{Int}

Get all node indices belonging to a set of cells.
Assumes bilinear (4-node) quad cells on Cartesian grid.

# Arguments
- `cells`: Cell indices
- `partition`: Grid cell counts (nx_cells, ny_cells)
"""
function cells_to_nodes(cells::Vector{Int}, partition::NTuple{2,Int})
    nx_cells, ny_cells = partition
    nx_nodes = nx_cells + 1
    
    # Convert cell index to its 4 node indices
    function cell_nodes(cell)
        j = div(cell - 1, nx_cells) + 1
        i = mod(cell - 1, nx_cells) + 1
        # Node indices (column-major)
        n1 = i + (j - 1) * nx_nodes
        n2 = (i + 1) + (j - 1) * nx_nodes
        n3 = i + j * nx_nodes
        n4 = (i + 1) + j * nx_nodes
        return (n1, n2, n3, n4)
    end
    
    nodes = Set{Int}()
    for cell in cells
        for n in cell_nodes(cell)
            push!(nodes, n)
        end
    end
    return nodes
end

# =============================================================================
# Finite Difference Curvature Computation  
# =============================================================================

"""
    compute_curvature_fd(ϕ, band_nodes, info) -> Dict{Int,Float64}

Compute curvature at specified nodes using finite differences.
Only computes at nodes where the full stencil is available.

Uses the level set curvature formula:
    κ = (ϕ_xx * ϕ_y² - 2*ϕ_x*ϕ_y*ϕ_xy + ϕ_yy * ϕ_x²) / |∇ϕ|³
"""
function compute_curvature_fd(ϕ::AbstractVector{<:Real}, 
                               band_nodes::Set{Int}, 
                               info::CartesianGridInfo)
    nx, ny = info.dims
    Δx, Δy = info.spacing
    
    # Regularization to avoid division by zero
    ε = 1e-12
    
    κ = Dict{Int,Float64}()
    
    for node in band_nodes
        # Convert node to (i,j) - 1-indexed
        j = div(node - 1, nx) + 1
        i = mod(node - 1, nx) + 1
        
        # Skip boundary nodes where stencil unavailable
        if i < 2 || i > nx - 1 || j < 2 || j > ny - 1
            continue
        end
        
        # Helper to get ϕ at (i,j)
        idx(ii, jj) = ii + (jj - 1) * nx
        
        # First derivatives (central)
        ϕ_x = (ϕ[idx(i+1, j)] - ϕ[idx(i-1, j)]) / (2Δx)
        ϕ_y = (ϕ[idx(i, j+1)] - ϕ[idx(i, j-1)]) / (2Δy)
        
        # Second derivatives (central)
        ϕ_xx = (ϕ[idx(i+1, j)] - 2ϕ[idx(i, j)] + ϕ[idx(i-1, j)]) / (Δx^2)
        ϕ_yy = (ϕ[idx(i, j+1)] - 2ϕ[idx(i, j)] + ϕ[idx(i, j-1)]) / (Δy^2)
        ϕ_xy = (ϕ[idx(i+1, j+1)] - ϕ[idx(i-1, j+1)] - ϕ[idx(i+1, j-1)] + ϕ[idx(i-1, j-1)]) / (4Δx * Δy)
        
        # Gradient magnitude
        grad_mag = sqrt(ϕ_x^2 + ϕ_y^2)
        
        # Curvature formula
        if grad_mag > ε
            κ[node] = (ϕ_xx * ϕ_y^2 - 2ϕ_x * ϕ_y * ϕ_xy + ϕ_yy * ϕ_x^2) / (grad_mag^3)
        else
            κ[node] = 0.0
        end
    end
    
    return κ
end

"""
    compute_gradient_fd(ϕ, band_nodes, info) -> Dict{Int,SVector{2,Float64}}

Compute gradient ∇ϕ at specified nodes using central finite differences.
"""
function compute_gradient_fd(ϕ::AbstractVector{<:Real}, 
                              band_nodes::Set{Int}, 
                              info::CartesianGridInfo)
    nx, ny = info.dims
    Δx, Δy = info.spacing
    
    ∇ϕ = Dict{Int,SVector{2,Float64}}()
    
    for node in band_nodes
        j = div(node - 1, nx) + 1
        i = mod(node - 1, nx) + 1
        
        # Skip boundary nodes
        if i < 2 || i > nx - 1 || j < 2 || j > ny - 1
            continue
        end
        
        idx(ii, jj) = ii + (jj - 1) * nx
        
        ϕ_x = (ϕ[idx(i+1, j)] - ϕ[idx(i-1, j)]) / (2Δx)
        ϕ_y = (ϕ[idx(i, j+1)] - ϕ[idx(i, j-1)]) / (2Δy)
        
        ∇ϕ[node] = SVector(ϕ_x, ϕ_y)
    end
    
    return ∇ϕ
end

# =============================================================================
# High-Level API
# =============================================================================

# Default narrow band layers (2 provides buffer for curvature stencil)
const DEFAULT_BAND_LAYERS = 2

"""
    curvature_at_band(eg::EvolvingDiscreteGeometry; n_layers=2) -> (Vector{Float64}, Set{Int})

Compute curvature only in the narrow band around the interface.

# Returns
- `κ_full`: Full-length vector with curvature values (zeros outside band)
- `band_nodes`: Set of node indices where curvature was computed

# Example
```julia
κ, band = curvature_at_band(eg)
# κ is ready to create FE function: FEFunction(V, κ)
```
"""
function curvature_at_band(eg::EvolvingDiscreteGeometry; n_layers::Int=DEFAULT_BAND_LAYERS)
    # Get cut geometry (uses cached version)
    cut_geo = current_cut(eg)
    info = grid_info(eg)
    ϕ = current_levelset(eg)
    
    # Get cells in narrow band
    cut_cells = get_cut_cells(cut_geo)
    band_cells = expand_to_band(cut_cells, n_layers, info.cells)
    band_nodes = cells_to_nodes(band_cells, info.cells)
    
    # Compute curvature at band nodes
    κ_dict = compute_curvature_fd(ϕ, band_nodes, info)
    
    # Create full vector (zeros outside band)
    κ_full = zeros(Float64, prod(info.dims))
    for (node, val) in κ_dict
        κ_full[node] = val
    end
    
    return κ_full, band_nodes
end

"""
    interface_curvature(eg::EvolvingDiscreteGeometry) -> CellField

Return curvature as a CellField defined on the background mesh.
Can be evaluated on the interface Γ in variational forms.

Uses cached computation - only recomputes when geometry changes.

# Example
```julia
cut_geo = current_cut(eg)
Γ = EmbeddedBoundary(cut_geo)
n_Γ = get_normal_vector(Γ)
κ_Γ = interface_curvature(eg)

dΓ = Measure(Γ, 2)
σ = 0.072  # Surface tension
∫( σ * κ_Γ * (n_Γ ⋅ v) )dΓ
```
"""
function interface_curvature(eg::EvolvingDiscreteGeometry)
    # Check cache
    if !isnothing(eg.cached_curvature) && !eg.curvature_dirty
        return eg.cached_curvature
    end
    
    # Compute curvature in narrow band
    κ_values, _ = curvature_at_band(eg)
    
    # Create FE space on background mesh
    V = FESpace(eg.bg_model, ReferenceFE(lagrangian, Float64, 1); conformity=:H1)
    κ_h = FEFunction(V, κ_values)
    
    # Cache and return
    eg.cached_curvature = κ_h
    eg.curvature_dirty = false
    
    return κ_h
end

"""
    interface_normal(eg::EvolvingDiscreteGeometry) -> CellField

Return the normal vector on the embedded boundary Γ.
This is a thin wrapper around GridapEmbedded's `get_normal_vector`.

# Example
```julia
cut_geo = current_cut(eg)
Γ = EmbeddedBoundary(cut_geo)
n_Γ = interface_normal(eg)  # Equivalent to get_normal_vector(Γ)
```
"""
function interface_normal(eg::EvolvingDiscreteGeometry)
    cut_geo = current_cut(eg)
    Γ = EmbeddedBoundary(cut_geo)
    return get_normal_vector(Γ)
end

# =============================================================================
# Visualization Helpers (stubs - implemented in Makie extension)
# =============================================================================

"""
    plot_curvature(eg::EvolvingDiscreteGeometry; kwargs...)

Visualize curvature field with the interface contour overlay.
Requires GLMakie: `using GLMakie` before calling.
"""
function plot_curvature end

"""
    plot_curvature!(ax, eg::EvolvingDiscreteGeometry; kwargs...)

Add curvature visualization to an existing Makie axis.
"""
function plot_curvature! end
