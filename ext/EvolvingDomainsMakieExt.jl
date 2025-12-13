# =============================================================================
# EvolvingDomainsMakieExt — GLMakie Visualization Extension
# =============================================================================
# This extension provides interactive visualization for EvolvingDomains.jl
# It only loads when the user imports both EvolvingDomains and GLMakie.
#
# Usage:
#   using EvolvingDomains
#   using GLMakie
#   
#   # Quick contour plot
#   plot_levelset(eg)
#   
#   # Interactive window with slider
#   plot_levelset(eg; interactive=true)

module EvolvingDomainsMakieExt

using EvolvingDomains
using EvolvingDomains: CartesianGridInfo, EvolvingDiscreteGeometry
using EvolvingDomains: current_levelset, current_time, grid_info, advance!, reinitialize!

using GLMakie
using GLMakie: Figure, Axis, Colorbar, Observable, lift, record
using GLMakie: contour!, contourf!, heatmap!, scatter!, lines!
using GLMakie: linkaxes!, hidedecorations!, DataAspect

# =============================================================================
# Grid Coordinate Helpers
# =============================================================================

"""
    _grid_coords(info::CartesianGridInfo{2})

Generate x, y coordinate ranges from grid info.
"""
function _grid_coords(info::CartesianGridInfo{2})
    nx, ny = info.dims
    x = range(info.origin[1], step=info.spacing[1], length=nx)
    y = range(info.origin[2], step=info.spacing[2], length=ny)
    return x, y
end

"""
    _reshape_levelset(ϕ::Vector, info::CartesianGridInfo{2})

Reshape level set vector to 2D grid (column-major).
"""
function _reshape_levelset(ϕ::Vector, info::CartesianGridInfo{2})
    nx, ny = info.dims
    return reshape(ϕ, (nx, ny))
end

# =============================================================================
# Core Visualization Functions
# =============================================================================

"""
    EvolvingDomains.plot_levelset(eg::EvolvingDiscreteGeometry; kwargs...)

Create an interactive contour plot of the current level set.

# Keyword Arguments
- `colormap::Symbol = :RdBu` : Colormap for the level set values
- `show_zero::Bool = true` : Highlight the zero contour (domain boundary)
- `title::String = ""` : Plot title (auto-generated if empty)
- `filled::Bool = true` : Use filled contours for interior

# Returns
- `Figure` : The Makie figure object

# Example
```julia
using EvolvingDomains, GLMakie
# ... create eg ...
fig = plot_levelset(eg)
fig = plot_levelset(eg; colormap=:viridis, title="My Geometry")
```
"""
function EvolvingDomains.plot_levelset(eg::EvolvingDiscreteGeometry; 
                                        colormap::Symbol = :RdBu,
                                        show_zero::Bool = true,
                                        title::String = "",
                                        filled::Bool = true)
    info = grid_info(eg)
    ϕ = current_levelset(eg)
    t = current_time(eg)
    
    x, y = _grid_coords(info)
    ϕ_2d = _reshape_levelset(ϕ, info)
    
    # Auto title
    if isempty(title)
        title = "Level Set (t = $(round(t, digits=4)))"
    end
    
    # Create figure
    fig = Figure(size = (700, 600))
    ax = Axis(fig[1, 1], 
              title = title,
              xlabel = "x",
              ylabel = "y",
              aspect = DataAspect())
    
    # Plot level set
    if filled
        # Filled contours showing inside/outside
        hm = heatmap!(ax, x, y, ϕ_2d; colormap = colormap, colorrange = (-1, 1) .* maximum(abs, ϕ))
        Colorbar(fig[1, 2], hm, label = "ϕ")
    end
    
    # Zero contour (domain boundary) - always prominent
    if show_zero
        contour!(ax, x, y, ϕ_2d; levels = [0.0], linewidth = 3, color = :black)
    end
    
    # Additional contours for structure
    contour!(ax, x, y, ϕ_2d; levels = 10, linewidth = 0.5, color = :gray50, linestyle = :dash)
    
    return fig
end

"""
    EvolvingDomains.plot_levelset!(ax, eg::EvolvingDiscreteGeometry; kwargs...)

Add level set visualization to an existing axis.

# Example
```julia
fig = Figure()
ax = Axis(fig[1, 1])
plot_levelset!(ax, eg)
```
"""
function EvolvingDomains.plot_levelset!(ax, eg::EvolvingDiscreteGeometry;
                                         colormap::Symbol = :RdBu,
                                         show_zero::Bool = true,
                                         filled::Bool = true)
    info = grid_info(eg)
    ϕ = current_levelset(eg)
    
    x, y = _grid_coords(info)
    ϕ_2d = _reshape_levelset(ϕ, info)
    
    if filled
        heatmap!(ax, x, y, ϕ_2d; colormap = colormap, colorrange = (-1, 1) .* maximum(abs, ϕ))
    end
    
    if show_zero
        contour!(ax, x, y, ϕ_2d; levels = [0.0], linewidth = 3, color = :black)
    end
    
    contour!(ax, x, y, ϕ_2d; levels = 10, linewidth = 0.5, color = :gray50, linestyle = :dash)
    
    return ax
end

# Print message when extension loads
function __init__()
    @info "EvolvingDomainsMakieExt loaded: plot_levelset() available"
end

end # module

