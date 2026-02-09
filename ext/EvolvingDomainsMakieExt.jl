# =============================================================================
# EvolvingDomainsMakieExt — CairoMakie Visualization Extension
# =============================================================================

module EvolvingDomainsMakieExt

using EvolvingDomains
using CairoMakie

"""
    plot_levelset!(ax, geom::EvolvingDiscreteGeometry;
                   colormap=:RdBu, levels=20, show_zero=true, linecolor=:black,
                   ϕ_buffer=nothing)

Add level set visualization to an existing Makie axis.

# Arguments
- `ax`: Makie axis
- `geom`: EvolvingDiscreteGeometry to visualize
- `colormap`: Colormap for level set values (default: `:RdBu`)
- `levels`: Number of contour levels (default: 20)
- `show_zero`: Show the zero contour (interface) in bold (default: true)
- `linecolor`: Color for the zero contour (default: `:black`)
- `ϕ_buffer`: Optional pre-allocated buffer for level set values (for animation loops)

# Animation Loop Example
```julia
info = EvolvingDomains.grid_info(geom.model)
ϕ_buffer = zeros(prod(info.dims))  # Allocate once

record(fig, "output.gif", 1:n_steps) do step
    advance!(geom, dt)
    empty!(ax)
    plot_levelset!(ax, geom; ϕ_buffer=ϕ_buffer)  # Reuses buffer
end
```
"""
function EvolvingDomains.plot_levelset!(ax, geom::EvolvingDiscreteGeometry;
    colormap=:RdBu,
    levels::Int=20,
    show_zero::Bool=true,
    linewidth::Real=2,
    linecolor=:black,
    ϕ_buffer::Union{Nothing,Vector{Float64}}=nothing)
    info = EvolvingDomains.grid_info(geom.grid)
    nx, ny = info.dims

    # Use provided buffer or allocate new
    if ϕ_buffer !== nothing
        copyto!(ϕ_buffer, EvolvingDomains.current_levelset(geom))
        ϕ = ϕ_buffer
    else
        ϕ = EvolvingDomains.current_levelset(geom)
    end

    ϕ_2d = reshape(ϕ, nx, ny)

    # Build coordinate arrays
    ox, oy = info.origin
    dx, dy = info.spacing
    xs = range(ox, step=dx, length=nx)
    ys = range(oy, step=dy, length=ny)

    # Heatmap of level set values
    heatmap!(ax, xs, ys, ϕ_2d; colormap=colormap)

    # Zero contour (interface)
    if show_zero
        contour!(ax, xs, ys, ϕ_2d; levels=[0.0], color=linecolor, linewidth=linewidth)
    end

    return ax
end

"""
    plot_levelset(geom::EvolvingDiscreteGeometry; kwargs...)

Create a new figure with level set visualization.
"""
function EvolvingDomains.plot_levelset(geom::EvolvingDiscreteGeometry;
    figsize=(600, 600), kwargs...)
    fig = Figure(size=figsize)
    ax = Axis(fig[1, 1], aspect=1, title="Level Set")
    EvolvingDomains.plot_levelset!(ax, geom; kwargs...)
    return fig
end

end # module
