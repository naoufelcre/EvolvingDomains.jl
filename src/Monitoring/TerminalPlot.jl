module TerminalPlot

using ..Geometric: EvolvingDiscreteGeometry, grid_info
import TPlot: plot, plot_geometry, plot_curves, Geometry

export plot, plot_geometry, plot_curves

_terminal_phi(geom::EvolvingDiscreteGeometry) =
    reshape(geom.levelset, grid_info(geom.grid).dims)

"""
    plot_geometry(geom::EvolvingDiscreteGeometry; kwargs...)

Render the level-set interior with TPlot. Existing terminal plotting keywords apply.
"""
plot_geometry(geom::EvolvingDiscreteGeometry; kwargs...) =
    plot_geometry(_terminal_phi(geom); kwargs...)

plot(geom::EvolvingDiscreteGeometry; kwargs...) = plot_geometry(geom; kwargs...)

plot(geom::EvolvingDiscreteGeometry, t::AbstractVector,
    y::Union{AbstractVector,AbstractMatrix}; kwargs...) =
    plot(_terminal_phi(geom), t, y; kwargs...)

# The panel holds a view of the live level-set values, not a geometry copy.
Geometry(geom::EvolvingDiscreteGeometry; kwargs...) =
    Geometry(_terminal_phi(geom); kwargs...)

end # module
