module EvolvingDomainsMakieExt

using EvolvingDomains
using EvolvingDomains: CartesianGridInfo, EvolvingDiscreteGeometry
using EvolvingDomains: current_levelset, current_time, grid_info, advance!, reinitialize!
using EvolvingDomains: SimulationFrame, SimulationResult

using GLMakie
using GLMakie: Figure, Axis, Colorbar, Observable, lift, record, Slider, Label
using GLMakie: contour!, contourf!, heatmap!, scatter!, lines!, on
using GLMakie: linkaxes!, hidedecorations!, DataAspect
using GLMakie: colsize!, Fixed, Relative

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

# =============================================================================
# Simulation Viewer (Interactive Slider)
# =============================================================================

"""
    EvolvingDomains.viewer(result::SimulationResult; kwargs...) -> Figure

Open an interactive viewer with a time slider for cached simulation frames.

# Arguments
- `result::SimulationResult`: Cached simulation frames from `snapshot()` calls
- `colormap::Symbol = :RdBu`: Colormap for level set visualization
- `wait::Bool = false`: If true, block until Enter is pressed

# Returns
- `Figure`: The Makie figure (for further customization)

# Example
```julia
frames = SimulationFrame[]
for step in 1:100
    advance!(eg, Δt)
    push!(frames, snapshot(eg))
end
result = SimulationResult(grid_info(eg), frames)
viewer(result)
```
"""
function EvolvingDomains.viewer(result::SimulationResult;
                                 colormap::Symbol = :RdBu,
                                 wait::Bool = false)
    info = result.grid_info
    frames = result.frames
    nframes = length(frames)
    
    if nframes == 0
        error("SimulationResult contains no frames")
    end
    
    x, y = _grid_coords(info)
    nx, ny = info.dims
    
    # Compute color range from all frames
    ϕ_max = maximum(maximum(abs, f.ϕ) for f in frames)
    colorrange = (-ϕ_max, ϕ_max)
    
    # Create figure
    fig = Figure(size = (700, 650))
    
    # Frame index observable
    frame_idx = Observable(1)
    
    # Title based on current frame
    title_obs = lift(frame_idx) do idx
        t = frames[idx].t
        "t = $(round(t, digits=4)) (frame $idx/$nframes)"
    end
    
    # Axis
    ax = Axis(fig[1, 1],
              title = title_obs,
              xlabel = "x", ylabel = "y",
              aspect = DataAspect(),
              limits = (info.origin[1], info.origin[1] + info.spacing[1]*(nx-1),
                        info.origin[2], info.origin[2] + info.spacing[2]*(ny-1)))
    
    # Slider
    sl = Slider(fig[2, 1], range = 1:nframes, startvalue = 1)
    connect!(frame_idx, sl.value)
    
    # Helper to render a frame
    function render_frame!(ax, frame)
        empty!(ax)
        heatmap!(ax, x, y, frame.ϕ; colormap = colormap, colorrange = colorrange)
        contour!(ax, x, y, frame.ϕ; levels = [0.0], linewidth = 3, color = :black)
    end
    
    # Initial render
    render_frame!(ax, frames[1])
    
    # React to slider changes
    on(frame_idx) do idx
        render_frame!(ax, frames[idx])
    end
    
    # Colorbar
    Colorbar(fig[1, 2], colormap = colormap, colorrange = colorrange,
             label = "ϕ", width = 12, height = Relative(0.7))
    colsize!(fig.layout, 2, Fixed(30))
    
    # Info label
    Label(fig[3, 1], "Drag slider to navigate through time",
          fontsize = 12, color = :gray50)
    
    display(fig)
    
    if wait
        println("Press Enter to close...")
        readline()
    end
    
    return fig
end

# =============================================================================
# Animated Playback
# =============================================================================

"""
    EvolvingDomains.view_live!(result::SimulationResult; kwargs...)

Play back cached simulation frames as an animation.

# Arguments
- `result::SimulationResult`: Cached simulation frames
- `fps::Real = 30`: Frames per second for playback
- `loop::Bool = false`: If true, loop animation continuously
- `colormap::Symbol = :RdBu`: Colormap for visualization

This function blocks until playback completes (or Ctrl+C to stop).

# Example
```julia
result = SimulationResult(grid_info(eg), frames)
view_live!(result; fps=15)
```
"""
function EvolvingDomains.view_live!(result::SimulationResult;
                                     fps::Real = 30,
                                     loop::Bool = false,
                                     colormap::Symbol = :RdBu)
    info = result.grid_info
    frames = result.frames
    nframes = length(frames)
    
    if nframes == 0
        error("SimulationResult contains no frames")
    end
    
    x, y = _grid_coords(info)
    nx, ny = info.dims
    
    # Compute color range from all frames
    ϕ_max = maximum(maximum(abs, f.ϕ) for f in frames)
    colorrange = (-ϕ_max, ϕ_max)
    
    # Create figure
    fig = Figure(size = (700, 600))
    ax = Axis(fig[1, 1],
              xlabel = "x", ylabel = "y",
              aspect = DataAspect(),
              limits = (info.origin[1], info.origin[1] + info.spacing[1]*(nx-1),
                        info.origin[2], info.origin[2] + info.spacing[2]*(ny-1)))
    
    # Initial render
    frame = frames[1]
    heatmap!(ax, x, y, frame.ϕ; colormap = colormap, colorrange = colorrange)
    contour!(ax, x, y, frame.ϕ; levels = [0.0], linewidth = 3, color = :black)
    
    Colorbar(fig[1, 2], colormap = colormap, colorrange = colorrange,
             label = "ϕ", width = 12, height = Relative(0.7))
    colsize!(fig.layout, 2, Fixed(30))
    
    display(fig)
    
    sleep_time = 1.0 / fps
    
    try
        while true
            for (idx, frame) in enumerate(frames)
                empty!(ax)
                ax.title = "t = $(round(frame.t, digits=4)) (frame $idx/$nframes)"
                heatmap!(ax, x, y, frame.ϕ; colormap = colormap, colorrange = colorrange)
                contour!(ax, x, y, frame.ϕ; levels = [0.0], linewidth = 3, color = :black)
                sleep(sleep_time)
            end
            
            if !loop
                break
            end
        end
    catch e
        if e isa InterruptException
            println("\nPlayback stopped")
        else
            rethrow(e)
        end
    end
    
    return nothing
end

# Print message when extension loads
function __init__()
    @info "EvolvingDomainsMakieExt loaded: plot_levelset(), viewer(), view_live!() available"
end

end # module
