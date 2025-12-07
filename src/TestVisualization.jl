# =============================================================================
# TestVisualization — Thin Wrapper for EmbeddedViz
# =============================================================================
# Provides EvolvingDomains-specific convenience functions.
# Delegates all actual visualization to EmbeddedViz.jl

module TestVisualization

using Gridap

export show_geometry, setup_test_viz!, close_test_viz
export start_recording!, save_animation

# Try to load EmbeddedViz
const EMBEDDED_VIZ = Ref{Any}(nothing)
const VIZ_ENABLED = Ref{Bool}(false)
const RECORDER = Ref{Any}(nothing)

"""
    setup_test_viz!(; enabled=true)

Initialize visualization using EmbeddedViz.
"""
function setup_test_viz!(; enabled::Bool=true)
    VIZ_ENABLED[] = enabled
    if enabled && isnothing(EMBEDDED_VIZ[])
        try
            @eval using EmbeddedViz
            eviz = @eval EmbeddedViz
            eviz.setup_viz!()
            EMBEDDED_VIZ[] = eviz
            println("[TestViz] Ready (using EmbeddedViz)")
            return true
        catch e
            @warn "[TestViz] EmbeddedViz not available, falling back to basic mode" exception=e
            # Fallback: try loading Plots directly
            try
                @eval using Plots
                EMBEDDED_VIZ[] = :plots_fallback
                println("[TestViz] Ready (fallback mode with Plots)")
                return true
            catch
                VIZ_ENABLED[] = false
                return false
            end
        end
    end
    return VIZ_ENABLED[]
end

"""
    start_recording!()

Start recording frames for animation.
"""
function start_recording!()
    !VIZ_ENABLED[] && return nothing
    eviz = EMBEDDED_VIZ[]
    isnothing(eviz) && return nothing
    
    if eviz !== :plots_fallback
        RECORDER[] = eviz.Recorder()
    end
    println("[TestViz] Recording started")
    return nothing
end

"""
    show_geometry(eg, label="")

Visualize current geometry state from EvolvingDiscreteGeometry.
"""
function show_geometry(eg, label::String="")
    !VIZ_ENABLED[] && return nothing
    eviz = EMBEDDED_VIZ[]
    isnothing(eviz) && return nothing
    
    # Get data from EvolvingDomains
    parent = parentmodule(TestVisualization)
    ϕ = parent.current_levelset(eg)
    t = parent.current_time(eg)
    model = eg.bg_model
    
    if eviz !== :plots_fallback
        # Use EmbeddedViz
        grid = eviz.GridInfo(model)
        rec = RECORDER[]
        
        if !isnothing(rec)
            eviz.record!(rec, :levelset, ϕ, grid; t=t, label=label)
        else
            eviz.show_levelset(ϕ, grid; t=t, label=label)
        end
    else
        # Fallback: basic contour plot
        _show_geometry_fallback(ϕ, model, t, label)
    end
    
    return nothing
end

"""
    save_animation(filename="geometry_evolution.gif"; fps=5)

Save recorded animation as GIF.
"""
function save_animation(filename::String="geometry_evolution.gif"; fps::Int=5)
    eviz = EMBEDDED_VIZ[]
    rec = RECORDER[]
    (isnothing(eviz) || isnothing(rec)) && return nothing
    
    if eviz !== :plots_fallback
        eviz.save_gif(rec, :levelset, filename; fps=fps)
    end
    
    RECORDER[] = nothing
    return filename
end

"""Clean up."""
function close_test_viz()
    VIZ_ENABLED[] = false
    EMBEDDED_VIZ[] = nothing
    RECORDER[] = nothing
    println("[TestViz] Closed")
end

# =============================================================================
# Fallback (when EmbeddedViz not available)
# =============================================================================

function _show_geometry_fallback(ϕ, model, t, label)
    Plt = @eval Plots
    desc = Gridap.Geometry.get_cartesian_descriptor(model)
    origin = desc.origin
    sizes = desc.sizes
    partition = desc.partition
    
    nx, ny = partition
    x = range(origin[1], origin[1] + sizes[1], length=nx+1)
    y = range(origin[2], origin[2] + sizes[2], length=ny+1)
    ϕ_grid = reshape(ϕ, (nx+1, ny+1))
    
    title = isempty(label) ? "t=$t" : "$label (t=$t)"
    
    p = Plt.contour(x, y, ϕ_grid', levels=[0.0], linewidth=2, color=:black, title=title)
    Plt.contourf!(p, x, y, ϕ_grid', levels=[-1e10, 0.0], color=[:lightblue], alpha=0.4)
    Plt.display(p)
end

end # module
