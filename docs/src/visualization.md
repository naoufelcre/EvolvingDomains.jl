# Visualization

visualization extension that activates automatically when [GLMakie](https://docs.makie.org/stable/) is loaded

## Setup

To use visualization features, simply load GLMakie alongside EvolvingDomains:

```julia
using EvolvingDomains
using GLMakie

# Create your evolving geometry...
evolver = LevelSetMethodsEvolver(; bg_model=model, initial_ls=ϕ₀, velocity=u)
eg = EvolvingDiscreteGeometry(evolver, model)

# Visualization functions are now available
fig = plot_levelset(eg)
```

When the extension loads, you'll see:
```
[ Info: EvolvingDomainsMakieExt loaded: plot_levelset(), viewer(), view_live!() available
```

## Static Plotting

### `plot_levelset`

Create a static plot of the current level set state.

```julia
plot_levelset(eg::EvolvingDiscreteGeometry; kwargs...) -> Figure
```

**Keyword Arguments:**
- `colormap::Symbol = :RdBu` - Colormap for level set values
- `show_zero::Bool = true` - Highlight the zero contour (domain boundary)
- `title::String = ""` - Plot title (auto-generated if empty)
- `filled::Bool = true` - Use filled heatmap for the level set

**Example:**

```julia
# Basic usage
fig = plot_levelset(eg)

# Custom styling
fig = plot_levelset(eg;
    colormap = :viridis,
    title = "My Domain",
    show_zero = true
)
display(fig)
```

### `plot_levelset!`

Add level set visualization to an existing Makie axis.

```julia
plot_levelset!(ax, eg::EvolvingDiscreteGeometry; kwargs...) -> Axis
```

**Example:**

```julia
fig = Figure(size = (1200, 400))
ax1 = Axis(fig[1, 1], title = "Before")
ax2 = Axis(fig[1, 2], title = "After")

# Initial state
plot_levelset!(ax1, eg)

# ... advance simulation ...
EvolvingDomains.advance!(eg, 0.5)

# Final state
plot_levelset!(ax2, eg)
display(fig)
```

## Simulation Snapshotting

To visualize time-dependent simulations, you first need to capture snapshots during the simulation:

```julia
using EvolvingDomains

# Run simulation and capture frames
frames = SimulationFrame[]
for step in 1:100
    EvolvingDomains.advance!(eg, Δt)
    push!(frames, snapshot(eg))
end

# Create result container
result = SimulationResult(grid_info(eg), frames)
```

The `snapshot()` function captures:
- Current time
- Level set values (copied)

## Interactive Viewer

### `viewer`

Open an interactive viewer with a time slider to browse simulation frames.

```julia
viewer(result::SimulationResult; kwargs...) -> Figure
```

**Keyword Arguments:**
- `colormap::Symbol = :RdBu` - Colormap for visualization
- `wait::Bool = false` - If true, block until Enter is pressed

**Example:**

```julia
# Run simulation with snapshots
frames = SimulationFrame[]
for step in 1:100
    EvolvingDomains.advance!(eg, Δt)
    if step % 5 == 0
        EvolvingDomains.reinitialize!(eg)
    end
    push!(frames, snapshot(eg))
end

result = SimulationResult(grid_info(eg), frames)

# Open interactive viewer
viewer(result)  # Drag the slider to navigate through time
```

The viewer displays:
- Heatmap of level set values
- Zero contour (domain boundary) in black
- Current time and frame number
- Interactive slider for frame selection

## Animated Playback

### `view_live!`

Play back cached simulation frames as an animation.

```julia
view_live!(result::SimulationResult; kwargs...) -> Nothing
```

**Keyword Arguments:**
- `fps::Real = 30` - Frames per second for playback
- `loop::Bool = false` - If true, loop animation continuously
- `colormap::Symbol = :RdBu` - Colormap for visualization

**Example:**

```julia
# Play animation at 15 fps
view_live!(result; fps = 15)

# Loop continuously (Ctrl+C to stop)
view_live!(result; fps = 20, loop = true)
```

!!! note
    `view_live!` blocks until playback completes. Use `loop = false` (default) for single playback, or press Ctrl+C to stop a looping animation.

## Complete Example

```julia
using EvolvingDomains
using Gridap, GridapEmbedded
using LevelSetMethods
using GLMakie

# Setup
domain = (-1.0, 1.0, -1.0, 1.0)
partition = (100, 100)
model = CartesianDiscreteModel(domain, partition)

ϕ₀(x) = 0.2 - sqrt(x[1]^2 + x[2]^2)
u(x) = (0.5, 0.3)

evolver = LevelSetMethodsEvolver(;
    bg_model = model,
    initial_ls = ϕ₀,
    velocity = u,
    spatial_scheme = :WENO5,
    bc = :Neumann
)
eg = EvolvingDiscreteGeometry(evolver, model)

# Show initial state
fig = plot_levelset(eg; title = "Initial Configuration")
display(fig)

# Run simulation with snapshots
frames = SimulationFrame[]
Δt = 0.01
for step in 1:100
    EvolvingDomains.advance!(eg, Δt)
    if step % 5 == 0
        EvolvingDomains.reinitialize!(eg)
    end
    push!(frames, snapshot(eg))
end

result = SimulationResult(grid_info(eg), frames)

# Interactive exploration
viewer(result)

# Or watch as animation
view_live!(result; fps = 20)
```
