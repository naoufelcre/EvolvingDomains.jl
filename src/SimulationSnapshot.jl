# =============================================================================
# Simulation Snapshotting for Visualization
# =============================================================================

"""
    SimulationFrame

A snapshot of simulation state at a specific time.

# Fields
- `t::Float64`: Simulation time
- `ϕ::Matrix{Float64}`: Level set values reshaped to grid (nx × ny)

# Example
```julia
frame = snapshot(eg)
@show frame.t
heatmap(frame.ϕ)
```
"""
struct SimulationFrame
    t::Float64
    ϕ::Matrix{Float64}
end

"""
    SimulationResult{N}

Container for a sequence of simulation snapshots with grid metadata.

# Fields
- `grid_info::CartesianGridInfo{N}`: Grid metadata (origin, spacing, dims)
- `frames::Vector{SimulationFrame}`: Cached simulation frames

# Example
```julia
frames = [snapshot(eg) for _ in 1:10 if (advance!(eg, Δt); true)]
result = SimulationResult(grid_info(eg), frames)
viewer(result)  # Opens interactive viewer (requires GLMakie)
```
"""
struct SimulationResult{N}
    grid_info::CartesianGridInfo{N}
    frames::Vector{SimulationFrame}
end

"""
    snapshot(eg::EvolvingDiscreteGeometry) -> SimulationFrame

Capture the current simulation state as a `SimulationFrame`.

Returns a frame containing the current time and level set values 
reshaped to the grid dimensions.

# Example
```julia
frames = SimulationFrame[]

for step in 1:nsteps
    # Your physics here...
    advance!(eg, Δt)
    
    # Cache every 5 steps
    if step % 5 == 0
        push!(frames, snapshot(eg))
    end
end

result = SimulationResult(grid_info(eg), frames)
viewer(result)  # Requires GLMakie
```
"""
function snapshot(eg::EvolvingDiscreteGeometry)
    info = grid_info(eg)
    t = current_time(eg)
    ϕ = current_levelset(eg)
    ϕ_2d = reshape(copy(ϕ), info.dims)  # copy to avoid aliasing
    return SimulationFrame(t, ϕ_2d)
end

# Utility: number of frames
Base.length(result::SimulationResult) = length(result.frames)

# Utility: iteration
Base.iterate(result::SimulationResult, state=1) = 
    state > length(result.frames) ? nothing : (result.frames[state], state + 1)

# Utility: indexing
Base.getindex(result::SimulationResult, i::Int) = result.frames[i]
Base.firstindex(result::SimulationResult) = 1
Base.lastindex(result::SimulationResult) = length(result.frames)
