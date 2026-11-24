# EvolvingDomains

Julia library for composable and minimal orchestration of simulation components for writing evolving domain problems.

Its design aim for minimality and composability in the following sense:
- **Minimalism**: The package public API is centered around 3 structs that are basically just organized arrays. By doing so, we aim to have mathematical structures of package component clearly exposed.
- **Composable**: Because all that the package manipulate is bare-bones arrays, it's relatively easy to use the package as the orchestrator between different solvers. We provide a shared representation, but you can (and should) compose with methods of your choice.

## Overview

Let $d=2$ or $3$ and $\mathcal{O} \subset \mathbb{R}^d$. Our aim is to solve problems of the form

$
\begin{align}
	\frac{d}{d\,t}S(t) = \mathcal{G}(t,S(t)),
\end{align}
$

where for all $t$ in $[0,T]$: $S(t) = (\Omega(t),s(t))$ with $\Omega(t) \subset \mathcal{O}$ is the **evolving domain** and $s(t,\mathbf{x}) \in \mathcal{B}(\Omega(t))$ is a **state function**. We assume that $\Omega(t)$ can be described as the sub-zero level-set of a function $\phi:\mathcal{O}\to \mathbb{R}$ : $\Omega(t) = \{ \mathbf{x} \in \mathcal{O}, \phi(t,\mathbf{x}) < 0 \}$.

### API overview

The library is organized around the following three types:

```julia
Grid(axes::AbstractVector...)          # or Grid((v1, v2, ...)); Cartesian, strictly increasing axes
    .axes::NTuple{D,Vector{T}}   .mesh_id::UInt64    # fresh id per construction, survives copy

Discretization(axes::AbstractVector...)        # shorthand for Dict(:init => Grid(axes))
Discretization(levels::AbstractDict)           # :Symbol => Grid; all levels share one physical box
    .levels::Dict{Symbol,Grid{D,T}}

Field(d::Discretization, data::AbstractDict)   # :Symbol => dense Array; retained, not copied
    .discretization   .data::Dict{Symbol,Array{T,N}}
```
Schematically, those are nested structures.

### Typical structure of a program

A program written using `EvolvingDomains.jl` take the following form:

```julia
# Initialization
geometry = Field(grid, initial_levelset)
state    = Field(grid, initial_state)

config = (; geometry, state, parameters, time = 0.0)

function dynamics(config)
    # May involve external PDE solvers, geometric operators,
    # constitutive laws, constraints, or prescribed forcing.
    return rates
end

# Explicit simulation loop
while config.time < final_time

    rates = dynamics(config)
    velocity = velocityFromRates(rates, config)

    geometry = eulerianTransport(config.geometry, velocity, dt)
    
    state = lagrangianTransport(config.state, velocity, dt)  # optional
    state = evolve_state(state, geometry, rates, dt)        # optional
    
    config = update(config; geometry, state, time = config.time + dt)
end
```
In particular, we think verifiability of LLMs generated code is critical. To this end we provide an [Agent Skill](https://agentskills.io/home) in [`skills/evolvingdomains/SKILL.md`](skills/evolvingdomains/SKILL.md) for guiding LLMs toward this program structure and making their output checkable. The skill was tested using agents running in the [Pi coding agent](https://pi.dev) running [Deepseek V4.1-1017 flash](
https://doi.org/10.48550/arXiv.2609.19969). Broader coverage of harnesses-model would be needed and any feedback in that direction would be greatly appreciated !

The library exposes functions to handle some of those operations:

### Velocity from rates. 

### Eulerian transport.

### Lagrangian transport.

### Configuration update.

----------------------

## Composable 

Interfaces with the following libraries were tested: 

Trixi: for eulerianTransport and lagrangianTransport
Ferrite: for FEM solving
Gridap ecosystem: for FEM solving
FreeFem++: for FEM solving

## Examples


