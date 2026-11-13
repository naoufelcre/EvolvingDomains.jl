module WENO5

using ...Geometric: _GridMetrics
using StaticArrays: SVector

"""Reusable RK3 workspace for WENO5 evolution.

A cache is mutable scratch storage for one evolution sequence. Do not share the same
cache between concurrent calls; pass a separate cache per task or omit it.
"""
mutable struct WENO5Cache
    rhs::Union{Nothing,Vector{Float64}}
    stage::Union{Nothing,Vector{Float64}}
    phi0::Union{Nothing,Vector{Float64}}

    WENO5Cache() = new(nothing, nothing, nothing)
end

"""Jiang & Shu (1996) fifth-order reconstruction from five stencil values.

The five values are the quantities being reconstructed (here the interval averages of
`∂φ/∂x`) and the result is their value at the right face of the middle one.
"""
@inline function _w5(v1, v2, v3, v4, v5)
    ϵ = 1e-6
    S1 = 13/12*(v1 - 2v2 + v3)^2 + 1/4*(v1 - 4v2 + 3v3)^2
    S2 = 13/12*(v2 - 2v3 + v4)^2 + 1/4*(v2 - v4)^2
    S3 = 13/12*(v3 - 2v4 + v5)^2 + 1/4*(3v3 - 4v4 + v5)^2
    α1 = 0.1/(S1 + ϵ)^2; α2 = 0.6/(S2 + ϵ)^2; α3 = 0.3/(S3 + ϵ)^2
    s = α1 + α2 + α3
    q1 = v1/3 - 7v2/6 + 11v3/6
    q2 = -v2/6 + 5v3/6 + v4/3
    q3 = v3/3 + 5v4/6 - v5/6
    w1 = α1/s; w2 = α2/s; w3 = α3/s
    return w1*q1 + w2*q2 + w3*q3
end

# `φ` at `idx`, with every index clamped into the box: the ghost value outside the domain is
# the edge value, i.e. the discrete zero-normal-gradient condition on φ.
@inline function _at(φ, idx::NTuple{N,Int}, dims::NTuple{N,Int}) where {N}
    linear = 1
    stride = 1
    @inbounds for d in 1:N
        linear += (clamp(idx[d], 1, dims[d]) - 1) * stride
        stride *= dims[d]
    end
    return @inbounds φ[linear]
end

@inline _shift(idx::NTuple{N,Int}, dim::Int, m::Int) where {N} =
    ntuple(d -> d == dim ? idx[d] + m : idx[d], N)

# The exact interval average of ∂φ/∂x_dim over the cell whose right face is `idx + m e_dim`:
# (φ(idx+m e)-φ(idx+(m-1)e))/h. Five of these are the five values `_w5` reconstructs, whose
# right-face value is then ∂φ/∂x at the node `idx`.
@inline function _average(φ, idx, dim, dims, h, m)
    return (_at(φ, _shift(idx, dim, m), dims) - _at(φ, _shift(idx, dim, m - 1), dims)) / h
end

# Upwinded fifth-order ∂φ/∂x_dim: left-biased for positive velocity, mirrored for negative.
@inline function _wenoDerivative(φ, idx, dim, dims, h, velocity_positive::Bool)
    if velocity_positive
        return _w5(_average(φ, idx, dim, dims, h, -2), _average(φ, idx, dim, dims, h, -1),
                   _average(φ, idx, dim, dims, h, 0), _average(φ, idx, dim, dims, h, 1),
                   _average(φ, idx, dim, dims, h, 2))
    end
    return _w5(_average(φ, idx, dim, dims, h, 3), _average(φ, idx, dim, dims, h, 2),
               _average(φ, idx, dim, dims, h, 1), _average(φ, idx, dim, dims, h, 0),
               _average(φ, idx, dim, dims, h, -1))
end

"""Right-hand side `-v⋅∇φ` of the level-set advection, 2D or 3D."""
function weno5Rhs!(rhs::Vector{Float64}, phi::Vector{Float64},
                   info::_GridMetrics{N}, velocity::Vector{SVector{N,Float64}}) where {N}
    dims, spacing = info.dims, info.spacing
    Base.Threads.@threads for last in 1:dims[N]
        @inbounds for rest in CartesianIndices(dims[1:N-1])
            idx = (Tuple(rest)..., last)
            linear = 1
            stride = 1
            for d in 1:N
                linear += (idx[d] - 1) * stride
                stride *= dims[d]
            end
            v = velocity[linear]
            dphi = 0.0
            for d in 1:N
                dphi += v[d] * _wenoDerivative(phi, idx, d, dims, spacing[d], v[d] > 0)
            end
            rhs[linear] = -dphi
        end
    end
    return rhs
end

function _wenoStep!(phi::Vector{Float64}, grid::_GridMetrics,
                     velocity, dt::Float64, cache)
    n = length(phi)
    if isnothing(cache.rhs) || length(cache.rhs) != n
        cache.rhs = zeros(Float64, n)
        cache.stage = zeros(Float64, n)
        cache.phi0 = zeros(Float64, n)
    end
    rhs, stage, phi0 = cache.rhs, cache.stage, cache.phi0
    @. phi0 = phi
    weno5Rhs!(rhs, phi, grid, velocity)
    @. stage = phi + dt*rhs
    @. phi = stage
    weno5Rhs!(rhs, phi, grid, velocity)
    @. stage = 0.75phi0 + 0.25stage + 0.25dt*rhs
    @. phi = stage
    weno5Rhs!(rhs, phi, grid, velocity)
    @. phi = (1/3)phi0 + (2/3)stage + (2/3)dt*rhs
    return phi
end

function weno5Step!(phi::Vector{Float64}, grid::_GridMetrics{N},
                     velocity::Vector{SVector{N,Float64}}, dt::Float64,
                     cache::Union{WENO5Cache,Nothing}=nothing) where {N}
    work = isnothing(cache) ? WENO5Cache() : cache
    return _wenoStep!(phi, grid, velocity, dt, work)
end

end # module
