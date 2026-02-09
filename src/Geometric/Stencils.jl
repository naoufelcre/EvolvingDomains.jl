module Stencils

using ..CartesianField: CartesianMeshField, meshsize

export D⁺, D⁻, D⁰, D²

# =============================================================================
# Finite Difference Primitives
# =============================================================================

@inline function _offset(dim)
    CartesianIndex(ntuple(i -> i == dim ? 1 : 0, 2))
end

"""
    D⁺(ϕ, I, dim)

Forward finite difference: (ϕ[I+1] - ϕ[I]) / h
"""
@inline function D⁺(ϕ::CartesianMeshField, I, dim)
    h = meshsize(ϕ, dim)
    Ip = I + _offset(dim)
    return (ϕ[Ip] - ϕ[I]) / h
end

"""
    D⁻(ϕ, I, dim)

Backward finite difference: (ϕ[I] - ϕ[I-1]) / h
"""
@inline function D⁻(ϕ::CartesianMeshField, I, dim)
    h = meshsize(ϕ, dim)
    Im = I - _offset(dim)
    return (ϕ[I] - ϕ[Im]) / h
end

"""
    D⁰(ϕ, I, dim)

Centered finite difference: (ϕ[I+1] - ϕ[I-1]) / 2h
"""
@inline function D⁰(ϕ::CartesianMeshField, I, dim)
    h = meshsize(ϕ, dim)
    Ip = I + _offset(dim)
    Im = I - _offset(dim)
    return (ϕ[Ip] - ϕ[Im]) / (2*h)
end

"""
    D²(ϕ, I, dim)

Second order centered difference: (ϕ[I+1] - 2ϕ[I] + ϕ[I-1]) / h^2
"""
@inline function D²(ϕ::CartesianMeshField, I, dim)
    h = meshsize(ϕ, dim)
    Ip = I + _offset(dim)
    Im = I - _offset(dim)
    return (ϕ[Ip] - 2*ϕ[I] + ϕ[Im]) / h^2
end

# =============================================================================
# Interpolation Kernels (for Semi-Lagrangian)
# =============================================================================

"""
    quadratic_interpolation_weights(α, stencil::Symbol)

Compute quadratic Lagrange interpolation weights for a fraction α ∈ [0, 1).

# Stencils
- `:left`: nodes {i-1, i, i+1}
- `:right`: nodes {i, i+1, i+2} (relative to i)
"""
@inline function quadratic_interpolation_weights(α, stencil::Symbol)
    if stencil == :left
        # Weights for [i-1, i, i+1]
        w1 = 0.5 * α * (α - 1)
        w2 = 1 - α^2
        w3 = 0.5 * α * (α + 1)
        return (w1, w2, w3)
    elseif stencil == :right
        # Weights for [i, i+1, i+2]
        # Shifted: evaluate polynomial at (1+α) using nodes {-1, 0, 1}? 
        # No, standard definition relative to i:
        # P(x) passing through f(i), f(i+1), f(i+2) evaluated at x = i+α
        # Let's derive it:
        # L_0 = (x - x1)(x - x2) / (x0 - x1)(x0 - x2) -> (α-1)(α-2)/(-1)(-2) = 0.5(α-1)(α-2)
        # L_1 = (x - x0)(x - x2) / (x1 - x0)(x1 - x2) -> (α)(α-2)/(1)(-1)     = -α(α-2) = α(2-α)
        # L_2 = (x - x0)(x - x1) / (x2 - x0)(x2 - x1) -> (α)(α-1)/(2)(1)      = 0.5α(α-1)
        
        w1 = 0.5 * (α - 1) * (α - 2)
        w2 = α * (2 - α)
        w3 = 0.5 * α * (α - 1)
        return (w1, w2, w3)
    else
        throw(ArgumentError("Unknown stencil: $stencil"))
    end
end

export quadratic_interpolation_weights

end # module
