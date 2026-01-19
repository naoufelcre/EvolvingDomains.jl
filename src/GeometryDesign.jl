# =============================================================================
# GeometryDesign.jl — CSG Primitives for Level Set Initialization
# =============================================================================

using Gridap: VectorValue
using LinearAlgebra: norm

"""
    AbstractGeometry

Base type for CSG geometry primitives.
All geometries implement `signed_distance(geo, x) -> Float64`.
"""
abstract type AbstractGeometry end

"""
    signed_distance(geo::AbstractGeometry, x) -> Float64

Compute signed distance from point x to geometry surface.
Negative inside, positive outside.
"""
function signed_distance end

# Make geometries callable
(geo::AbstractGeometry)(x) = signed_distance(geo, x)

# =============================================================================
# Primitives
# =============================================================================

"""
    Circle(center::VectorValue, radius::Real)

2D circle primitive.
"""
struct Circle{T} <: AbstractGeometry
    center::VectorValue{2,T}
    radius::T
end

function signed_distance(c::Circle, x)
    return norm(x - c.center) - c.radius
end

"""
    Rectangle(pmin::VectorValue, pmax::VectorValue)

2D axis-aligned rectangle primitive.
"""
struct Rectangle{T} <: AbstractGeometry
    pmin::VectorValue{2,T}
    pmax::VectorValue{2,T}
end

function signed_distance(r::Rectangle, x)
    # Box SDF: max of distance to each face
    dx = max(r.pmin[1] - x[1], x[1] - r.pmax[1])
    dy = max(r.pmin[2] - x[2], x[2] - r.pmax[2])
    
    # Inside: max of negative distances
    # Outside: distance to corner or edge
    if dx <= 0 && dy <= 0
        return max(dx, dy)  # Inside: negative
    elseif dx <= 0
        return dy  # Above/below
    elseif dy <= 0
        return dx  # Left/right
    else
        return sqrt(dx^2 + dy^2)  # Corner
    end
end

# =============================================================================
# CSG Operations
# =============================================================================

struct GeometryUnion{A,B} <: AbstractGeometry
    a::A
    b::B
end

struct GeometryIntersection{A,B} <: AbstractGeometry
    a::A
    b::B
end

struct GeometryDifference{A,B} <: AbstractGeometry
    a::A
    b::B
end

signed_distance(g::GeometryUnion, x) = min(signed_distance(g.a, x), signed_distance(g.b, x))
signed_distance(g::GeometryIntersection, x) = max(signed_distance(g.a, x), signed_distance(g.b, x))
signed_distance(g::GeometryDifference, x) = max(signed_distance(g.a, x), -signed_distance(g.b, x))

# Standard Julia set operations
Base.union(a::AbstractGeometry, b::AbstractGeometry) = GeometryUnion(a, b)
Base.intersect(a::AbstractGeometry, b::AbstractGeometry) = GeometryIntersection(a, b)
Base.setdiff(a::AbstractGeometry, b::AbstractGeometry) = GeometryDifference(a, b)

# =============================================================================
# Transformations
# =============================================================================

"""
    Translate(geo::AbstractGeometry, offset::VectorValue)

Translate geometry by offset vector.
"""
struct Translate{G,T} <: AbstractGeometry
    geo::G
    offset::VectorValue{2,T}
end

function signed_distance(t::Translate, x)
    return signed_distance(t.geo, x - t.offset)
end
