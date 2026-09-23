module CIP

using ...Geometric: CartesianGridInfo, CartesianMeshField
using Gridap.TensorValues: VectorValue

export CIPCache, cip_advect!

"""
    CIPCache(source::CartesianMeshField{Float64})

Optional workspace for intensive CIP transport. Besides reusable scratch arrays, the
cache carries the interpolation profile `(cx, cy, cxy)` between calls. A copy of the
last transported values is retained so that external changes to `source.data`
automatically rebuild the profile before the next step.

Using no cache is valid, but reconstructs the profile from nodal values every call.
"""
mutable struct CIPCache
    grid::CartesianGridInfo
    cx::Vector{Float64}
    cy::Vector{Float64}
    cxy::Vector{Float64}
    tmp::Vector{Float64}
    tmp_cx::Vector{Float64}
    tmp_cy::Vector{Float64}
    tmp_cxy::Vector{Float64}
    snapshot::Vector{Float64}
end

function CIPCache(source::CartesianMeshField{Float64})
    _check_field_size(source)
    _check_grid(source.grid)
    n = length(source.data)
    cache = CIPCache(
        source.grid,
        zeros(n), zeros(n), zeros(n),
        zeros(n), zeros(n), zeros(n), zeros(n),
        copy(source.data),
    )
    _reconstruct_profile!(cache, source.data)
    return cache
end

@inline _index(i, j, nx) = i + (j - 1) * nx

@inline function _same_grid(a::CartesianGridInfo, b::CartesianGridInfo)
    return a.origin == b.origin && a.spacing == b.spacing &&
           a.dims == b.dims && a.cells == b.cells
end

function _check_field_size(field::CartesianMeshField)
    expected = prod(field.grid.dims)
    length(field.data) == expected || throw(DimensionMismatch(
        "CartesianMeshField has $(length(field.data)) values; grid requires $expected."))
    return nothing
end

function _check_grid(grid::CartesianGridInfo)
    nx, ny = grid.dims
    nx >= 2 && ny >= 2 || throw(ArgumentError(
        "CIP requires at least two nodes in each grid direction."))
    all(h -> isfinite(h) && h > 0, grid.spacing) || throw(ArgumentError(
        "CIP requires finite positive grid spacing."))
    all(1:2) do d
        origin = grid.origin[d]
        spacing = grid.spacing[d]
        n = grid.dims[d]
        upper = origin + (n - 1) * spacing
        return isfinite(origin) && isfinite(upper) && origin + spacing > origin &&
               upper > origin + (n - 2) * spacing
    end || throw(ArgumentError(
        "CIP grid nodes must be finite and distinctly representable."))
    return nothing
end

@inline function _derivative_x(values, i, j, nx, dx)
    if i == 1
        nx > 2 && return (4 * (values[_index(2, j, nx)] - values[_index(1, j, nx)]) -
                              (values[_index(3, j, nx)] - values[_index(1, j, nx)])) / (2dx)
        return (values[_index(2, j, nx)] - values[_index(1, j, nx)]) / dx
    elseif i == nx
        nx > 2 && return (4 * (values[_index(nx, j, nx)] - values[_index(nx - 1, j, nx)]) -
                              (values[_index(nx, j, nx)] - values[_index(nx - 2, j, nx)])) / (2dx)
        return (values[_index(nx, j, nx)] - values[_index(nx - 1, j, nx)]) / dx
    end
    return (values[_index(i + 1, j, nx)] - values[_index(i - 1, j, nx)]) / (2dx)
end

@inline function _derivative_y(values, i, j, nx, ny, dy)
    if j == 1
        ny > 2 && return (4 * (values[_index(i, 2, nx)] - values[_index(i, 1, nx)]) -
                              (values[_index(i, 3, nx)] - values[_index(i, 1, nx)])) / (2dy)
        return (values[_index(i, 2, nx)] - values[_index(i, 1, nx)]) / dy
    elseif j == ny
        ny > 2 && return (4 * (values[_index(i, ny, nx)] - values[_index(i, ny - 1, nx)]) -
                              (values[_index(i, ny, nx)] - values[_index(i, ny - 2, nx)])) / (2dy)
        return (values[_index(i, ny, nx)] - values[_index(i, ny - 1, nx)]) / dy
    end
    return (values[_index(i, j + 1, nx)] - values[_index(i, j - 1, nx)]) / (2dy)
end

function _reconstruct_profile!(cache::CIPCache, values::Vector{Float64})
    nx, ny = cache.grid.dims
    dx, dy = cache.grid.spacing

    @inbounds for j in 1:ny, i in 1:nx
        k = _index(i, j, nx)
        cache.cx[k] = _derivative_x(values, i, j, nx, dx)
        cache.cy[k] = _derivative_y(values, i, j, nx, ny, dy)
    end

    # Average the two equivalent discrete definitions of the mixed derivative.
    @inbounds for j in 1:ny, i in 1:nx
        k = _index(i, j, nx)
        cache.cxy[k] = 0.5 * (
            _derivative_y(cache.cx, i, j, nx, ny, dy) +
            _derivative_x(cache.cy, i, j, nx, dx)
        )
    end
    copyto!(cache.snapshot, values)
    _profile_isfinite(cache) || throw(ArgumentError(
        "CIP could not reconstruct a finite interpolation profile from the source."))
    return cache
end

@inline _cache_arrays(cache::CIPCache) = (
    cache.cx, cache.cy, cache.cxy, cache.tmp, cache.tmp_cx, cache.tmp_cy,
    cache.tmp_cxy, cache.snapshot,
)

function _check_cache(cache::CIPCache, n)
    arrays = _cache_arrays(cache)
    all(values -> length(values) == n, arrays) || throw(DimensionMismatch(
        "CIP cache storage does not match the Cartesian grid."))
    for i in eachindex(arrays), j in (i + 1):length(arrays)
        Base.mightalias(arrays[i], arrays[j]) && throw(ArgumentError(
            "CIP cache arrays must not alias each other."))
    end
    return nothing
end

@inline function _profile_isfinite(cache::CIPCache)
    return all(isfinite, cache.cx) && all(isfinite, cache.cy) &&
           all(isfinite, cache.cxy)
end

@inline function _temporary_profile_isfinite(cache::CIPCache)
    return all(isfinite, cache.tmp) && all(isfinite, cache.tmp_cx) &&
           all(isfinite, cache.tmp_cy) && all(isfinite, cache.tmp_cxy)
end

@inline function _velocity_derivative_x(velocity, component, i, j, nx, dx)
    if i == 1
        nx > 2 && return (4 * (velocity[_index(2, j, nx)][component] -
                               velocity[_index(1, j, nx)][component]) -
                              (velocity[_index(3, j, nx)][component] -
                               velocity[_index(1, j, nx)][component])) / (2dx)
        return (velocity[_index(2, j, nx)][component] -
                velocity[_index(1, j, nx)][component]) / dx
    elseif i == nx
        nx > 2 && return (4 * (velocity[_index(nx, j, nx)][component] -
                               velocity[_index(nx - 1, j, nx)][component]) -
                              (velocity[_index(nx, j, nx)][component] -
                               velocity[_index(nx - 2, j, nx)][component])) / (2dx)
        return (velocity[_index(nx, j, nx)][component] -
                velocity[_index(nx - 1, j, nx)][component]) / dx
    end
    return (velocity[_index(i + 1, j, nx)][component] -
            velocity[_index(i - 1, j, nx)][component]) / (2dx)
end

@inline function _velocity_derivative_y(velocity, component, i, j, nx, ny, dy)
    if j == 1
        ny > 2 && return (4 * (velocity[_index(i, 2, nx)][component] -
                               velocity[_index(i, 1, nx)][component]) -
                              (velocity[_index(i, 3, nx)][component] -
                               velocity[_index(i, 1, nx)][component])) / (2dy)
        return (velocity[_index(i, 2, nx)][component] -
                velocity[_index(i, 1, nx)][component]) / dy
    elseif j == ny
        ny > 2 && return (4 * (velocity[_index(i, ny, nx)][component] -
                               velocity[_index(i, ny - 1, nx)][component]) -
                              (velocity[_index(i, ny, nx)][component] -
                               velocity[_index(i, ny - 2, nx)][component])) / (2dy)
        return (velocity[_index(i, ny, nx)][component] -
                velocity[_index(i, ny - 1, nx)][component]) / dy
    end
    return (velocity[_index(i, j + 1, nx)][component] -
            velocity[_index(i, j - 1, nx)][component]) / (2dy)
end

@inline function _velocity_mixed_derivative(velocity, component, i, j, nx, ny, dx, dy)
    if j == 1
        ny > 2 && return (
            4 * (_velocity_derivative_x(velocity, component, i, 2, nx, dx) -
                 _velocity_derivative_x(velocity, component, i, 1, nx, dx)) -
                (_velocity_derivative_x(velocity, component, i, 3, nx, dx) -
                 _velocity_derivative_x(velocity, component, i, 1, nx, dx))
        ) / (2dy)
        return (_velocity_derivative_x(velocity, component, i, 2, nx, dx) -
                _velocity_derivative_x(velocity, component, i, 1, nx, dx)) / dy
    elseif j == ny
        ny > 2 && return (
            4 * (_velocity_derivative_x(velocity, component, i, ny, nx, dx) -
                 _velocity_derivative_x(velocity, component, i, ny - 1, nx, dx)) -
                (_velocity_derivative_x(velocity, component, i, ny, nx, dx) -
                 _velocity_derivative_x(velocity, component, i, ny - 2, nx, dx))
        ) / (2dy)
        return (_velocity_derivative_x(velocity, component, i, ny, nx, dx) -
                _velocity_derivative_x(velocity, component, i, ny - 1, nx, dx)) / dy
    end
    return (_velocity_derivative_x(velocity, component, i, j + 1, nx, dx) -
            _velocity_derivative_x(velocity, component, i, j - 1, nx, dx)) / (2dy)
end

"""Hermite value and its first two derivatives inside one grid interval."""
@inline function _hermite(f0, g0, f1, g1, h, distance)
    t = distance / h
    t2 = t * t
    t3 = t2 * t
    delta = f1 - f0

    value = f0 + (-2t3 + 3t2) * delta + (t3 - 2t2 + t) * h * g0 +
            (t3 - t2) * h * g1
    first = ((-6t2 + 6t) / h) * delta + (3t2 - 4t + 1) * g0 +
            (3t2 - 2t) * g1
    second = ((-12t + 6) / h^2) * delta + ((6t - 4) / h) * g0 +
             ((6t - 2) / h) * g1
    return value, first, second
end

@inline function _departure_cell(x, origin, spacing, n)
    upper = origin + (n - 1) * spacing
    coordinate_ulp = max(eps(abs(origin)), eps(abs(upper)))
    tolerance = min(4coordinate_ulp, 1e-8 * spacing)
    if x < origin - tolerance || x > upper + tolerance
        throw(DomainError(x,
            "CIP departure point lies outside [$origin, $upper]; outer-boundary inflow is not supported."))
    end

    bounded = clamp(x, origin, upper)
    coordinate = (bounded - origin) / spacing
    left = min(floor(Int, coordinate) + 1, n - 1)
    distance = bounded - (origin + (left - 1) * spacing)
    return left, distance
end

function _check_departures(velocity, grid::CartesianGridInfo, dt)
    nx, ny = grid.dims
    x0, y0 = grid.origin
    dx, dy = grid.spacing
    @inbounds for j in 1:ny
        previous = -Inf
        for i in 1:nx
            k = _index(i, j, nx)
            x = x0 + (i - 1) * dx
            departure = x - dt * velocity[k][1]
            _departure_cell(departure, x0, dx, nx)
            departure > previous || throw(DomainError(dt,
                "CIP x-characteristics cross; reduce the time step."))
            1 - dt * _velocity_derivative_x(velocity, 1, i, j, nx, dx) > 0 ||
                throw(DomainError(dt,
                    "CIP x-characteristic Jacobian is non-positive; reduce the time step."))
            previous = departure
        end
    end
    @inbounds for i in 1:nx
        previous = -Inf
        for j in 1:ny
            k = _index(i, j, nx)
            y = y0 + (j - 1) * dy
            departure = y - dt * velocity[k][2]
            _departure_cell(departure, y0, dy, ny)
            departure > previous || throw(DomainError(dt,
                "CIP y-characteristics cross; reduce the time step."))
            1 - dt * _velocity_derivative_y(velocity, 2, i, j, nx, ny, dy) > 0 ||
                throw(DomainError(dt,
                    "CIP y-characteristic Jacobian is non-positive; reduce the time step."))
            previous = departure
        end
    end
    return nothing
end

function _x_sweep!(cache::CIPCache, source, velocity, dt)
    nx, ny = cache.grid.dims
    x0, _ = cache.grid.origin
    dx, dy = cache.grid.spacing

    @inbounds for j in 1:ny, i in 1:nx
        k = _index(i, j, nx)
        x = x0 + (i - 1) * dx
        departure = x - dt * velocity[k][1]
        left_i, distance = _departure_cell(departure, x0, dx, nx)
        left = _index(left_i, j, nx)
        right = left + 1

        value, cx_departure, cxx_departure = _hermite(
            source[left], cache.cx[left], source[right], cache.cx[right], dx, distance)
        cy_departure, cxy_departure, _ = _hermite(
            cache.cy[left], cache.cxy[left], cache.cy[right], cache.cxy[right], dx, distance)

        ux = _velocity_derivative_x(velocity, 1, i, j, nx, dx)
        uy = _velocity_derivative_y(velocity, 1, i, j, nx, ny, dy)
        uxy = _velocity_mixed_derivative(velocity, 1, i, j, nx, ny, dx, dy)
        departure_x = 1 - dt * ux
        departure_y = -dt * uy
        departure_xy = -dt * uxy

        cache.tmp[k] = value
        cache.tmp_cx[k] = cx_departure * departure_x
        cache.tmp_cy[k] = cy_departure + cx_departure * departure_y
        cache.tmp_cxy[k] =
            (cxx_departure * departure_y + cxy_departure) * departure_x +
            cx_departure * departure_xy
    end
    return nothing
end

function _y_sweep!(target, cache::CIPCache, velocity, dt)
    nx, ny = cache.grid.dims
    _, y0 = cache.grid.origin
    dx, dy = cache.grid.spacing

    @inbounds for j in 1:ny, i in 1:nx
        k = _index(i, j, nx)
        y = y0 + (j - 1) * dy
        departure = y - dt * velocity[k][2]
        lower_j, distance = _departure_cell(departure, y0, dy, ny)
        lower = _index(i, lower_j, nx)
        upper = lower + nx

        value, cy_departure, cyy_departure = _hermite(
            cache.tmp[lower], cache.tmp_cy[lower],
            cache.tmp[upper], cache.tmp_cy[upper], dy, distance)
        cx_departure, cxy_departure, _ = _hermite(
            cache.tmp_cx[lower], cache.tmp_cxy[lower],
            cache.tmp_cx[upper], cache.tmp_cxy[upper], dy, distance)

        vx = _velocity_derivative_x(velocity, 2, i, j, nx, dx)
        vy = _velocity_derivative_y(velocity, 2, i, j, nx, ny, dy)
        vxy = _velocity_mixed_derivative(velocity, 2, i, j, nx, ny, dx, dy)
        departure_x = -dt * vx
        departure_y = 1 - dt * vy
        departure_xy = -dt * vxy

        target[k] = value
        cache.cx[k] = cx_departure + cy_departure * departure_x
        cache.cy[k] = cy_departure * departure_y
        cache.cxy[k] =
            (cyy_departure * departure_x + cxy_departure) * departure_y +
            cy_departure * departure_xy
    end
    return nothing
end

"""
    cip_advect!(target, source, velocity, dt; cache=nothing)

Advect an intensive scalar satisfying `∂c/∂t + v⋅∇c = 0` with first-order,
x-then-y Lie-split CIP. Each directional characteristic uses the Euler departure
`x - dt*v(x)`. `velocity` is a frozen nodal vector field for this step. Transport
is performed on the complete Cartesian grid, independently of any active geometry mask.

When `cache` is omitted, the interpolation profile is reconstructed from `source`
on every call. A supplied `CIPCache` carries the profile between calls and avoids
scratch allocations. If `source` differs from the cache's last output, its profile
is reconstructed automatically.

Departure points must remain inside the Cartesian grid. Standard CIP is not
monotone and can overshoot near unresolved or discontinuous data. Steps whose
discrete directional characteristic maps fold are rejected and must be subdivided.
"""
function cip_advect!(
    target::CartesianMeshField{Float64},
    source::CartesianMeshField{Float64},
    velocity::AbstractVector{<:VectorValue{2}},
    dt::Real;
    cache::Union{Nothing,CIPCache}=nothing,
)
    _check_field_size(target)
    _check_field_size(source)
    _check_grid(source.grid)
    _same_grid(target.grid, source.grid) || throw(ArgumentError(
        "CIP source and target must use the same Cartesian grid."))
    Base.mightalias(target.data, source.data) && throw(ArgumentError(
        "CIP source and target must not alias."))
    (Base.mightalias(target.data, velocity) || Base.mightalias(source.data, velocity)) &&
        throw(ArgumentError("CIP scalar fields must not alias velocity storage."))

    n = length(source.data)
    length(velocity) == n || throw(DimensionMismatch(
        "CIP velocity has $(length(velocity)) values; grid requires $n."))
    Base.require_one_based_indexing(velocity)
    all(isfinite, source.data) || throw(ArgumentError("CIP source contains non-finite values."))
    all(v -> isfinite(v[1]) && isfinite(v[2]), velocity) ||
        throw(ArgumentError("CIP velocity contains non-finite values."))
    isfinite(dt) || throw(ArgumentError("CIP time step must be finite."))
    dt >= 0 || throw(ArgumentError("CIP time step must be non-negative."))

    if !isnothing(cache)
        _same_grid(cache.grid, source.grid) || throw(ArgumentError(
            "CIP cache belongs to a different Cartesian grid."))
        _check_cache(cache, n)
        any(values -> Base.mightalias(target.data, values) || Base.mightalias(source.data, values),
            _cache_arrays(cache)) && throw(ArgumentError(
            "CIP source and target must not alias cache storage."))
    end

    stationary = iszero(dt) || all(v -> iszero(v[1]) && iszero(v[2]), velocity)
    if stationary && isnothing(cache)
        copyto!(target.data, source.data)
        return target
    end

    stationary || _check_departures(velocity, source.grid, Float64(dt))
    work = isnothing(cache) ? CIPCache(source) : cache
    if work.snapshot != source.data || !_profile_isfinite(work)
        try
            _reconstruct_profile!(work, source.data)
        catch
            work.snapshot[1] = NaN
            rethrow()
        end
    end

    if stationary
        copyto!(target.data, source.data)
        copyto!(work.snapshot, source.data)
        return target
    end

    _x_sweep!(work, source.data, velocity, Float64(dt))
    _temporary_profile_isfinite(work) || throw(ErrorException(
        "CIP produced non-finite values or profile moments; reduce the time step."))
    try
        _y_sweep!(work.snapshot, work, velocity, Float64(dt))
    catch
        work.snapshot[1] = NaN
        rethrow()
    end
    if !all(isfinite, work.snapshot) || !_profile_isfinite(work)
        work.snapshot[1] = NaN
        throw(ErrorException(
            "CIP produced non-finite values or profile moments; reduce the time step."))
    end
    copyto!(target.data, work.snapshot)
    return target
end

function cip_advect!(
    target::CartesianMeshField{Float64},
    source::CartesianMeshField{Float64},
    velocity::CartesianMeshField{T},
    dt::Real;
    cache::Union{Nothing,CIPCache}=nothing,
) where {T<:VectorValue{2}}
    _check_field_size(velocity)
    _same_grid(source.grid, velocity.grid) || throw(ArgumentError(
        "CIP velocity and scalar fields must use the same Cartesian grid."))
    return cip_advect!(target, source, velocity.data, dt; cache=cache)
end

end # module
