module CIP

using ...Geometric: _GridMetrics, _gridMetrics, _checkmetrics, samemetrics
using ...EvolvingDomains: Grid
using StaticArrays: SVector

"""
    CIPCache(source::Vector{Float64}, grid::Grid)

Optional workspace for intensive CIP transport. `grid` must be a uniform 2D or 3D
Cartesian Grid. Besides reusable scratch arrays, the
cache carries the interpolation profile between calls (`cx`, `cy`, `cxy` in 2D and
all first and mixed derivatives in 3D). A copy of the last transported values is
retained so external changes to `source` rebuild the profile.

Using no cache is valid, but reconstructs the profile from nodal values every call.
"""
mutable struct CIPCache{N}
    grid::_GridMetrics{N}
    cx::Vector{Float64}
    cy::Vector{Float64}
    cxy::Vector{Float64}
    tmp::Vector{Float64}
    tmp_cx::Vector{Float64}
    tmp_cy::Vector{Float64}
    tmp_cxy::Vector{Float64}
    snapshot::Vector{Float64}
    cz::Vector{Float64}
    cxz::Vector{Float64}
    cyz::Vector{Float64}
    cxyz::Vector{Float64}
    tmp_cz::Vector{Float64}
    tmp_cxz::Vector{Float64}
    tmp_cyz::Vector{Float64}
    tmp_cxyz::Vector{Float64}
end

CIPCache(source::Vector{Float64}, grid::Grid) = CIPCache(source, _gridMetrics(grid))

function CIPCache(source::Vector{Float64}, info::_GridMetrics{N}) where {N}
    N in (2, 3) || throw(ArgumentError("CIP supports only 2D or 3D grids."))
    length(source) == prod(info.dims) || throw(DimensionMismatch(
        "CIP source has $(length(source)) values; grid requires $(prod(info.dims))."))
    _checkmetrics(info)
    n = length(source)
    # The z-family exists only for N == 3. Allocating it empty in 2D halves a 2D cache
    # instead of reserving workspace no 2D code path can reach.
    z = N == 3 ? ntuple(_ -> zeros(n), 8) : ntuple(_ -> Float64[], 8)
    cache = CIPCache{N}(
        info,
        zeros(n), zeros(n), zeros(n),
        zeros(n), zeros(n), zeros(n), zeros(n),
        copy(source),
        z...,
    )
    _reconstructProfile!(cache, source, info)
    return cache
end

@inline _index(i, j, nx) = i + (j - 1) * nx
@inline _index3(i, j, k, nx, ny) = i + (j - 1) * nx + (k - 1) * nx * ny

# Nodal ∂/∂x_axis: central difference inside, the three-point one-sided (4Δ₁ - Δ₂)/(2h) at
# the two edge nodes, and the plain two-node difference when the axis has only two nodes.
# The two-point neighbours are read only where the three-point stencil exists.
@inline _stride(dims, axis) = axis == 1 ? 1 : prod(dims[1:axis - 1])

function _derivative(values, axis, grid::_GridMetrics{N}) where {N}
    dims = grid.dims
    n = dims[axis]
    h = grid.spacing[axis]
    stride = _stride(dims, axis)
    out = similar(values, Float64)
    @inbounds for I in CartesianIndices(dims)
        index = 1
        s = 1
        for d in 1:N
            index += (I[d] - 1) * s
            s *= dims[d]
        end
        a = values[index]
        if I[axis] == 1
            b = values[index + stride]
            out[index] = n > 2 ? (4 * (b - a) - (values[index + 2stride] - a)) / (2h) :
                (b - a) / h
        elseif I[axis] == n
            b = values[index - stride]
            out[index] = n > 2 ? (4 * (a - b) - (a - values[index - 2stride])) / (2h) :
                (a - b) / h
        else
            out[index] = (values[index + stride] - values[index - stride]) / (2h)
        end
    end
    return out
end

# The interpolation profile of the source: nodal values plus the derivatives the sweeps
# transport. Mixed derivatives differentiate twice; `D_a` and `D_b` are Kronecker products on
# disjoint axes, so the order of the two applications is the identity, not a choice.
function _reconstructProfile!(cache::CIPCache{N}, values, grid::_GridMetrics{N}) where {N}
    copyto!(cache.snapshot, values)
    cache.cx .= _derivative(values, 1, grid)
    cache.cy .= _derivative(values, 2, grid)
    cache.cxy .= _derivative(cache.cx, 2, grid)
    if N == 3
        cache.cz .= _derivative(values, 3, grid)
        cache.cxz .= _derivative(cache.cx, 3, grid)
        cache.cyz .= _derivative(cache.cy, 3, grid)
        cache.cxyz .= _derivative(cache.cxy, 3, grid)
    end
    _profileIsFinite(cache) || throw(ArgumentError(
        "CIP could not reconstruct a finite interpolation profile from the source."))
    return cache
end

# Only the arrays that exist for the cache's dimension. The z-family is allocated empty in 2D,
# so a single all-16 view would make every 2D cache look inconsistent.
@inline _cacheArrays(cache::CIPCache{2}) = (
    cache.cx, cache.cy, cache.cxy, cache.tmp, cache.tmp_cx, cache.tmp_cy, cache.tmp_cxy,
    cache.snapshot,
)
@inline _cacheArrays(cache::CIPCache{3}) = (
    cache.cx, cache.cy, cache.cxy, cache.cz, cache.cxz, cache.cyz, cache.cxyz,
    cache.tmp, cache.tmp_cx, cache.tmp_cy, cache.tmp_cxy,
    cache.tmp_cz, cache.tmp_cxz, cache.tmp_cyz, cache.tmp_cxyz, cache.snapshot,
)

function _checkCache(cache::CIPCache, n)
    arrays = _cacheArrays(cache)
    all(values -> length(values) == n, arrays) || throw(DimensionMismatch(
        "CIP cache storage does not match the Cartesian grid."))
    for i in eachindex(arrays), j in (i + 1):length(arrays)
        Base.mightalias(arrays[i], arrays[j]) && throw(ArgumentError(
            "CIP cache arrays must not alias each other."))
    end
    return nothing
end

@inline function _profileIsFinite(cache::CIPCache)
    two_d = all(isfinite, cache.cx) && all(isfinite, cache.cy) && all(isfinite, cache.cxy)
    return cache.grid isa _GridMetrics{2} ? two_d :
        two_d && all(isfinite, cache.cz) && all(isfinite, cache.cxz) &&
        all(isfinite, cache.cyz) && all(isfinite, cache.cxyz)
end

@inline function _temporaryProfileIsFinite(cache::CIPCache)
    two_d = all(isfinite, cache.tmp) && all(isfinite, cache.tmp_cx) &&
            all(isfinite, cache.tmp_cy) && all(isfinite, cache.tmp_cxy)
    return cache.grid isa _GridMetrics{2} ? two_d :
        two_d && all(isfinite, cache.tmp_cz) && all(isfinite, cache.tmp_cxz) &&
        all(isfinite, cache.tmp_cyz) && all(isfinite, cache.tmp_cxyz)
end

@inline function _velocityDerivativeX(velocity, component, i, j, nx, dx)
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

@inline function _velocityDerivativeY(velocity, component, i, j, nx, ny, dy)
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

@inline function _velocityMixedDerivative(velocity, component, i, j, nx, ny, dx, dy)
    if j == 1
        ny > 2 && return (
            4 * (_velocityDerivativeX(velocity, component, i, 2, nx, dx) -
                 _velocityDerivativeX(velocity, component, i, 1, nx, dx)) -
                (_velocityDerivativeX(velocity, component, i, 3, nx, dx) -
                 _velocityDerivativeX(velocity, component, i, 1, nx, dx))
        ) / (2dy)
        return (_velocityDerivativeX(velocity, component, i, 2, nx, dx) -
                _velocityDerivativeX(velocity, component, i, 1, nx, dx)) / dy
    elseif j == ny
        ny > 2 && return (
            4 * (_velocityDerivativeX(velocity, component, i, ny, nx, dx) -
                 _velocityDerivativeX(velocity, component, i, ny - 1, nx, dx)) -
                (_velocityDerivativeX(velocity, component, i, ny, nx, dx) -
                 _velocityDerivativeX(velocity, component, i, ny - 2, nx, dx))
        ) / (2dy)
        return (_velocityDerivativeX(velocity, component, i, ny, nx, dx) -
                _velocityDerivativeX(velocity, component, i, ny - 1, nx, dx)) / dy
    end
    return (_velocityDerivativeX(velocity, component, i, j + 1, nx, dx) -
            _velocityDerivativeX(velocity, component, i, j - 1, nx, dx)) / (2dy)
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

@inline function _hermite3(f0, g0, f1, g1, h, distance)
    value, first, second = _hermite(f0, g0, f1, g1, h, distance)
    third = (-12 / h^3) * (f1 - f0) + (6 / h^2) * (g0 + g1)
    return value, first, second, third
end

@inline function _departureCell(x, origin, spacing, n)
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

function _checkDepartures(velocity, grid::_GridMetrics, dt)
    nx, ny = grid.dims
    x0, y0 = grid.origin
    dx, dy = grid.spacing
    @inbounds for j in 1:ny
        previous = -Inf
        for i in 1:nx
            k = _index(i, j, nx)
            x = x0 + (i - 1) * dx
            departure = x - dt * velocity[k][1]
            _departureCell(departure, x0, dx, nx)
            departure > previous || throw(DomainError(dt,
                "CIP x-characteristics cross; reduce the time step."))
            1 - dt * _velocityDerivativeX(velocity, 1, i, j, nx, dx) > 0 ||
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
            _departureCell(departure, y0, dy, ny)
            departure > previous || throw(DomainError(dt,
                "CIP y-characteristics cross; reduce the time step."))
            1 - dt * _velocityDerivativeY(velocity, 2, i, j, nx, ny, dy) > 0 ||
                throw(DomainError(dt,
                    "CIP y-characteristic Jacobian is non-positive; reduce the time step."))
            previous = departure
        end
    end
    return nothing
end

function _xSweep!(cache::CIPCache, source, velocity, dt)
    nx, ny = cache.grid.dims
    x0, _ = cache.grid.origin
    dx, dy = cache.grid.spacing

    @inbounds for j in 1:ny, i in 1:nx
        k = _index(i, j, nx)
        x = x0 + (i - 1) * dx
        departure = x - dt * velocity[k][1]
        left_i, distance = _departureCell(departure, x0, dx, nx)
        left = _index(left_i, j, nx)
        right = left + 1

        value, cx_departure, cxx_departure = _hermite(
            source[left], cache.cx[left], source[right], cache.cx[right], dx, distance)
        cy_departure, cxy_departure, _ = _hermite(
            cache.cy[left], cache.cxy[left], cache.cy[right], cache.cxy[right], dx, distance)

        ux = _velocityDerivativeX(velocity, 1, i, j, nx, dx)
        uy = _velocityDerivativeY(velocity, 1, i, j, nx, ny, dy)
        uxy = _velocityMixedDerivative(velocity, 1, i, j, nx, ny, dx, dy)
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

function _ySweep!(target, cache::CIPCache, velocity, dt)
    nx, ny = cache.grid.dims
    _, y0 = cache.grid.origin
    dx, dy = cache.grid.spacing

    @inbounds for j in 1:ny, i in 1:nx
        k = _index(i, j, nx)
        y = y0 + (j - 1) * dy
        departure = y - dt * velocity[k][2]
        lower_j, distance = _departureCell(departure, y0, dy, ny)
        lower = _index(i, lower_j, nx)
        upper = lower + nx

        value, cy_departure, cyy_departure = _hermite(
            cache.tmp[lower], cache.tmp_cy[lower],
            cache.tmp[upper], cache.tmp_cy[upper], dy, distance)
        cx_departure, cxy_departure, _ = _hermite(
            cache.tmp_cx[lower], cache.tmp_cxy[lower],
            cache.tmp_cx[upper], cache.tmp_cxy[upper], dy, distance)

        vx = _velocityDerivativeX(velocity, 2, i, j, nx, dx)
        vy = _velocityDerivativeY(velocity, 2, i, j, nx, ny, dy)
        vxy = _velocityMixedDerivative(velocity, 2, i, j, nx, ny, dx, dy)
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

@inline _bit(axis) = 1 << (axis - 1)

function _velocityPartials3(velocity, grid::_GridMetrics{3}, component)
    values = [Float64(v[component]) for v in velocity]
    parts = Dict{Tuple,Vector{Float64}}()
    for axis in 1:3
        parts[(axis,)] = _derivative(values, axis, grid)
    end
    for a in 1:3, b in a:3
        parts[(a, b)] = _derivative(parts[(a,)], b, grid)
    end
    for a in 1:3, b in a:3, c in b:3
        if a == 1 && b == 2 && c == 3
            orders = ((1, 2, 3), (1, 3, 2), (2, 1, 3), (2, 3, 1), (3, 1, 2), (3, 2, 1))
            # Differentiate each second-order permutation once in the remaining axis.
            vals = [begin
                first = _derivative(values, o[1], grid)
                second = _derivative(first, o[2], grid)
                _derivative(second, o[3], grid)
            end for o in orders]
            parts[(1, 2, 3)] = reduce(+, vals) ./ length(vals)
        else
            parts[(a, b, c)] = _derivative(parts[(a, b)], c, grid)
        end
    end
    return parts
end

function _qPartials3(velocity, grid::_GridMetrics{3}, axis, dt)
    parts = _velocityPartials3(velocity, grid, axis)
    for values in values(parts)
        values .*= -dt
    end
    parts[(axis,)] .+= 1.0
    return parts
end

function _profileArrays3(cache, temporary)
    return temporary ?
        (cache.tmp, cache.tmp_cx, cache.tmp_cy, cache.tmp_cxy, cache.tmp_cz,
         cache.tmp_cxz, cache.tmp_cyz, cache.tmp_cxyz) :
        (cache.snapshot, cache.cx, cache.cy, cache.cxy, cache.cz,
         cache.cxz, cache.cyz, cache.cxyz)
end

function _sweep3!(cache, velocity, dt, axis, grid::_GridMetrics{3})
    nx, ny, nz = grid.dims
    dims = grid.dims
    strides = (1, nx, nx * ny)
    h = grid.spacing[axis]
    origin = grid.origin[axis]
    b, c = filter(!=(axis), 1:3)
    bit_a, bit_b, bit_c = _bit(axis), _bit(b), _bit(c)
    bit_ab, bit_ac, bit_bc = bit_a | bit_b, bit_a | bit_c, bit_b | bit_c
    bit_abc = 7
    input = _profileArrays3(cache, false)
    output = _profileArrays3(cache, true)
    qparts = _qPartials3(velocity, grid, axis, dt)
    qfirst = ntuple(d -> qparts[(d,)], 3)
    qthird = qparts[(1,2,3)]
    stride = strides[axis]

    # Check each directional characteristic before any cache or target write.
    @inbounds for k in 1:nz, j in 1:ny, i in 1:nx
        coord = (i, j, k)[axis]
        index = _index3(i, j, k, nx, ny)
        x = origin + (coord - 1) * h
        departure = x - dt * velocity[index][axis]
        _departureCell(departure, origin, h, dims[axis])
        qa = qfirst[axis][index]
        qa > 0 || throw(DomainError(dt,
            "CIP $axis-characteristics have non-positive Jacobian; reduce the time step."))
        if coord > 1
            prior = index - stride
            prior_coord = origin + (coord - 2) * h
            prior_departure = prior_coord - dt * velocity[prior][axis]
            departure > prior_departure || throw(DomainError(dt,
                "CIP $axis-characteristics cross; reduce the time step."))
        end
    end

    @inbounds for k in 1:nz, j in 1:ny, i in 1:nx
        index = _index3(i, j, k, nx, ny)
        coord = (i, j, k)[axis]
        x = origin + (coord - 1) * h
        departure = x - dt * velocity[index][axis]
        left_coord, distance = _departureCell(departure, origin, h, dims[axis])
        left_index = index + (left_coord - coord) * stride
        right_index = left_index + stride

        # Each transverse derivative is itself a Hermite profile along this axis.
        function evaluate(mask)
            base = input[mask + 1]
            slope = input[(mask | bit_a) + 1]
            return _hermite3(base[left_index], slope[left_index],
                             base[right_index], slope[right_index], h, distance)
        end
        F = evaluate(0)
        Fb_profile = evaluate(bit_b)
        Fc_profile = evaluate(bit_c)
        Fbc_profile = evaluate(bit_bc)
        Fa, Faa, Faaa = F[2], F[3], F[4]
        Fb, Fab, Faab = Fb_profile[1], Fb_profile[2], Fb_profile[3]
        Fc, Fac, Faac = Fc_profile[1], Fc_profile[2], Fc_profile[3]
        Fbc, Fabc = Fbc_profile[1], Fbc_profile[2]
        q_a = qfirst[axis][index]
        q_b = qfirst[b][index]
        q_c = qfirst[c][index]
        q_ab = qparts[minmax(axis, b)][index]
        q_ac = qparts[minmax(axis, c)][index]
        q_bc = qparts[minmax(b, c)][index]
        q_abc = qthird[index]
        value = F[1]
        output[1][index] = value
        output[bit_a + 1][index] = Fa * q_a
        output[bit_b + 1][index] = Fb + Fa * q_b
        output[bit_c + 1][index] = Fc + Fa * q_c
        output[bit_ab + 1][index] = q_a * (Faa * q_b + Fab) + Fa * q_ab
        output[bit_ac + 1][index] = q_a * (Faa * q_c + Fac) + Fa * q_ac
        output[bit_bc + 1][index] = Fbc + Fab * q_c + Fac * q_b +
            Faa * q_b * q_c + Fa * q_bc
        output[bit_abc + 1][index] = q_a * (Faaa * q_b * q_c + Faab * q_c +
            Faac * q_b + Fabc + Faa * q_bc) +
            q_ab * (Faa * q_c + Fac) + q_ac * (Faa * q_b + Fab) +
            Fa * q_abc
    end
    return nothing
end

function _lagrangianTransport3D!(target::Vector{Float64}, velocity, dt, cache)
    grid = cache.grid
    for axis in 1:3
        _sweep3!(cache, velocity, dt, axis, grid)
        _temporaryProfileIsFinite(cache) || throw(ErrorException(
            "CIP produced non-finite values or profile moments; reduce the time step."))
        copyto!(cache.snapshot, cache.tmp)
        copyto!(cache.cx, cache.tmp_cx)
        copyto!(cache.cy, cache.tmp_cy)
        copyto!(cache.cz, cache.tmp_cz)
        copyto!(cache.cxy, cache.tmp_cxy)
        copyto!(cache.cxz, cache.tmp_cxz)
        copyto!(cache.cyz, cache.tmp_cyz)
        copyto!(cache.cxyz, cache.tmp_cxyz)
    end
    if !all(isfinite, cache.snapshot) || !_profileIsFinite(cache)
        cache.snapshot[1] = NaN
        throw(ErrorException("CIP produced non-finite values or profile moments; reduce the time step."))
    end
    copyto!(target, cache.snapshot)
    return target
end

"""Advect an intensive scalar with directionally split CIP on flat grid data."""
function cipLagrangianTransport!(target::Vector{Float64}, source::Vector{Float64},
        velocity::AbstractVector{<:SVector{N}}, info::_GridMetrics{N},
        dt::Real; cache::Union{Nothing,CIPCache}=nothing) where {N}
    N in (2, 3) || throw(ArgumentError("CIP supports only 2D or 3D grids."))
    length(target) == prod(info.dims) || throw(DimensionMismatch(
        "CIP target has $(length(target)) values; grid requires $(prod(info.dims))."))
    length(source) == prod(info.dims) || throw(DimensionMismatch(
        "CIP source has $(length(source)) values; grid requires $(prod(info.dims))."))
    _checkmetrics(info)
    Base.mightalias(target, source) && throw(ArgumentError(
        "CIP source and target must not alias."))
    (Base.mightalias(target, velocity) || Base.mightalias(source, velocity)) &&
        throw(ArgumentError("CIP scalar fields must not alias velocity storage."))
    n = length(source)
    length(velocity) == n || throw(DimensionMismatch(
        "CIP velocity has $(length(velocity)) values; grid requires $n."))
    Base.require_one_based_indexing(velocity)
    all(isfinite, source) || throw(ArgumentError("CIP source contains non-finite values."))
    all(v -> all(isfinite, v), velocity) ||
        throw(ArgumentError("CIP velocity contains non-finite values."))
    isfinite(dt) || throw(ArgumentError("CIP time step must be finite."))
    dt >= 0 || throw(ArgumentError("CIP time step must be non-negative."))
    if !isnothing(cache)
        cache isa CIPCache{N} || throw(ArgumentError(
            "CIP cache is $(length(cache.grid.dims))D but the grid is $N-D."))
        samemetrics(cache.grid, info) || throw(ArgumentError(
            "CIP cache belongs to a different Cartesian grid."))
        _checkCache(cache, n)
        any(values -> Base.mightalias(target, values) || Base.mightalias(source, values),
            _cacheArrays(cache)) && throw(ArgumentError(
            "CIP source and target must not alias cache storage."))
    end
    step = Float64(dt)
    stationary = iszero(step) || all(iszero, velocity)
    if stationary && isnothing(cache)
        copyto!(target, source); return target
    end
    N == 2 && !stationary && _checkDepartures(velocity, info, step)
    work = isnothing(cache) ? CIPCache(source, info) : cache
    if work.snapshot != source || !_profileIsFinite(work)
        try _reconstructProfile!(work, source, info) catch
            work.snapshot[1] = NaN; rethrow()
        end
    end
    if stationary
        copyto!(target, source); copyto!(work.snapshot, source); return target
    end
    if N == 2
        _xSweep!(work, source, velocity, step)
        _temporaryProfileIsFinite(work) || throw(ErrorException(
            "CIP produced non-finite values or profile moments; reduce the time step."))
        try _ySweep!(target, work, velocity, step) catch
            work.snapshot[1] = NaN; rethrow()
        end
        # The sweep left the transported derivatives in `cx/cy/cxy`; the snapshot must describe
        # the same field, or the next call reuses a profile that belongs to another time level.
        copyto!(work.snapshot, target)
        if !all(isfinite, target) || !_profileIsFinite(work)
            work.snapshot[1] = NaN
            throw(ErrorException(
                "CIP produced non-finite values or profile moments; reduce the time step."))
        end
    else
        try _lagrangianTransport3D!(target, velocity, step, work) catch
            work.snapshot[1] = NaN; rethrow()
        end
    end
    return target
end

end # module
