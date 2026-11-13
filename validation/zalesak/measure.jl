# Measuring a computed level set against the analytic one.
#
# Three plain measurements, no configuration objects. Every number a variant reports comes
# from here, so the sampling assumptions are stated once.

"""
    phase_measure(u, axes, exact, t; samples=4)

Negative-phase measure of `u` and its symmetric difference from the analytic level set, by
midpoint sampling of the multilinear cell reconstruction. `axes` is the tuple of
nodal axes, `exact` the analytic level set `(x..., t) -> value`.

Midpoint sampling is a diagnostic: it is not a conservative projection and not exact cut
geometry, so it is reported separately from native cut volumes.
"""
function phase_measure(u, axes, exact, t; samples=4)
    samples isa Integer && samples > 0 || throw(ArgumentError("samples must be positive"))
    D = length(axes)
    side = length(axes[1]) - 1
    all(axis -> length(axis) == side + 1, axes) ||
        throw(DimensionMismatch("expected a cubic nodal grid"))
    size(u) == ntuple(_ -> side + 1, D) || throw(DimensionMismatch("expected an x-first grid"))
    weight = prod(axis -> (axis[end] - axis[1])/side/samples, axes)
    offsets = ntuple(_ -> ((i - 0.5)/samples for i in 1:samples), D)
    inside = difference = 0
    for cell in Iterators.product(ntuple(_ -> 1:side, D)...)
        for w in Iterators.product(offsets...)
            phi = _interpolate(u, cell, w, Val(D))
            x = ntuple(d -> axes[d][cell[d]] +
                             w[d]*(axes[d][cell[d]+1] - axes[d][cell[d]]), D)
            truth = exact(x..., t) < 0
            inside += phi < 0
            difference += (phi < 0) != truth
        end
    end
    (measure=inside*weight, symmetric_difference=difference*weight)
end

function _interpolate(u, cell::NTuple{D,Int}, w::NTuple{D,Float64}, ::Val{D}) where D
    total = 0.0
    for corner in Iterators.product(ntuple(_ -> 0:1, D)...)
        weight = prod(d -> corner[d] == 0 ? 1 - w[d] : w[d], 1:D)
        total += weight * u[ntuple(d -> cell[d] + corner[d], D)...]
    end
    total
end

"""
    field_norms(u, exact, axes, center, radius)

Mean, RMS and maximum nodal error inside the ball of `radius` about `center`, where the
characteristics stay inside the box. Outside it, boundary treatment dominates.
"""
function field_norms(u, exact, axes, center, radius)
    errors = Float64[]
    for I in CartesianIndices(u)
        x = ntuple(d -> axes[d][I[d]], length(axes))
        hypot((x .- center)...) <= radius || continue
        push!(errors, abs(u[I] - exact[I]))
    end
    isempty(errors) && error("the norm region is empty")
    (l1=sum(errors)/length(errors), l2=sqrt(sum(abs2, errors)/length(errors)),
     linf=maximum(errors))
end

"""Edge zero crossings of a 2D level set: a sampled contour, not an exact Hausdorff bound."""
function interface_points(u, axis)
    points = Tuple{Float64,Float64}[]
    crossing(a, b, p, q) = if iszero(a)
        push!(points, p)
    elseif iszero(b)
        push!(points, q)
    elseif signbit(a) != signbit(b)
        r = a/(a - b)
        push!(points, (p[1] + r*(q[1]-p[1]), p[2] + r*(q[2]-p[2])))
    end
    for j in eachindex(axis), i in 1:length(axis)-1
        crossing(u[i,j], u[i+1,j], (axis[i], axis[j]), (axis[i+1], axis[j]))
    end
    for j in 1:length(axis)-1, i in eachindex(axis)
        crossing(u[i,j], u[i,j+1], (axis[i], axis[j]), (axis[i], axis[j+1]))
    end
    points
end

"Symmetric nearest-sample distance between two contours."
function interface_distance(a, b)
    (isempty(a) || isempty(b)) && return Inf
    directed(p, q) = maximum(x -> minimum(y -> hypot(x[1]-y[1], x[2]-y[2]), q), p)
    max(directed(a, b), directed(b, a))
end

"""One acceptance gate: print `name = value`, fail loudly if it leaves `[lower, upper]`."""
function within(value, name, lower, upper)
    ok = lower <= value <= upper
    println(rpad(string(name), 24), " = ", value, ok ? "" : "   OUT OF [$lower, $upper]")
    ok || error("gate failed: $name = $value ∉ [$lower, $upper]")
    value
end
