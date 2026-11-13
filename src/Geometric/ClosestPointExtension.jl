module ClosestPointExtension

using ..GridInterpolation: getInterpolator
using ..Geometric: _gridMetrics, _checklevel
using ....EvolvingDomains: Field, _levelgrid, _levelvalues

"""
    extendClosestPoint!(values, levelset::Field, level::Symbol) -> values

Extend nodal `values` from the bulk (φ ≤ 0) into the void (φ > 0) along the level set's own
normal, in place, on the grid of `level` of `levelset`:

    u_ext(x) = u(x - φ(x) ∇φ(x) / |∇φ(x)|)

Bulk nodes keep their value; only void nodes are overwritten. `∇φ` is a central difference
clamped into the box at the edge nodes. The interpolant is taken over a **copy** of
`values`, so reading ahead of the writes is not a hazard.

**This is a projection along the level-set normal, not a distance solve.** It is the
closest-point map only where `|∇φ| ≈ 1`; on a level set that is not a distance function the
projected point misses the interface by `φ (|∇φ| - 1)/|∇φ|` and the extension is wrong by
that much. Run `reinitialize!` first. It is also first-order: one projection, so a curved
front is offset by roughly the truncation error of `∇φ`, and around a thin feature the
projected point can land beyond the bulk on the far side, where the interpolant returns
that side's value instead. Both are the method's, not the implementation's.

`values` must be finite everywhere, void included, even where the value is meaningless
before extension: interpolation weights every corner of the target cell. Only 2D and 3D on a
uniformly spaced level grid are supported, matching `reinitialize!`.
"""
function extendClosestPoint!(values::AbstractArray, levelset::Field{T,N},
                             level::Symbol) where {T,N}
    grid = _levelgrid(levelset.discretization, level)
    phi = _levelvalues(levelset, level)
    info = _gridMetrics(grid)
    _checklevel(phi, nothing, info, "Level-set Field values")
    _checklevel(values, nothing, info, "extended values")
    itp = getInterpolator(vec(values), info)   # copies: reads never see the writes below
    dims, origin, spacing = info.dims, info.origin, info.spacing

    for index in CartesianIndices(phi)
        distance = phi[index]
        distance <= 0 && continue                    # bulk: keep the original value
        x = ntuple(d -> Float64(origin[d] + (index[d] - 1) * spacing[d]), N)
        gradient = ntuple(N) do d
            c = index[d]
            lo, hi = max(c - 1, 1), min(c + 1, dims[d])
            neighbour = v -> CartesianIndex(ntuple(k -> k == d ? v : index[k], N))
            (phi[neighbour(hi)] - phi[neighbour(lo)]) / ((hi - lo) * spacing[d])
        end
        gnorm = hypot(gradient...)
        # A vanishing gradient carries no direction: leave the node alone rather than
        # project with amplified noise.
        x_projected = gnorm > 1e-12 ?
            ntuple(d -> x[d] - distance * gradient[d] / gnorm, N) : x
        values[index] = itp(x_projected...)
    end
    return values
end

end # module
