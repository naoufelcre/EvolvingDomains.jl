module ClosestPointExtension

using ..GridInterpolation: _getInterpolator
using ..Geometric: _gridMetrics, _checklevel
using ....EvolvingDomains: Field, _levelgrid, _levelvalues

"""
    extendClosestPoint!(values::AbstractArray, levelset::Field{T,N}, level::Symbol) where {T,N}

Inputs:
- `values`: finite rank-`N` Float64 nodal samples `u` matching the selected level grid; mutated only in the positive phase `φ > 0`.
- `levelset`: finite rank-`N` Float64 nodal values of `φ` on a uniform 2D or 3D grid.
- `level`: attached level identifier shared by the values and level set.

Returns:
- The mutated nodal array, with `u_ext(x) = u(x − φ∇φ/|∇φ|)` projected from the interface into `φ > 0`; negative-phase values are unchanged.

Requires: `|∇φ| ≈ 1` where projection is used; the level and array shapes must match the uniform grid.
Use: Extend nodal data from the negative phase by one level-set-normal projection.
Algorithm: Central differences estimate `∇φ`; clamped multilinear interpolation samples a copy of `values`.
"""
function extendClosestPoint!(values::AbstractArray, levelset::Field{T,N},
                             level::Symbol) where {T,N}
    grid = _levelgrid(levelset.discretization, level)
    phi = _levelvalues(levelset, level)
    info = _gridMetrics(grid)
    _checklevel(phi, nothing, info, "Level-set Field values")
    _checklevel(values, nothing, info, "extended values")
    itp = _getInterpolator(vec(values), info)   # copies: reads never see the writes below
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
