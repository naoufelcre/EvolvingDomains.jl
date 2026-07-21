using Gridap: FESpace, FEFunction, CellField, interpolate, Point, get_triangulation, get_free_dof_values
using ..Geometric: CartesianGridInfo, CartesianMeshField, get_interpolator
import TransferOperator
import Gridap

@inline _tmag(v::Real) = abs(v)
@inline _tmag(v) = sqrt(sum(abs2, v.data))

"""
    GridMeshTransfer

Operator to transfer fields between a Cartesian background grid and a 
finite element space (typically on a CutFEM mesh).
"""
struct GridMeshTransfer
    grid_info::CartesianGridInfo
    target_space::FESpace
    active_indices::Vector{Int}
    levelset::Union{Nothing,Vector{Float64}}
end

"""
Two constructors, two different operations. Choose by asking whether the FE function
means anything outside the physical domain Ω.

  * **3-arg (no level set)** — `prolong` evaluates at every active node, including
    nodes with φ > 0. Correct when `u_mesh` is a global interpolant of a field defined
    everywhere: its values outside Ω are meaningful, and `prolong` inverts `restrict`.
    This is what the transfer round-trip tests rely on.

  * **4-arg (with level set)** — `prolong` evaluates nodes outside Ω at their closest
    point on Γ instead. Correct when `u_mesh` solves a PDE posed only on Ω, where its
    values outside are an uncontrolled polynomial extension. Note this deliberately
    breaks the round-trip identity outside Ω, because there is nothing to round-trip
    to.

Picking the wrong one is silent: the 3-arg form on a PDE solution returned values 33x
the physical scale for Q2 with no error anywhere downstream.
"""
GridMeshTransfer(info::CartesianGridInfo, space::FESpace, active::Vector{Int}) =
    GridMeshTransfer(info, space, active, nothing)

"""
    restrict(op::GridMeshTransfer, u_grid::CartesianMeshField)

Map data from Eulerian Grid to FE Mesh (FESpace).
Leverages the grid's interpolator to provide a continuous CellField for Gridap.
"""
function TransferOperator.restrict(op::GridMeshTransfer, u_grid::CartesianMeshField)
    itp = get_interpolator(u_grid)
    # Wrap the interpolator in a CellField that Gridap can sample at quadrature points
    trian = get_triangulation(op.target_space)
    cf = CellField(x -> itp(x[1], x[2]), trian)
    
    # Project/Interpolate the CellField into the FE space
    return interpolate(cf, op.target_space)
end

"""
    prolong(op::GridMeshTransfer, u_mesh::FEFunction)

Map data from FE Mesh (FEFunction) to Eulerian Grid.
Optimized to use direct index mapping (Glue) if available, falling back to 
batch evaluation otherwise.
"""
function TransferOperator.prolong(op::GridMeshTransfer, u_mesh::FEFunction)
    nx, ny = op.grid_info.dims

    # With a level set available we always take the batch path: it is the only one
    # that can choose *where* to evaluate, which is the whole point of the fix.
    isnothing(op.levelset) || return _prolong_batch(op, u_mesh)
    
    # 1. Attempt High-Efficiency Path (Direct Index Mapping)
    # This avoids all geometric searches and point evaluations.
    trian = get_triangulation(u_mesh)
    
    # BodyFittedTriangulation from a GridPortion contains node_to_parent_node
    # which maps: active_node_idx -> background_node_idx
    if hasproperty(trian, :grid) && hasproperty(trian.grid, :node_to_parent_node)
        mapping = trian.grid.node_to_parent_node
        
        # Apply constraints if in AgFEMSpace to get actual nodal values
        V = u_mesh.fe_space
        u_vals = if hasproperty(V, :space)
            Vstd = V.space
            u_std = interpolate(u_mesh, Vstd)
            get_free_dof_values(u_std)
        else
            get_free_dof_values(u_mesh)
        end
        
        # Determine if we need to reinterpret (for VectorValue fields)
        # We compare the number of DOFs to the number of nodes.
        n_nodes = Gridap.Geometry.num_nodes(trian)
        n_dofs = length(u_vals)
        
        T_dof = eltype(u_vals)
        
        if n_dofs == n_nodes
            # Scalar case
            new_data = zeros(T_dof, nx * ny)
            new_data[mapping] .= u_vals
            return CartesianMeshField(new_data, op.grid_info)
        elseif n_dofs > n_nodes && n_dofs % n_nodes == 0
            # Vector case: DOFs are interleaved [u1, v1, u2, v2...]
            # We determine the dimension from the ratio
            dim = n_dofs ÷ n_nodes
            
            # Since we only support 2D currently, we can assume VectorValue{dim, T_dof}
            # To be robust, we create the type dynamically.
            T_vec = Gridap.TensorValues.VectorValue{dim, T_dof}
            u_vals_view = reinterpret(T_vec, u_vals)
            
            new_data = zeros(T_vec, nx * ny)
            new_data[mapping] .= u_vals_view
            return CartesianMeshField(new_data, op.grid_info)
        else
            # Unexpected DOF layout (e.g. Higher order or complex constraints)
            # Fallback to Batch Evaluation
            return _prolong_batch(op, u_mesh)
        end
    end

    return _prolong_batch(op, u_mesh)
end

function _prolong_batch(op::GridMeshTransfer, u_mesh::FEFunction)
    nx, ny = op.grid_info.dims
    origin = op.grid_info.origin
    dx, dy = op.grid_info.spacing
    
    if isempty(op.active_indices)
        error("Cannot prolong field: active_indices is empty. " *
              "Ensure setup_transfer(geom, space) was called after the levelset was set, " *
              "and that the geometry intersects the background grid.")
    end

    # Where to evaluate the FE function.
    #
    # `active_indices` holds the nodes of IN *and CUT* cells, so a good share of them
    # lie at φ > 0 — outside the physical domain. An AgFEM function is only an
    # approximation on Ω; outside it is the root cell's polynomial extended, which the
    # error estimate says nothing about and which grows like the Chebyshev bound
    # |T_k(1 + 2d/h)| with the polynomial degree. Measured on a clean Hele-Shaw
    # geometry: Q1 overshoots the physical scale by 1x, Q2 by 33x.
    #
    # So for a node outside Ω we do not ask the FE function about that node. We ask it
    # about the closest point on Γ,
    #
    #     x_Γ = x - φ(x) ∇φ/|∇φ|,
    #
    # which lies on ∂Ω where the trace *is* controlled. This is also exactly the
    # constant-along-normals velocity extension the level set wants (the solution of
    # ∇u·∇φ = 0), so it is the right object rather than a workaround. The two branches
    # agree as φ → 0, so no discontinuity is introduced across Γ.
    φ = op.levelset
    grad = isnothing(φ) ? nothing : _compute_fd_gradient(op.grid_info, φ)

    points = Vector{Point{2,Float64}}(undef, length(op.active_indices))
    outside = falses(length(op.active_indices))
    for (k, idx) in enumerate(op.active_indices)
        i = (idx - 1) % nx + 1
        j = (idx - 1) ÷ nx + 1
        px = origin[1] + (i-1)*dx
        py = origin[2] + (j-1)*dy
        if !isnothing(φ) && φ[idx] > 0
            g = grad[idx]
            gn = hypot(g[1], g[2])
            if gn > 1e-8                      # |∇φ| ≈ 0: no usable normal, leave as is
                # Step no further than the band. If φ is not a true distance function
                # (no reinitialisation, or drift between calls) the projection can
                # overshoot past the active region, and the FE evaluation then throws
                # "Point ... is not inside any active cell". Capping the step keeps the
                # query inside the cells the FE function is defined on; the value is then
                # a near-Γ trace rather than an exact one, which is still far better than
                # reading the extrapolant at the node.
                step = min(φ[idx], 4 * max(dx, dy))
                px -= step * g[1] / gn
                py -= step * g[2] / gn
                outside[k] = true
            end
        end
        points[k] = Point(px, py)
    end

    # Batch evaluate on the FEFunction.
    #
    # The projection can still land outside the active region when φ has drifted from a
    # distance function, and Gridap then asserts "Point ... is not inside any active
    # cell". Capping the step is not a guarantee, so fall back to evaluating at the
    # nodes themselves: that reintroduces the extrapolation this projection exists to
    # avoid, but a degraded field beats an exception, and the warning says which.
    vals = try
        u_mesh(points)
    catch err
        isnothing(φ) && rethrow()
        @warn "prolong: closest-point projection left the active region; \
               falling back to evaluation at nodes for this call" maxlog = 3
        fill!(outside, false)
        for (k, idx) in enumerate(op.active_indices)
            i = (idx - 1) % nx + 1
            j = (idx - 1) ÷ nx + 1
            points[k] = Point(origin[1] + (i-1)*dx, origin[2] + (j-1)*dy)
        end
        u_mesh(points)
    end

    # Guard. After projection every value is a trace on Γ, so the void nodes cannot
    # exceed the interior ones by much. This exact failure stayed silent through a
    # passing test suite and a validated benchmark, because nothing downstream looks
    # at the magnitude — hence the check lives in the operator, not in a test.
    if !isnothing(φ) && any(outside) && !all(outside)
        m_out = maximum(_tmag, @view vals[outside])
        m_in = maximum(_tmag, @view vals[.!outside])
        if m_out > 2 * m_in + 1e-12
            @warn "prolong: values outside Ω exceed interior values" m_out m_in maxlog = 3
        end
    end

    T = eltype(vals)
    new_data = zeros(T, nx * ny)
    new_data[op.active_indices] .= vals

    return CartesianMeshField(new_data, op.grid_info)
end
