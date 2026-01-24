
module IncrementalIntegration

using Gridap
using GridapIncremental
using GridapEmbedded
using ..EvolvingDomains: EvolvingDiscreteGeometry, current_cut

export IncrementalIntegrator
export update_integrator!
export measure_Ω, measure_Γ

"""
    mutable struct IncrementalIntegrator

Orchestrates the hybrid integration strategy.
Requires explicit caches for Scalar/Vector integrals (RHS) and Matrix integrals (LHS)
due to strict typing in CachedMeasure.
"""
mutable struct IncrementalIntegrator
    bg_measure_matrix::CachedMeasure # Matrix{Float64}
    bg_model::DiscreteModel
    prev_cell_status::Vector{Int8}
    degree::Int
    last_updated_step::Int
end

"""
    IncrementalIntegrator(model::DiscreteModel, degree::Int)

Constructs an integrator with a default background measure of given degree.
Optimized for Matrix Assembly (Left-Hand Side).
"""
function IncrementalIntegrator(model::DiscreteModel, degree::Int)
    trian = Triangulation(model)
    dΩ = Measure(trian, degree)

    cm_matrix = CachedMeasure(dΩ, Matrix{Float64})

    n_cells = num_cells(model)
    status = fill(Int8(-1), n_cells)
    return IncrementalIntegrator(cm_matrix, model, status, degree, -1)
end

"""
    update_integrator!(integrator::IncrementalIntegrator, geom::EvolvingDiscreteGeometry)

Updates the integrator state based on geometry classification on the grid changes.
We update the status history and let the `CachedMeasure` lazily fill any missing cells.
No invalidation is required for `IN` cells.
"""
function update_integrator!(integrator::IncrementalIntegrator, geom::EvolvingDiscreteGeometry)
    cut = current_cut(geom)

    # Access status safely
    if hasfield(typeof(cut), :ls_to_bgcell_to_inoutcut)
        # Assuming single level set for now
        current_status = cut.ls_to_bgcell_to_inoutcut[1]
    else
        # Fallback or error
        # Try cut.geo if it exists
         if hasfield(typeof(cut), :geo) && hasfield(typeof(cut.geo), :bgcell_to_inoutcut)
             current_status = cut.geo.bgcell_to_inoutcut
         else
             error("IncrementalIntegrator requires :ls_to_bgcell_to_inoutcut or :bgcell_to_inoutcut.")
         end
    end

    update_measure!(integrator.bg_measure_matrix, Int[])

    # Update History
    integrator.prev_cell_status .= current_status
    integrator.last_updated_step = geom.step

    return integrator
end

# ===================================================================

import Gridap.CellData: integrate

"""
    measure_Ω(integrator::IncrementalIntegrator, geom::EvolvingDiscreteGeometry)

Compute Volume Measure using Hybrid Strategy:
1. Fresh Cut Contribution (`dΩ_cut` on CUT cells).
2. Cached Bulk Contribution (`dΩ_bg` on IN cells, filtered).
"""
function measure_Ω(integrator::IncrementalIntegrator, geom::EvolvingDiscreteGeometry)
    # SAFETY CHECK
    if integrator.last_updated_step != geom.step
        error("IncrementalIntegrator Out-of-Sync! Integrator step: $(integrator.last_updated_step), Geometry step: $(geom.step). Call update_integrator!(integrator, geom) before integrating.")
    end
    return IncrementalVolumeMeasure(integrator, geom)
end

struct IncrementalVolumeMeasure <: Gridap.CellData.Measure
    integrator::IncrementalIntegrator
    geom::EvolvingDiscreteGeometry
end

# Overload integrate for this custom measure
function Gridap.CellData.integrate(f, m::IncrementalVolumeMeasure)
    integrator = m.integrator
    geom = m.geom

    cut = current_cut(geom)
    if hasfield(typeof(cut), :ls_to_bgcell_to_inoutcut)
        status = cut.ls_to_bgcell_to_inoutcut[1]
    else
        status = cut.geo.bgcell_to_inoutcut
    end

    # SAFETY CHECK
    if integrator.last_updated_step != geom.step
        error("IncrementalIntegrator Out-of-Sync! Integrator step: $(integrator.last_updated_step), Geometry step: $(geom.step). Call update_integrator!(integrator, geom) before integrating.")
    end

    IN = GridapEmbedded.Interfaces.IN
    CUT = GridapEmbedded.Interfaces.CUT
    PHYSICAL = GridapEmbedded.Interfaces.PHYSICAL

    degree = integrator.degree

    # 1. CUT contribution (Fresh, Lazy)
    # Optimization: Construct the view FIRST, then integrate.
    # This avoids generic integration over the full physical domain.
    trian_phys = Triangulation(cut, PHYSICAL)
    
    # Identify indices of CUT cells within the physical triangulation
    val_glue = Gridap.Geometry.get_glue(trian_phys, Val(2))
    active_to_bg = val_glue.tface_to_mface
    cut_indices = findall(bg_id -> status[bg_id] == CUT, active_to_bg)

    # Wrap in LazyTriangulation to avoid expensive connectivity maps
    # This view only contains CUT cells.
    lazy_trian_cut = LazyTriangulation(Triangulation(trian_phys, cut_indices))
    
    # Create Measure on the restricted view
    dΩ_cut = Measure(lazy_trian_cut, degree)

    # Integrate: Computes quadratures and evaluates f ONLY on CUT cells
    contrib_final = integrate(f, dΩ_cut)

    # 2. Bulk contribution (Cached)
    # Use the Matrix Cache
    contrib_bg_all = integrate(f, integrator.bg_measure_matrix)

    # Restrict to IN
    in_mask = findall(status .== IN)
    bg_trian = integrator.bg_measure_matrix.measure.quad.trian
    vals_bg = Gridap.CellData.get_contribution(contrib_bg_all, bg_trian)

    # Optimization: Only select needed values by lazy reindex
    pruned_vals_bg = Gridap.Arrays.lazy_map(Gridap.Arrays.Reindex(vals_bg), in_mask)
    
    # Use LazyTriangulation for IN contribution key as well
    trian_in = LazyTriangulation(Triangulation(bg_trian, in_mask))

    Gridap.CellData.add_contribution!(contrib_final, trian_in, pruned_vals_bg)

    return contrib_final
end

"""
    measure_Γ(integrator::IncrementalIntegrator, geom::EvolvingDiscreteGeometry)

Compute Boundary Measure (Fresh).
"""
function measure_Γ(integrator::IncrementalIntegrator, geom::EvolvingDiscreteGeometry)
    cut = current_cut(geom)
    # Boundary is the Interface
    Γ = EmbeddedBoundary(cut)

    # SAFETY CHECK
    if integrator.last_updated_step != geom.step
         error("IncrementalIntegrator Out-of-Sync! Integrator step: $(integrator.last_updated_step), Geometry step: $(geom.step). Call update_integrator!(integrator, geom) before integrating.")
    end

    degree = integrator.degree
    return Measure(Γ, degree)
end

end
