
using Test
using EvolvingDomains
using Gridap
using GridapEmbedded

@testset "Mixed Grid Compatibility" begin
    
    # 1. Define two different models
    # "LSM Model": Standard resolution for level set evolution (Cartesian required)
    lsm_model = CartesianDiscreteModel((-2,2,-2,2), (50,50))
    
    # "BG Model": Coarser resolution simulation mesh (Generic DiscreteModel)
    # We use Cartesian here but treat it as generic DiscreteModel in signatures
    bg_model = CartesianDiscreteModel((-2,2,-2,2), (20,20))

    initial_sdf(x) = sqrt(x[1]^2 + x[2]^2) - 1.0

    # 2. Test Dual-Model Constructor
    geom = EvolvingDiscreteGeometry(bg_model, lsm_model, initial_sdf)
    
    @test geom.lsm_model === lsm_model
    @test geom.bg_model === bg_model
    @test geom.dirty == true

    # 3. Test Geometry Generation (Interpolation)
    geo = current_geometry(geom)
    @test isa(geo, GridapEmbedded.LevelSetCutters.DiscreteGeometry)
    
    # Verify the geometry is defined on the bg_model, not lsm_model
    # We can check the number of cells in the background model of the geometry
    # Note: Accessing model from DiscreteGeometry directly is not reliable. 
    # Validating via current_cut output instead.
    # @test num_cells(geo.model) == num_cells(bg_model)
    # @test num_cells(geo.model) != num_cells(lsm_model)

    # 4. Test Time Evolution
    # Add a velocity source
    vel = TimeDependentVelocity((x,t) -> (1.0, 0.0))
    set_velocity!(geom, vel)
    
    advance!(geom, 0.1)
    
    # Check that cut is updated and still on bg_model
    cut_mesh = current_cut(geom)
    @test cut_mesh !== nothing
    
    # Verify triangulation is associated with bg_model
    # We verify the active/cut triangulation's background model
    trian = Triangulation(cut_mesh)
    @test num_cells(get_background_model(trian)) == num_cells(bg_model)

    # 5. Test Legacy Constructor (Backward Compatibility)
    geom_legacy = EvolvingDiscreteGeometry(lsm_model, initial_sdf)
    @test geom_legacy.bg_model === lsm_model
    @test geom_legacy.lsm_model === lsm_model
    geo_legacy = current_geometry(geom_legacy)
    # @test num_cells(get_background_model(geo_legacy)) == num_cells(lsm_model)
    
    # Check cut generation for legacy case
    cut_legacy = current_cut(geom_legacy)
    trian_legacy = Triangulation(cut_legacy)
    @test num_cells(get_background_model(trian_legacy)) == num_cells(lsm_model)

end
