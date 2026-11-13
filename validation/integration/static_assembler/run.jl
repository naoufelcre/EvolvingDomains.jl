using Test

# This optional check deliberately loads the redesigned source module rather than
# `using EvolvingDomains`, which may resolve an older globally installed package.
include(joinpath(@__DIR__, "..", "..", "..", "src", "EvolvingDomains.jl"))
const ED = Main.EvolvingDomains

const STATIC_ASSEMBLER_PATH = let
    isempty(ARGS) ? get(ENV, "STATICASSEMBLER_PATH", "") : first(ARGS)
end
isempty(STATIC_ASSEMBLER_PATH) && error(
    "Set STATICASSEMBLER_PATH to the downstream StaticAssembler directory, or pass it as the first argument")
include(joinpath(abspath(STATIC_ASSEMBLER_PATH), "StaticAssembler.jl"))
include(joinpath(@__DIR__, "descriptors.jl"))
const SA = Main.StaticAssembler

struct _MassOperator <: SA.AbstractDifferentialOperator end
SA._op_pair(::_MassOperator, ::_MassOperator, row::Int, col::Int,
            sh_t, gx_t, gy_t, sh_c, gx_c, gy_c) = sh_t[row] * sh_c[col]

function _downstream_records(mesh_descriptor, cut_descriptor)
    boundary_descriptor = cut_descriptor.boundary
    mesh = SA.Mesh{Float64}(mesh_descriptor.mesh_id,
        mesh_descriptor.nx, mesh_descriptor.ny,
        mesh_descriptor.origin_x, mesh_descriptor.origin_y,
        mesh_descriptor.hx, mesh_descriptor.hy)
    boundary = SA.BoundaryData{Float64}(
        boundary_descriptor.facet_to_bgcell, boundary_descriptor.point_ids,
        boundary_descriptor.point_pointers, boundary_descriptor.normals,
        boundary_descriptor.physical_points, boundary_descriptor.reference_points)
    cut = SA.CutData{Float64}(cut_descriptor.mesh_id, cut_descriptor.bg_states,
        cut_descriptor.subcell_states, cut_descriptor.cell_to_bgcell,
        cut_descriptor.point_ids, cut_descriptor.point_pointers,
        cut_descriptor.physical_points, cut_descriptor.reference_points, boundary)
    return mesh, cut
end

function _phase_areas(cut)
    areas = Dict(Int8(-1) => 0.0, Int8(1) => 0.0)
    for subcell in eachindex(cut.subcell_states)
        state = cut.subcell_states[subcell]
        state in keys(areas) || continue
        ids = Int(cut.point_pointers[subcell]):(Int(cut.point_pointers[subcell + 1]) - 1)
        p, q, r = cut.physical_points[Int.(cut.point_ids[ids])]
        areas[state] += abs((q[1] - p[1]) * (r[2] - p[2]) -
                            (q[2] - p[2]) * (r[1] - p[1])) / 2
    end
    return areas
end

@testset "redesigned Cut descriptors assemble in downstream StaticAssembler" begin
    scenarios = (
        (name="cut-only", xs=[0.0, 1.0], xcut=0.5,
         cut_in_area=0.5, cut_out_area=0.5, active_area=0.5, boundary_measure=1.0),
        (name="mixed full/cut/outside", xs=collect(0.0:1.0:4.0), xcut=1.5,
         cut_in_area=2.0, cut_out_area=2.0, active_area=6.0, boundary_measure=4.0),
        (name="CutCells aligned grid", xs=[0.0, 0.5, 1.0], xcut=0.5,
         cut_in_area=0.5, cut_out_area=0.0, active_area=0.5, boundary_measure=1.0),
    )

    mass = _MassOperator()
    slot = SA.BilinearSlot(mass, mass, SA.ConstantCoefficient(1.0))
    form = SA.ElementVariationalForm(
        SA.FiniteElement(SA.Q1ReferenceQuadrilateral(), 1), (slot,), (), 0;
        boundary_bilinear=(slot,))

    for scenario in scenarios
        xs = scenario.xs
        grid = ED.Grid((xs, copy(xs)))
        levelset = [x - scenario.xcut for x in xs, y in xs]
        field = ED.Field(ED.Discretization(Dict(:init => grid)), Dict(:init => levelset))
        snapshot = ED.cut(field, :init)
        @testset "$(scenario.name): CutCells" begin
                mesh_descriptor, cut_descriptor = to_descriptors(grid, snapshot)
                mesh, cut = _downstream_records(mesh_descriptor, cut_descriptor)
                @test mesh isa SA.Mesh{Float64}
                @test cut isa SA.CutData{Float64}
                @test cut.boundary isa SA.BoundaryData{Float64}
                areas = _phase_areas(cut)
                @test isapprox(areas[Int8(-1)], scenario.cut_in_area; atol=1e-12)
                @test isapprox(areas[Int8(1)], scenario.cut_out_area; atol=1e-12)
                @test !isempty(cut.boundary.facet_to_bgcell)
                @test all(bgcell -> cut.bg_states[Int(bgcell)] == 0,
                          cut.cell_to_bgcell)
                @test all(n -> isapprox(n[1], 1.0; atol=1e-12) &&
                               isapprox(n[2], 0.0; atol=1e-12), cut.boundary.normals)

                if scenario.name == "mixed full/cut/outside"
                    @test all(state -> state in cut.bg_states, (Int8(-1), Int8(0), Int8(1)))
                elseif scenario.name == "CutCells aligned grid"
                    @test count(==(Int8(0)), cut.bg_states) == 2
                    @test count(==(Int8(1)), cut.bg_states) == 2
                    @test !any(==(Int8(-1)), cut.bg_states)
                    @test all(==(Int8(-1)), cut.subcell_states)
                    @test length(cut.boundary.facet_to_bgcell) == 2
                else
                    @test all(==(Int8(0)), cut.bg_states)
                end

                roots = Int32.(1:(mesh.nx * mesh.ny))
                plan = SA.StaticSystem(mesh, cut, roots; form)
                matrix, rhs = SA.assemble!(plan, (); threaded=false)
                @test size(matrix, 1) == size(matrix, 2)
                @test length(rhs) == size(matrix, 1)
                @test all(iszero, rhs)
                # Q1 partition of unity integrates active volume plus interface measure.
                @test isapprox(sum(matrix),
                               scenario.active_area + scenario.boundary_measure; atol=1e-11)
                @test isapprox(sum(plan.refresh.quadrature.weight),
                               scenario.cut_in_area; atol=1e-11)
                boundary_quadrature = plan.refresh.boundary_quadrature
                @test isapprox(sum(boundary_quadrature.weight),
                               scenario.boundary_measure; atol=1e-11)
                @test all(isapprox(nx, 1.0; atol=1e-12) && isapprox(ny, 0.0; atol=1e-12)
                          for (nx, ny) in zip(boundary_quadrature.normal_x,
                                              boundary_quadrature.normal_y))
            end
        end
    end
