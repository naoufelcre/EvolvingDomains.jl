using Test
using EvolvingDomains
using Gridap, GridapEmbedded
using LevelSetMethods

@testset "Simulation Snapshot API" begin
    # Setup minimal geometry
    domain = (-1.0, 1.0, -1.0, 1.0)
    partition = (20, 20)
    model = CartesianDiscreteModel(domain, partition)
    
    R = 0.3
    ϕ₀(x) = R - sqrt(x[1]^2 + x[2]^2)
    u(x) = (0.1, 0.0)
    
    evolver = LevelSetMethodsEvolver(;
        bg_model = model,
        initial_ls = ϕ₀,
        velocity = u,
        spatial_scheme = :WENO5,
        bc = :Neumann
    )
    eg = EvolvingDiscreteGeometry(evolver, model)
    info = grid_info(eg)
    
    @testset "snapshot() returns SimulationFrame" begin
        frame = snapshot(eg)
        @test frame isa SimulationFrame
        @test frame.t == EvolvingDomains.current_time(eg)
        @test size(frame.ϕ) == info.dims
    end
    
    @testset "SimulationFrame fields" begin
        frame = snapshot(eg)
        @test hasfield(SimulationFrame, :t)
        @test hasfield(SimulationFrame, :ϕ)
        @test frame.t isa Float64
        @test frame.ϕ isa Matrix{Float64}
    end
    
    @testset "snapshot captures current state" begin
        frame1 = snapshot(eg)
        advance!(eg, 0.01)
        frame2 = snapshot(eg)
        
        @test frame2.t > frame1.t
        @test frame1.ϕ != frame2.ϕ  # Level set should have changed
    end
    
    @testset "snapshot copies data (no aliasing)" begin
        frame1 = snapshot(eg)
        ϕ_before = copy(frame1.ϕ)
        advance!(eg, 0.01)
        frame2 = snapshot(eg)
        
        # Frame1 should not have changed
        @test frame1.ϕ == ϕ_before
    end
    
    @testset "SimulationResult construction" begin
        frames = SimulationFrame[]
        for _ in 1:5
            push!(frames, snapshot(eg))
            advance!(eg, 0.01)
        end
        
        result = SimulationResult(info, frames)
        @test result isa SimulationResult
        @test length(result) == 5
        @test result.grid_info === info
    end
    
    @testset "SimulationResult iteration" begin
        frames = [snapshot(eg) for _ in 1:3 if (advance!(eg, 0.01); true)]
        result = SimulationResult(info, frames)
        
        # Test iteration
        count = 0
        for frame in result
            @test frame isa SimulationFrame
            count += 1
        end
        @test count == 3
        
        # Test indexing
        @test result[1] isa SimulationFrame
        @test result[end] === frames[end]
    end
end
