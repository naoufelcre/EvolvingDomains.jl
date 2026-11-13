using Test
using EvolvingDomains
using EvolvingDomains.Geometric
using Gridap
using Gridap.TensorValues
using LinearAlgebra
using CairoMakie

function compute_grad_mag_field(phi, info)
    nx, ny = info.dims
    dx, dy = info.spacing
    mag = zeros(nx, ny)
    for j in 1:ny, i in 1:nx
        im, ip = max(1, i-1), min(nx, i+1)
        jm, jp = max(1, j-1), min(ny, j+1)
        gx = (phi[ip + (j-1)*nx] - phi[im + (j-1)*nx]) / ((ip-im)*dx)
        gy = (phi[i + (jp-1)*nx] - phi[i + (jm-1)*nx]) / ((jp-jm)*dy)
        mag[i, j] = sqrt(gx^2 + gy^2)
    end
    return mag
end

function save_reinit_plots(geom, name, output_dir; phi_ref=nothing)
    info = grid_info(geom.grid)
    nx, ny = info.dims
    ox, oy = info.origin
    dx, dy = info.spacing
    
    xs = range(ox, step=dx, length=nx)
    ys = range(oy, step=dy, length=ny)
    
    phi = current_levelset(geom)
    mag = compute_grad_mag_field(phi, info)
    
    num_cols = isnothing(phi_ref) ? 2 : 3
    fig = Figure(size=(500 * num_cols, 450))
    
    ax1 = Axis(fig[1, 1], aspect=1, title="$name: Level Set")
    contourf!(ax1, xs, ys, reshape(phi, nx, ny), levels=20, colormap=:RdBu)
    contour!(ax1, xs, ys, reshape(phi, nx, ny), levels=[0.0], color=:black, linewidth=2)
    
    ax2 = Axis(fig[1, 2], aspect=1, title="$name: |grad phi|")
    hm = heatmap!(ax2, xs, ys, mag, colormap=:viridis, colorrange=(0.7, 1.3))
    Colorbar(fig[2, 2], hm, vertical=false, label="Gradient Magnitude")
    
    if !isnothing(phi_ref)
        ax3 = Axis(fig[1, 3], aspect=1, title="$name: Error (vs GT)")
        err = reshape(phi .- vec(collect(phi_ref)), nx, ny)
        hm_err = heatmap!(ax3, xs, ys, err, colormap=:RdBu, colorrange=(-0.05, 0.05))
        Colorbar(fig[2, 3], hm_err, vertical=false, label="Signed Error")
    end
    
    save(joinpath(output_dir, "$name.png"), fig)
end

@testset "Reinitialization Test - Distorted Circle" begin
    output_dir = joinpath(@__DIR__, "output_reinit")
    mkpath(output_dir)

    n = 100
    grid = CartesianDiscreteModel((-1, 1, -1, 1), (n, n))
    geom = EvolvingDiscreteGeometry(grid)
    
    # Ground Truth: Unit circle SDF
    pts = Gridap.Geometry.get_node_coordinates(grid)
    phi_gt = vec(collect(map(x -> norm(x) - 0.5, pts)))
    
    # "Dirty" Initialization: Massive deformation of the distance field
    # We keep the zero level set at r=0.5, but distort the gradient elsewhere
    phi_init = vec(collect(map(pts) do x
        d = norm(x) - 0.5
        distortion = 2.0 + 0.5 * sin(3*x[1]) * cos(3*x[2])
        return d * distortion
    end))
    set_levelset!(geom, phi_init)
    
    save_reinit_plots(geom, "distorted_before", output_dir, phi_ref=phi_gt)
    
    # Reinitialize
    reinitialize!(geom)
    
    save_reinit_plots(geom, "distorted_after", output_dir, phi_ref=phi_gt)
    
    phi_new = current_levelset(geom)
    err = phi_new .- phi_gt
    l2_err = sqrt(sum(err.^2) / length(err))
    linf_err = maximum(abs.(err))
    
    println("Distorted Circle - L2 Error: $l2_err, Linf Error: $linf_err")
    
    # Assertions: For a 100x100 grid (dx=0.02), we expect error ~ O(dx)
    dx = 2.0 / n
    @test l2_err < 1.0 * dx
    @test linf_err < 2.5 * dx
    
    # Check gradient magnitude in the interior
    info = grid_info(grid)
    mag = compute_grad_mag_field(phi_new, info)
    avg_mag = sum(mag[5:end-5, 5:end-5]) / length(mag[5:end-5, 5:end-5])
    @test abs(avg_mag - 1.0) < 0.05
end

@testset "Merging Circles - Distorted" begin
    output_dir = joinpath(@__DIR__, "output_reinit")
    
    n = 100
    grid = CartesianDiscreteModel((-1, 1, -1, 1), (n, n))
    geom = EvolvingDiscreteGeometry(grid)
    
    pts = Gridap.Geometry.get_node_coordinates(grid)
    phi_gt = vec(collect(map(pts) do x
        d1 = norm(x - VectorValue(-0.3, 0.0)) - 0.25
        d2 = norm(x - VectorValue(0.3, 0.0)) - 0.25
        return min(d1, d2)
    end))
    
    # Distort the merging circles field
    phi_init = vec(collect(map(pts) do x
        d1 = norm(x - VectorValue(-0.3, 0.0)) - 0.25
        d2 = norm(x - VectorValue(0.3, 0.0)) - 0.25
        d = min(d1, d2)
        return d * (1.5 + 0.5 * sin(pi*x[1]))
    end))
    set_levelset!(geom, phi_init)
    
    save_reinit_plots(geom, "merging_before", output_dir, phi_ref=phi_gt)
    reinitialize!(geom)
    save_reinit_plots(geom, "merging_after", output_dir, phi_ref=phi_gt)
    
    phi_new = current_levelset(geom)
    err = phi_new .- phi_gt
    l2_err = sqrt(sum(err.^2) / length(err))
    @test l2_err < 0.05
end