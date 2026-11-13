using Test
using Main.EvolvingDomains
using Main.EvolvingDomains.Geometric.Reinitialization: reinitialize!

function _sphere_reinitialization_error(n)
    axis=collect(range(0.0,1.0,length=n)); grid=EvolvingDomains.Grid((axis,axis,axis))
    phi=Array{Float64}(undef,n,n,n)
    for k in 1:n,j in 1:n,i in 1:n
        phi[i,j,k]=sqrt((axis[i]-.5)^2+(axis[j]-.5)^2+(axis[k]-.5)^2)-.25
    end
    cells=n-1; states=Vector{Int8}(undef,cells^3); c=0
    for k in 1:cells,j in 1:cells,i in 1:cells
        vals=(phi[i,j,k],phi[i+1,j,k],phi[i,j+1,k],phi[i+1,j+1,k],
              phi[i,j,k+1],phi[i+1,j,k+1],phi[i,j+1,k+1],phi[i+1,j+1,k+1])
        lo,hi=extrema(vals); c+=1
        states[c]=lo<0<=hi ? 0 : (hi<0 ? -1 : 1)
    end
    facets=fill((Int8(1),Int8(1),Int8(1),Int8(1),Int8(1),Int8(1)),cells^3)
    cut=EvolvingDomains.Cut(grid,:sphere;bg_states=states,cell_facet_states=facets)
    field=EvolvingDomains.Field(EvolvingDomains.Discretization(Dict(:sphere=>grid)),Dict(:sphere=>phi))
    reinitialize!(field,cut)
    i=round(Int,.875*(n-1))+1; j=round(Int,.5*(n-1))+1; k=j
    exact=sqrt((axis[i]-.5)^2+(axis[j]-.5)^2+(axis[k]-.5)^2)-.25
    abs(phi[i,j,k]-exact)
end

@testset "3D field/cut reinitialization" begin
    grid=EvolvingDomains.Grid((collect(range(0.0,1.0,length=9)),collect(range(0.0,1.0,length=7)),collect(range(0.0,1.0,length=5))))
    d=EvolvingDomains.Discretization(Dict(:fine=>grid))
    phi=Array{Float64}(undef,9,7,5)
    for k in axes(phi,3),j in axes(phi,2),i in axes(phi,1); phi[i,j,k]=grid.axes[1][i]-0.53; end
    states=Int8[]
    for k in 1:4,j in 1:6,i in 1:8
        push!(states, i==5 ? 0 : (i<5 ? -1 : 1))
    end
    facets=fill((Int8(1),Int8(1),Int8(1),Int8(1),Int8(1),Int8(1)),8*6*4)
    cut=EvolvingDomains.Cut(grid,:fine;bg_states=states,cell_facet_states=facets)
    field=EvolvingDomains.Field(d,Dict(:fine=>phi))
    before=copy(phi)
    reinitialize!(field,cut)
    @test all(isfinite,phi)
    @test sign.(phi)==sign.(before)
    # A plane's signed distance is exact away from the frozen cut-node band.
    @test maximum(abs.(abs.(phi[8,:,:]) .- 0.345)) < 0.08
    @test maximum(abs.(phi[5,:,:] .- before[5,:,:])) < 1e-12

    # No CUT cell means a well-formed level set is unchanged.
    unchanged=copy(phi)
    nocut=EvolvingDomains.Cut(grid,:fine;bg_states=fill(Int8(-1),8*6*4),cell_facet_states=facets)
    @test_logs (:warn,r"single-signed") reinitialize!(field,nocut)
    @test phi==unchanged
    @test_throws ArgumentError reinitialize!(field,EvolvingDomains.Cut(grid,:other;bg_states=states,cell_facet_states=facets))
    shortfield=EvolvingDomains.Field(d,Dict(:fine=>zeros(8,7,5)))
    @test_throws DimensionMismatch reinitialize!(shortfield,cut)
    floatfield=EvolvingDomains.Field(d,Dict(:fine=>zeros(Float32,9,7,5)))
    @test_throws ArgumentError reinitialize!(floatfield,cut)
    othergrid=EvolvingDomains.Grid((collect(range(0.0,1.0,length=9)),collect(range(0.0,1.0,length=7)),[0.0,.2,.45,.7,1.0]))
    otherd=EvolvingDomains.Discretization(Dict(:fine=>othergrid))
    otherfield=EvolvingDomains.Field(otherd,Dict(:fine=>zeros(9,7,5)))
    othercut=EvolvingDomains.Cut(othergrid,:fine;bg_states=states,cell_facet_states=facets)
    @test_throws ArgumentError reinitialize!(otherfield,othercut)
    phi[1]=NaN
    @test_throws ArgumentError reinitialize!(field,cut)

    coarse=_sphere_reinitialization_error(9)
    fine=_sphere_reinitialization_error(17)
    @test fine <= coarse + 1e-12
    @test fine < 0.04
end
