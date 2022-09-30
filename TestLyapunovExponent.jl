using DynamicalSystems, Plots
#using CairoMakie

ds = Systems.lorenz(σ=16,ρ=45.92,β=4)
# create a timeseries of 1 dimension
Δt = 0.01
N = 5000
time = 0.0:Δt:N*Δt
x = trajectory(ds, N*Δt; Δt)[:, 1]

ks2 = 0:4:200

#fig = Figure()
#ax = Axis(fig[1,1]; xlabel="k (0.05×t)", ylabel="E - E(0)")
ntype = NeighborNumber(5) #5 nearest neighbors of each state

begin
    #x = sin.(0:1e-2:4*pi)
    f = plot()
    J = getConstructionDelay(x)
    for kend in 1:100:2000
        ks2 = 0:1:kend

        
        for d in 8, τ in J
            r = embed(x, d, τ)

            E2 = lyapunov_from_data(r, ks2; ntype)
            λ2 = linear_region(ks2 .* Δt, E2)[2]
            #display(E2[end])
            plot!(ks2, E2; label = "kend=$(kend), λ=$(round(λ2, digits = 3))",legend=false) #, d=$(d), τ=$(τ)
            #plot([0,200],[-1])
        end
    end
    f
end

using StatsBase


function getConstructionDelay(x)
    acf = autocor(x)
    #bar(acf)
    #plot!([0,length(acf)],[1-1/exp(1),1-1/exp(1)])
    return findfirst(acf .< (1-1/exp(1)))
end

