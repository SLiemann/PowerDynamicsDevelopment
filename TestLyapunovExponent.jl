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
    #for kend in 1:100:2000
        ks2 = 0:4:Int64(floor(length(x)/10))

        for d in [1], N in 1:30#, τ in 1:3:9
            r = embed(x, d, J)
            ntype =  NeighborNumber(N)
            E2 = lyapunov_from_data(r, ks2;ntype)
            λ2 = linear_region(ks2 .* Δt, E2)[2]
            #display(E2[end])
            if ~isnan(λ2)
                #plot!(ks2, E2; label = "N=$(N),d=$(d), λ=$(round(λ2, digits = 3))",legend=:topleft) #, d=$(d), τ=$(τ)
                #plot([0,200],[-1])
                bar!([N],[λ2],legend=false)
            end
        end
    #end
    f
end

using StatsBase


function getConstructionDelay(x;demean=true)
    acf = autocor(x,collect(1:length(x)-1),demean=demean)
    #bar(acf)
    #plot!([0,length(acf)],[1-1/exp(1),1-1/exp(1)])
    J  = findfirst(acf .< (1-1/exp(1)))
    if isnothing(J)
        @warn "Autocorrelation higher then 0.63212... \n setting J to 1"
        return 1
    else
        return J
    end
end

using CSV, DataFrames
 test = DataFrame(CSV.File("C:\\Users\\liemann\\dissertation\\Dissertation\\document\\images\\231_LVRT curves and cases\\unstable.csv"; header=false, delim=';', type=Float64))

 x = Float64.(test[3:end,2])
 getConstructionDelay(x)