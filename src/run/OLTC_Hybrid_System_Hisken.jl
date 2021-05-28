using PowerDynamics
using DifferentialEquations
using Plots
using ModelingToolkit

#data1 = sqrt.((sensi[2][7, :] .^ 2) .+ (sensi[2][8, :] .^ 2))
#plot(pgsol.dqsol.t[1:end-1], sensi[2][3, :])

begin
    include("C:/Users/liemann/github/PowerDynamicsDevelopment/src/grids/OLTC_Hybrid_Sensis.jl")
    include("C:/Users/liemann/github/PowerDynamicsDevelopment/src/sensitivity_analyses\\Local_Sensitivity.jl")
end


mtk = GetMTKOLTCSystem()
pg, ic = GetInitializedOLTCHisken()
tspan = tspanOLTCHisken()
params = GetParametersOLTCHisken()
pgsol,evr = SimulateOLTCHIsken()
sol = pgsol.dqsol
#sensi = TimeDomainSensitivies(mtk, tspan, ic, params, [1], [3], pgsol.dqsol)
function CalcHybridTrajectorySensitivity(mtk,sol,evr,u0_sensi,p_sensi)
    sensis = Vector{Array{Float64}}(undef, len_sens)
    for i = 1:length(sensi)
      sensis[i] = Array{Float64}(
        undef,
        size(D_states)[1] + size(A_states)[1],
        size(sol)[2] - 1,
      )
    ic = sol.prob.uo
    p  = evr[1,2:end-2]
    ind_sol = vcat(1,setdiff(indexin(evr[:,1],sol.t),[nothing]),indexin(sol.t[end],sol.t))
    xx0_k, yx0_k, state, A_states, D_states, M, N, O, symp, Δt,len_sens, J =
        InitTrajectorySensitivity(mtk, ic, params, u0_sensi, p_sensi, sol)
    fx,fy,gx,gy = J

    for i = 1:length(ind_sol)
        sol_part = sol[ind_sol[i]:ind_sol[i+1]]
        sensi_part = ContinuousSensitivity(
                                    sol,
                                    xx0_k,
                                    yx0_k,
                                    state,
                                    A_states,
                                    D_states,
                                    M,
                                    N,
                                    O,
                                    symp,
                                    Δt,
                                    len_sens,
                                )
        for j = 1:length(sensi_part)
            sensis[j][:,ind_sol[i]:ind_sol[i+1]-1] = sensi_part[j]
        end
        if ind_sol[i+1] ≠ ind_sol[end]
            # Hier weitermachen mit dem Sprung!
        else
            return sensis
        end

    end
end

function test()
    a = 1
    for i = 1:3
        display(a)
        if a > -1
            a += 1
        end
    end
end
test()
