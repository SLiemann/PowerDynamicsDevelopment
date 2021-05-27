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
#sensi = TimeDomainSensitivies(mtk, tspan, ic, params, [1], [3], pgsol.dqsol)

xx0_k, yx0_k, state, A_states, D_states, M, N, O, symp, Δt,len_sens =
    InitTrajectorySensitivity(mtk, tspan, ic, params, [1], [3], pgsol.dqsol)

pgsol.dqsol.t

sol = deepcopy(pgsol.dqsol)
ind_sol = vcat(1,indexin(evr[:,1],sol.t),indexin(sol.t[end],sol.t))

for i = 1:3

    sensi = ContinuousSensitivity(
                                sensi,
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
end
