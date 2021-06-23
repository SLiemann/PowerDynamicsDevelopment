using Plots
using DifferentialEquations
using ModelingToolkit


mtk_normal = GetMTKLTVSSystem(pg_state = "gfc_normal")
mtk_fault = GetMTKLTVSSystem(pg_state = "gfc_fault")
mtk_postfault = GetMTKLTVSSystem(pg_state = "gfc_postfault")
mtk = [mtk_normal; mtk_fault; mtk_postfault]



s = GetTriggCondsLTVS(mtk_normal)
h = GetStateResFunLTVS(mtk_normal)
p_pre = GFC_LTVS_params()

#@time toll = CalcHybridTrajectorySensitivity(mtk,pgsol.dqsol,p_pre,evr,s,h,[1],[13])
u0_sensi = [1]
p_sensi = [13]
mtk0 = mtk[1] # it is assumed that the first element is the initial system

ic = pgsol.dqsol.prob.u0
xx0_k, yx0_k, sym_states,sym_params, A_states, D_states, M, N, O, symp, Δt,len_sens, f, g, J =
    InitTrajectorySensitivity(mtk0, ic, p_pre, u0_sensi, p_sensi)
xx0 = [i[1] for i in xx0_k]
yx0 = [i[1] for i in yx0_k]
f_all,g_all,J_all,M_all,N_all,O_all = GetEqsJacobianSensMatrices(mtk,xx0,yx0,u0_sensi,p_sensi)
Fx_all, Fy_all, Gx_all, Gy_all = J_all
hx,hy,sx,sy = CalcTriggerAndStateResetJacobians(mtk0,s,h)
sensis = Vector{Array{Float64}}(undef, len_sens)
for i = 1:length(sensis)
  sensis[i] = Array{Float64}(
    undef,
    size(D_states)[1] + size(A_states)[1],
    size(pgsol.dqsol)[2] - 1,
  )
end
ind_sol = vcat(1,setdiff(indexin(evr[:,1],pgsol.dqsol.t).+1,[nothing]),length(pgsol.dqsol.t))
#ind_sol = [1]
#for i in evr[:,1] # DifferentialEquations.jl has multiple time points
#    ind_sol = vcat(ind_sol,findall(x->x==i,sol.t)[end])
#end
#ind_sol = vcat(ind_sol,length(sol.t))

Fx_pre = Fx_all[1]
Fy_pre = Fy_all[1]
Gx_pre = Gx_all[1]
Gy_pre = Gy_all[1]
J_pre = Fx_pre, Fy_pre, Gx_pre, Gy_pre

f_pre = f_all[1]
g_pre = g_all[1]
