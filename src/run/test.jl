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

u0_sens = u0_sensi
p_sens = p_sensi

f = Vector{Array{Num,1}}(undef,length(mtk))
g = similar(f)

Fx = Vector{Array{Num,2}}(undef,length(mtk))
Fy = similar(Fx)
Gx = similar(Fx)
Gy = similar(Fx)
Fp = similar(Fx)
Gp = similar(Fx)

M = similar(Fx)
N = similar(Fx)
O = similar(Fx)

len_sens = size(u0_sens)[1] + size(p_sens)[1]
sensis_u0 = states(mtk[1])[u0_sens]
sensis_p = parameters(mtk[1])[u0_sens]
Diff_u0 = Differential.(u0_sens)
Diff_p = Differential.(p_sens)

for (ind,val) in enumerate(mtk)
    fulleqs = equations(val)
    symstates = states(val)
    f_tmp, g_tmp, x, y = GetSymbolicEquationsAndStates(fulleqs, symstates)
    Fx[ind],Fy[ind],Gx[ind],Gy[ind] = GetSymbolicFactorizedJacobian(f_tmp, g_tmp, x, y)

    Fp[ind] = Num(0).*zeros(length(x),len_sens) #init with zeros
    Gp[ind] = Num(0).*zeros(length(y),len_sens)
    for (ind2, val2) in enumerate(Diff_p)
      Fp[ind][:, ind2+size(Diff_u0)[1]] =
        Num.(expand_derivatives.(map(val2, my_rhs.(f_tmp))))
      Gp[ind][:, ind2+size(Diff_u0)[1]] =
        Num.(expand_derivatives.(map(val2, my_rhs.(g_tmp))))
    end

    M[ind],N[ind],O[ind] = TrajectorySensitivityMatrices([Fx[ind],Fy[ind],Gx[ind],Gy[ind]],Fp[ind],Gp[ind],xx0,yx0,g_tmp,y,len_sens)
    f[ind] = my_rhs.(f_tmp)
    g[ind] = my_rhs.(g_tmp)
end
