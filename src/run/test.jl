#using Plots
using PowerDynamics
using DifferentialEquations
#using DiffEqSensitivity
#using LinearAlgebra
using ModelingToolkit


include("C:/Users/liemann/github/PowerDynamicsDevelopment/src/include_costum_nodes_lines_utilities.jl")
include("C:/Users/liemann/github/PowerDynamicsDevelopment/src/grids/GFC_Test_Grid.jl")
include("C:/Users/liemann/github/PowerDynamicsDevelopment/src/sensitivity_analyses/Local_Sensitivity.jl")
################################################################################
#function InitMyGFC()
pg = GFC_Test_Grid()
U,δ,ic0 = PowerFlowClassic(pg, iwamoto = true, max_tol = 1e-7)

pg1 ,ic = InitializeInternalDynamics(pg,ic0)
params = GFC_params()
prob = ODEProblem(rhs(pg1),ic,(0.0,3),params)
pgsol,evr = simGFC(prob)

mtk_normal = GetMTKSystem(pg,(0.0,1.0),params)
mtk_fault = GetMTKSystem(GFC_Test_Grid(y_new = yfault()),(0.0,1.0),params)
mtk = [mtk_normal,mtk_fault]
s = GetTriggCondsGFCTest(mtk_normal)
h = GetStateResFunGFCTest(mtk_normal)
#@time sens = CalcHybridTrajectorySensitivity(mtk,pgsol.dqsol,params,evr,s,h,[],collect(1:15))
#############################################################################
#CalcHybridTrajectorySensitivity(mtk::Vector{ODESystem},sol,p_pre,evr,s,h,u0_sensi,p_sensi)
sol = pgsol.dqsol
p_pre = params
u0_sensi = []
p_sensi = collect(1:15)
#############################################################################
mtk0 = mtk[1] # it is assumed that the first element is the initial system
ic = sol.prob.u0
sensis_u0_pd = u0_sensi
sensis_p_pd = p_sensi
p = p_pre
#############################################################################
@time xx0_k, yx0_k, sym_states,sym_params, A_states, D_states, M, N, O, symp, Δt,len_sens, f, g, J =
   InitTrajectorySensitivity(mtk0, ic, p_pre, u0_sensi, p_sensi);
mtsys = mtk0
fulleqs = equations(mtsys)
sym_states = states(mtsys)
sym_params = parameters(mtsys)
eqs, aeqs, D_states, A_states = GetSymbolicEquationsAndStates(fulleqs, sym_states)
sensis_u0 = sym_states[sensis_u0_pd]
#sensis_p_pd is here a list with indices of the parameters p
sensis_p = sym_params[sensis_p_pd]

#dict from states and parameters with their starting values
symu0 = sym_states .=> ic
symp = sym_params .=> p

Fx, Fy, Gx, Gy = GetSymbolicFactorizedJacobian(eqs, aeqs, D_states, A_states)

len_sens = size(sensis_u0)[1] + size(sensis_p)[1]
Diff_u0 = Differential.(sensis_u0)
Diff_p = Differential.(sensis_p)
Fp = Array{Num}(undef, size(eqs)[1], len_sens)
Gp = Array{Num}(undef, size(aeqs)[1], len_sens)

Fp[:, 1:size(Diff_u0)[1]] .= Num(0)
Gp[:, 1:size(Diff_u0)[1]] .= Num(0)

for (ind, val) in enumerate(Diff_p)
  Fp[:, ind+size(Diff_u0)[1]] =
    Num.(expand_derivatives.(map(val, my_rhs.(eqs))))
  Gp[:, ind+size(Diff_u0)[1]] =
    Num.(expand_derivatives.(map(val, my_rhs.(aeqs))))
end

@parameters Δt
@parameters xx0[1:size(D_states)[1], 1:len_sens] #xx0 are the sensitivities regarding differential states
@parameters yx0[1:size(A_states)[1], 1:len_sens] #yx0 are the sensitivities regarding algebraic states
xx0 = Symbolics.scalarize(xx0)
yx0 = Symbolics.scalarize(yx0)
#end

@time tmp = InitMyGFC()
display(length.(string.(equations(tmp))))
