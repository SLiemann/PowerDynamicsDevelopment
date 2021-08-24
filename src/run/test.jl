using Plots
using PowerDynamics
using DifferentialEquations
using DiffEqSensitivity
using LinearAlgebra
using ModelingToolkit

function f(du,u,p,t)
  du[1] = dx = p[1]*u[1] - p[2]*u[1]*u[2]
  du[2] = dy = -p[3]*u[2] + u[1]*u[2]
  du[3] = dz = u[1] + u[2]
end

ode_f = ODEFunction(f,mass_matrix = Diagonal([1.;1.;1.0]))
p = [1.5,1.0,3.0]
prob = ODEForwardSensitivityProblem(ode_f,[1.0;1.0;0.0],(0.0,10.0),p,sensealg = ForwardDiffSensitivity();)
sol = solve(prob,DP8())

################################################################################
################################################################################
################################################################################
################################################################################
begin
    include("C:/Users/liemann/github/PowerDynamicsDevelopment/src/include_costum_nodes_lines_utilities.jl")
    include("C:/Users/liemann/github/PowerDynamicsDevelopment/src/grids/GFC_Test_Grid.jl")
    include("C:/Users/liemann/github/PowerDynamicsDevelopment/src/sensitivity_analyses/Local_Sensitivity.jl")

    pg = GFC_Test_Grid()
    U,δ,ic0 = PowerFlowClassic(pg, iwamoto = true, max_tol = 1e-7)

    pg1 ,ic = InitializeInternalDynamics(pg,ic0)
    params = GFC_params()
    prob = ODEProblem(rhs(pg1),ic,(0.0,3),params)
    pgsol,evr = simGFC(prob)
end
plot(pgsol,collect(keys(pg.nodes))[2:3],:v)
plot(pgsol,["bus3"],:i_abs, label = "Iabs")

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
p_sensi = collect(1:16)
#############################################################################
mtk0 = mtk[1] # it is assumed that the first element is the initial system
ic = sol.prob.u0
#############################################################################
#xx0_k, yx0_k, sym_states,sym_params, A_states, D_states, M, N, O, symp, Δt,len_sens, f, g, J =
#    InitTrajectorySensitivity(mtk0, ic, p_pre, u0_sensi, p_sensi)
#InitTrajectorySensitivity(
#  mtsys::ODESystem,
#  ic::Array{Float64,1},
#  p::Array{Float64,1},
#  sensis_u0_pd::Array{Any,1}, #array of indices
#  sensis_p_pd::Array{Int64,1}, #array of indices
#)
mtsys = mtk0
p = p_pre
sensis_u0_pd = u0_sensi
sensis_p_pd = p_sensi
##############################
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

@time Fx, Fy, Gx, Gy = GetSymbolicFactorizedJacobian(eqs, aeqs, D_states, A_states)

len_sens = size(sensis_u0)[1] + size(sensis_p)[1]
Diff_u0 = Differential.(sensis_u0)
Diff_p = Differential.(sensis_p)
Fp = Array{Num}(undef, size(eqs)[1], len_sens)
Gp = Array{Num}(undef, size(aeqs)[1], len_sens)

Fp[:, 1:size(Diff_u0)[1]] .= Num(0)
Gp[:, 1:size(Diff_u0)[1]] .= Num(0)

@time for (ind, val) in enumerate(Diff_p)
  Fp[:, ind+size(Diff_u0)[1]] =
    #Num.(expand_derivatives.(map(val, my_rhs.(eqs))))
  Gp[:, ind+size(Diff_u0)[1]] =
    Num.(expand_derivatives.(map(val, my_rhs.(aeqs))))
end

@parameters Δt
@parameters xx0[1:size(D_states)[1], 1:len_sens] #xx0 are the sensitivities regarding differential states
@parameters yx0[1:size(A_states)[1], 1:len_sens] #yx0 are the sensitivities regarding algebraic states
@time M,N,O = TrajectorySensitivityMatrices([Fx, Fy, Gx, Gy],Fp,Gp,xx0,yx0,aeqs,A_states,len_sens);

#Initialisierung: xx0 enthält die Sensis für x0 und p bezüglich Differentialzustände
xx0_k = xx0 .=> 0.0
xx0_f = zeros(size(xx0)[1], len_sens)
ind = setdiff(indexin(sensis_u0, D_states),[nothing])
for i = 1:length(ind)
  xx0_k[i, ind[i]] = xx0_k[i, ind[i]][1] => 1.0
  xx0_f[i, ind[i]] = 1.0
end
# Bei den Sensis für y werden zuerst die dy/x0 Sensi initialisiert
@time Gy_float = Substitute.(Gy, [symu0; symp])
for i = 1:16
  @time tmp = substitute.(Gp[:,i], ([symu0; symp],))
end
tmp = 2*tmp
@time tmp = Symbolics.value.(tmp)
@time Gp_float = Substitute(Gp, [symu0; symp])

yx0_t0 = -inv(Gy_float) * (Gx * xx0_f[:, 1:size(sensis_u0)[1]])
yp_t0 =
  -inv(Gy_float) *
  (Gp * vcat(zeros(size(sensis_u0)[1], size(sensis_p)[1]), I))
@time yx0_k = yx0 .=> Substitute([yx0_t0 yp_t0], [symu0; symp])

eqs = Num.(my_rhs.(eqs))
aeqs = Num.(my_rhs.(aeqs))

b = false
a =
  b ? 1 : 0
