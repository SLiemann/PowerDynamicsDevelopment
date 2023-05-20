using DifferentialEquations
using Plots
using SciMLSensitivity

function f(dx,x,p,t)
    Ud,w,L,C,R = p
    dx[1] = (Ud*cos(w*t) - x[2])/L * x[3]
    dx[2] = (x[1] - x[2]/R)/C
    dx[3] = 0.0
    dx[4] = x[4] - (x[1]*x[3] - C*dx[2]) 
end
M = [1. 0  0  0
     0  1. 0  0
     0  0  1. 0
     0  0  0  0];

s1(u, t, integrator) = u[1]
h1(integrator) = integrator.u[3] = 0.0
cb1 = ContinuousCallback(s1, h1)

function s2(u, t, integrator)
    Ud,w,L,C,R = integrator.p
    (Ud*cos(w*t) - integrator.u[2]) > 0.0
end
h2(integrator) = integrator.u[3] = 1.0
cb2 = DiscreteCallback(s2, h2)

u0 = [0.0,0.0,1.0,0.0]
ic = u0
p = [sqrt(2),100*pi,5e-3,1e-3,50]
ode_fun = ODEFunction(f,mass_matrix = M)
prob = ODEProblem(ode_fun, u0, (0.0, 0.25), p)
sol = solve(prob, Rodas4(), callback = CallbackSet(cb1,cb2),dtmax=1e-3)
plot(sol,idxs=(4))

###################################
sens_prob = ODEForwardSensitivityProblem(f, u0, (0.0, 0.06), p,callback = CallbackSet(cb1,cb2),sensealg=ForwardDiffSensitivity(;chunk_size = 0,convert_tspan =true))
sens_prob = ODEForwardSensitivityProblem(f, u0, (0.0, 0.06), p,callback = CallbackSet(cb1,cb2),sensealg=ReverseDiffAdjoint())
sol = solve(sens_prob, Tsit5(), callback = CallbackSet(cb1,cb2),dtmax=1e-3);
x, dp = extract_local_sensitivities(sol);
da = dp[3];
plot(sol.t, da', lw = 3)
#################################
include("C:/Users/liemann/github/PowerDynamicsDevelopment/src/sensitivity_analyses/Local_Sensitivity.jl")
mtsys = modelingtoolkitize(prob)

fulleqs = equations(mtsys);
sym_states = states(mtsys);
sym_params = parameters(mtsys);
eqs, aeqs, D_states, A_states = GetSymbolicEquationsAndStates(fulleqs, sym_states);
D_states = vcat(D_states,sym_params);

hs = Num.([D_states D_states])
hs[3,1] = Num(0.0)
hs[3,2] = Num(1.0)

s = Num.([D_states[2], D_states[3]-1.0])

mtk0 = mtsys
ic = sol.prob.u0
xx0_k, yx0_k, sym_states,sym_params, A_states, D_states, M, N, Δt,len_sens, f1, g, J =
    InitTrajectorySensitivity(mtk0, ic, p);
xx0 = [i[1] for i in xx0_k];
yx0 = [i[1] for i in yx0_k];
f_all,g_all,J_all,M_all,N_all =  GetEqsJacobianSensMatrices([mtk0],xx0,yx0); #anpassen wenn nicht alle Variablen gerechnet werdne sollen

Fx_all, Fy_all, Gx_all, Gy_all = J_all;
hx,hy,sx,sy = CalcTriggerAndStateResetJacobians(s,hs,D_states,A_states);

sensis = Vector{Array{Float64}}(undef, len_sens)
for i = 1:length(sensis)
  sensis[i] = Array{Float64}(
    undef,
    size(D_states)[1] + size(A_states)[1],
    size(sol)[2] - 1,
  )
end
ind_sol = vcat(1,setdiff(indexin(evr[:,1],sol.t),[nothing]),length(sol.t))