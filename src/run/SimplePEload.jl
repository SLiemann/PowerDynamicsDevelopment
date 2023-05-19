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

hs = Num.([D_states D_states])
hs[3,1] = Num(0.0)
hs[3,2] = Num(1.0)

s = [D_states[2], D_states[3]-1.0]

hx = Array{Array{Num}}(undef,size(hs)[2],1)
hy = similar(hx)
sx = Array{Array{Num}}(undef,length(s),1)
sy = similar(sx)


#HIER WEITER MACHE!!!!!
for i=eachindex(hs)
    hx[i] = GetJacobian(hs[:,i],D_states)
    hy[i] = GetJacobian(hs[:,i],A_states)
end 
for i=eachindex(s)
    sx[i] = GetJacobian([s[i]],D_states)
    sy[i] = GetJacobian([s[i]],A_states)
end