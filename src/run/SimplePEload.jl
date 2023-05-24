using DifferentialEquations
using Plots
using SciMLSensitivity
using ForwardDiff
include("C:/Users/liemann/github/PowerDynamicsDevelopment/src/sensitivity_analyses/Local_Sensitivity.jl")

function simPEload(theta)
  function f(dx,x,p,t)
      Ud,w,L,C,R = p
      dx[1] = (Ud*cos(w*t) - x[2])/L * x[3]
      dx[2] = (x[1] - x[2]/R)/C
      dx[3] = 0.0
      dx[4] = x[4] - Ud*cos(w*t)
  end
  #dx[4] = x[4] - (x[1]*x[3] - C*dx[2]) 
  M = [1. 0  0  0
       0  1. 0  0
       0  0  1. 0
       0  0  0  0];

  evr = Array{Float64}(undef,0,4)

  s1(u, t, integrator) = u[1] #+ integrator.u[3] - 1
  function h1(integrator) 
    integrator.u[3] = 0.0
    evr = vcat(evr, [integrator.t 1 1 1])
  end
  cb1 = ContinuousCallback(s1,nothing, affect_neg! = h1)

  function s2(u, t, integrator)
      Ud,w,L,C,R = integrator.p
      -((Ud*cos(w*t) - integrator.u[2]) + integrator.u[3] - 1e-2)
  end
  function h2(integrator) 
    integrator.u[3] = 1.0
    evr = vcat(evr, [integrator.t 1 2 2])
  end
  cb2 = ContinuousCallback(s2, nothing,affect_neg! =h2)

  u0 = [0.0,0.0,1.0,sqrt(2)]
  p = [sqrt(2),100*pi,5e-3,1e-3,50]
  ode_fun = ODEFunction(f,mass_matrix = M)
  prob = ODEProblem(ode_fun, u0, (0.0, 0.025), p)
  _prob = remake(prob,p=theta)
  prob = ODEForwardSensitivityProblem(ode_fun, u0, (0.0, 0.03), p,ForwardDiffSensitivity();)
  sol = Array(solve(_prob, Rodas4(), callback = CallbackSet(cb1,cb2),dtmax=1e-4))' #,abstol=1e-10, 
  #sol = solve(_prob, Rodas4(), callback = CallbackSet(cb1,cb2),dtmax=1e-4)
  #return sol, evr
end

p = [sqrt(2),100*pi,5e-3,1e-3,50]
res = SortFDResults(ForwardDiff.jacobian(simPEload,p),4);
plot(res[2][1,:].*100)
plot!(hybrid_sen[5][1,:].*100)
plot!(sol[3,:])


sol,evr = simPEload(p);
plot(sol.t,sol[1,:])
plot!(sol.t,sol[2,:])
plot!(sol.t,sol[3,:])
plot!(sol.t,sol[4,:])
plot(sol)


######
mtk = modelingtoolkitize(sol.prob)
sym_states = states(mtk)
symp = parameters(mtk)

hs = Vector{Vector{Num}}(undef,0)
ident = Num.(vcat(sym_states[1:3],symp))
ident[3] = 0.0
push!(hs,ident)
ident =  Num.(vcat(sym_states[1:3],symp))
ident[3] = 1.0
push!(hs,ident)

@parameters t
s = Vector{Num}(undef,0)
push!(s,sym_states[1])
push!(s,-((symp[1]*cos(symp[2]*t) - sym_states[2])+sym_states[3]-1e-2))

hybrid_sen = CalcHybridTrajectorySensitivity([mtk],sol,evr,s,hs);
plot(sol.t[1:end-1],hybrid_sen[4][9,:])

###################################
sens_prob = ODEForwardSensitivityProblem(f, u0, (0.0, 0.06), p,callback = CallbackSet(cb1,cb2),sensealg=ForwardDiffSensitivity(;chunk_size = 0,convert_tspan =true))
sens_prob = ODEForwardSensitivityProblem(f, u0, (0.0, 0.06), p,callback = CallbackSet(cb1,cb2),sensealg=ReverseDiffAdjoint())
sol = solve(sens_prob, Tsit5(), callback = CallbackSet(cb1,cb2),dtmax=1e-3);
x, dp = extract_local_sensitivities(sol);
da = dp[3];
plot(sol.t, da', lw = 3)


#################################
include("C:/Users/liemann/github/PowerDynamicsDevelopment/src/sensitivity_analyses/Local_Sensitivity.jl")

function lv(du, u, p, t)
  du[1] = dx = p[1] * u[1] - p[2] * u[1] * u[2]
  du[2] = dy = -p[3] * u[2] + u[1] * u[2]
  du[3] = dz = u[3] - (u[1]+u[2])
end
Mlv = [1. 0  0 
0  1. 0  
0  0  0];

p = [1.5, 1.0, 3.0];
ode_fun = ODEFunction(lv,mass_matrix=Mlv);
prob = ODEProblem(ode_fun,[1.0; 1.0;2.0], (0.0, 10.0), p)
sol = solve(prob,Rodas4(),dtmax=1e-2);

mtk = modelingtoolkitize(prob)
sen_sol = CalcContinuousSensitivity([mtk],sol);
plot(sol.t[1:end-1],sen_sol[5][[1,2,6],:]')


prob2 = ODEForwardSensitivityProblem(ode_fun, [1.0; 1.0;2.0], (0.0, 10.0), p,ForwardDiffSensitivity();)
sol2 = solve(prob2, Rodas4(),dtmax=1e-2);
x, dp = extract_local_sensitivities(sol2);
da = dp[3]

plot(sol2.t, da') #,size=(2000,2000)