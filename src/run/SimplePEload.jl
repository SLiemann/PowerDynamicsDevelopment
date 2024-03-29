using DifferentialEquations
using Plots
using SciMLSensitivity
using ForwardDiff
using MATLAB
include("C:/Users/liemann/github/PowerDynamicsDevelopment/src/sensitivity_analyses/Local_Sensitivity.jl")

function simPEload(;u0=[10*sqrt(2),0.0,0.0,1.0],p=[10*sqrt(2),100*pi,5e-3,1e-3,50],dt_max=1e-5)
  function f(dx,x,p,t)
      Ud,w,L,C,R = p
      dx[1] = x[1] - Ud*cos(w*t)
      dx[2] = (x[1] - x[3])/L * x[4]
      dx[3] = (x[2] - x[3]/R)/C
      dx[4] = 0.0
  end
  #dx[4] = x[4] - (x[1]*x[3] - C*dx[2]) 
  M = [0  0  0  0
       0  1. 0  0
       0  0  1. 0
       0  0  0  1.];

  evr = Array{Float64}(undef,0,4)

  s1(u, t, integrator) = u[2] + integrator.u[4] - 1
  function h1(integrator) 
    integrator.u[4] = 0.0
    evr = vcat(evr, [integrator.t 1 1 1])
  end
  cb1 = ContinuousCallback(s1,nothing, affect_neg! = h1)

  function s2(u, t, integrator)
      Ud,w,L,C,R = integrator.p
      -((Ud*cos(w*t) - integrator.u[3]) + integrator.u[4] - 1e-2)
  end
  function h2(integrator) 
    integrator.u[4] = 1.0
    evr = vcat(evr, [integrator.t 1 2 2])
  end
  cb2 = ContinuousCallback(s2, nothing,affect_neg! =h2)

  ode_fun = ODEFunction(f,mass_matrix = M)
  prob = ODEProblem(ode_fun, u0, (0.0, 0.16), p)
  #_prob = remake(prob,p=theta)
  #prob = ODEForwardSensitivityProblem(ode_fun, u0, (0.0, 0.03), p,ForwardDiffSensitivity();)
  #sol = Array(solve(_prob, Rodas4(), callback = CallbackSet(cb1,cb2),dtmax=1e-4))' #,abstol=1e-10, 
  sol = solve(prob, Rodas4P(), callback = CallbackSet(cb1,cb2),dtmax=dt_max)#,abstol=1e-6,reltol=1e-6)
  return sol, evr
end

sol,evr = simPEload();
plot(sol.t,sol[1,:])
plot!(sol.t,sol[3,:])
plot!(sol.t,sol[2,:])
plot!(sol.t,sol[4,:])

#p = [sqrt(2),100*pi,5e-3,1e-3,50]
#res = SortFDResults(ForwardDiff.jacobian(simPEload,p),4);
#plot!(res[2][1,:])
#plot!(hybrid_sen[5][1,:])
#plot!(sol[3,:])


######
mtk = modelingtoolkitize(sol.prob)
sym_states = states(mtk)
symp = parameters(mtk)

hs = Vector{Vector{Num}}(undef,0)
ident = Num.(vcat(sym_states[2:4],symp))
ident[3] = 0.0
push!(hs,ident)
ident =  Num.(vcat(sym_states[2:4],symp))
ident[3] = 1.0
push!(hs,ident)

@parameters t
s = Vector{Num}(undef,0)
push!(s,sym_states[2])
push!(s,-((symp[1]*cos(symp[2]*t) - sym_states[3])+sym_states[4]-1e-2))

hybrid_sen,Δτ = CalcHybridTrajectorySensitivity([mtk],sol,evr,s,hs);


ProfileView.@profview for i=1:5 CalcHybridTrajectorySensitivity([mtk],sol,evr,s,hs) end

sensi_ind = 7
Δp = -0.3*1e-3
sol_appr = ApproximatedTrajectory(sol,hybrid_sen[sensi_ind],Δp);
sol_perc,evr_per = simPEload(p=[10*sqrt(2),100*pi,5e-3,1e-3+Δp,50.0],dt_max=1e-5);
sol_refin = TrajectoryRefinement([mtk],sol,evr,hybrid_sen,Δτ,sensi_ind,Δp);

x=3
xl=(0.075,0.09)
yl=(-0.01,0.2)
plot(sol.t,sol[x,:])
plot!(sol_perc.t,sol_perc[x,:])
plot!(sol.t,sol_appr[x,:])
plot!(sol.t,sol_refin[x,:])

plot!(sol.t,sol_refin[x,:] ,xlims=xl,ylims=yl)

begin
  sol_perc,evr_per = simPEload(p=[10*sqrt(2),100*pi,5e-3,1e-3+Δp,50.0],dt_max=1e-5);
  plot(sol.t,sol[x,:])
  plot!(sol_perc.t,sol_perc[x,:])
  display(plot!(sol.t,sol_appr[x,:]))
end

####Plotten der Sensis für Diss
plot(sol.t,hybrid_sen[7][3,:].*1e-3) #von C auf V
plot!(sol.t,hybrid_sen[6][3,:].*5e-3) #von L auf V

plot(sol.t,hybrid_sen[7][2,:].*1e-3) #von C auf I
plot!(sol.t,hybrid_sen[6][2,:].*5e-3) #von L auf I


###Write to MATLAB
write_matfile("C:\\Users\\liemann\\dissertation\\Dissertation\\document\\images\\332_Example_Trajectory_Sensitivity\\sensis_example.mat"; t=sol.t,tperc=sol_perc.t,tperi=sol_peri.t,sol=sol[:,:],sol_perc=sol_perc[:,:],sol_peri=sol_peri[:,:], vc=hybrid_sen[7][3,:], il=hybrid_sen[6][2,:]) 

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



ind1 = vcat(2,setdiff(indexin(evr[:,1],sol.t),[nothing]).+2)
ind2 = vcat(setdiff(indexin(evr[:,1],sol.t),[nothing]),length(sol.t))
ind_sol = Int.(hcat(ind1,ind2))

ind_sol = vcat(2,setdiff(indexin(evr[:,1],sol.t),[nothing]),length(sol.t))


tmp = deepcopy(hybrid_sen)
len_params = 5
len_d0 = 3
D0_indices = [1,2,4]
A_indices = [3]
len_states =  length(states([mtk][1]))
tmp2 = zeros(Float64,len_states,length(sol))
for k in 1:length(sensis)
  tmp[D0_indices,:] = sensis[k][1:len_d0,:]
  tmp[A_indices,:] = sensis[k][len_d0+len_params+1:end,:]
  sensis[k] = deepcopy(tmp)
end

###Write to MATLAB
write_matfile("ts_example_res.mat"; t=sol.t,sol = sol[:,:]) 
