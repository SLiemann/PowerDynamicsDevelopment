using PowerDynamics: SlackAlgebraic, VoltageDependentLoad, SlackAlgebraicParam
using PowerDynamics
import PowerDynamics: PiModel
using OrderedCollections: OrderedDict
using DifferentialEquations
using DiffEqSensitivity
using Plots

function fiip(du,u,p,t)
  du[1] = p[1]*u[1] - p[2]*u[1]*u[2]
  du[2] =-p[3]*u[2] + u[1]*u[2]
  out[3] = u[3]*u[1] + u[2] -0
end

function fiip2(du,u,p,t)
  du[1] = p[1]*u[1] - p[2]*u[1]*u[2]
  du[2] =-p[3]*u[2] + u[1]*u[2]
  du[3] = u[1] + u[2] - u[3]
end

begin
  include("operationpoint/Local_Sensitivity_test.jl")
  p  = [1.5,1.0,3.0]
  u0 = [1.0,1.0,2.0]
  tspan = (0.0,20.0)
  odefun = ODEFunction(fiip2,mass_matrix = Diagonal([1,1,0]))
  ODEProb = ODEProblem(odefun,u0,tspan,p)
  sol = solve(ODEProb, Rodas5(),dt = 1e-3, adaptive = false)
  sensis = TimeDomainSensitivies(ODEProb,tspan,u0,p,[1],[1,2],sol)

  sol_u1 = Vector{Float64}(undef,length(sol.u))
  for (index,value) in enumerate(sol.u)
    sol_u1[index] = value[1]
  end

  Δu0 = [0.0,0.0]
  Δp  =  [0.03,0.0,0.0]
  newprob = ODEProblem(fiip,u0+Δu0,tspan,p+Δp)
  pert_sol = solve(newprob, Rodas4(),dt = 1e-3, adaptive = false)

  sol_pert_u1 = Vector{Float64}(undef,length(pert_sol.u))
  for (index,value) in enumerate(pert_sol.u)
    sol_pert_u1[index] = value[1]
  end
end
begin
  plot(sol.t,sol_u1)
  plot!(pert_sol.t,sol_pert_u1)
  plot!(sol.t[1:end-1],sol_u1[1:end-1]+sensis[2]'[:,1]*Δp[1])
  xlims!((10.0,20))
end

begin
  ref_prob = ODEForwardSensitivityProblem(fiip,u0,tspan,p)
  ref_sens = solve(ref_prob,DP8(),dt = 1e-3, adaptive = false)
  x,dp = extract_local_sensitivities(ref_sens)
end
begin#Vergleich der Sensitivitäten
  plot(ref_sens.t,dp[1]')
  plot!(sol.t[1:end-1],sensis[2]'[:,1],lw = 3)
end

begin
  plot(ref_sens.t[2:end],x'[2:end,1]+dp[1]'[2:end,1]*Δp[1])
  plot!(sol.t[1:end-1],sol_u1[1:end-1]+sensis[2]'[:,1]*Δp[1])
  xlims!((5.0,7.6))
  ylims!((7.49,7.6))
end
