using DifferentialEquations
using Plots

function myfun(du,u,p,t)
  du[1] = (p[1]-u[1])/5.0
end

function myfun2(du,u,p,t)
  du[1] = p[4]*(-1.0/p[2]*u[1] + (1.0 - p[3]/p[2])*p[1])
  du[2] = u[1]/p[2] + (p[3]/p[2])*p[1]
end
#dx1 = K * (-1 / T1 * x1 + (1 - T2 /  T1) * ω) # Block Input

#dP = (1 / T1) * x1 + (T2 / T1) * ω

function myfun3(du,u,p,t)
  du[1] = (p[1]*p[4]-u[1])/p[2]
  du[2] = u[1] + du[1]*p[3] -u[2]
end

u0 = [0.0,0.0]
p  = [0.0,20.0,10.0,70.0]
tspan = (0.0,180.0)
M = [1.0 0; 0 0]
f = ODEFunction(myfun3,mass_matrix=M)
prob = ODEProblem(f,u0,tspan,p) #
#prob = ODEProblem(myfun2,u0,tspan,p)

jumptime = [30.0]
condition(u,t,integrator) = t ∈ jumptime
function affect!(integrator)
   integrator.p[1] = 1
 end

cb = DiscreteCallback(condition,affect!)
sol = solve(prob,Rodas4(), callback = cb, tstops=jumptime, dt=1e-3, adaptive = false)
#sol = solve(prob,Rodas4())
plot(sol)
