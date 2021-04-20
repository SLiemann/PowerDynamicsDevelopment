using DifferentialEquations
using ModelingToolkit
using Plots
using IfElse
using CSV #read PF DataFrames
using DataFrames #for CSV

function LowAntiWindupIntegrator(du,u,p,t)
  du[1] = IfElse.ifelse(u[1]<=p[2],IfElse.ifelse(p[1]<0.0,0.0,p[1]),p[1])
end

function HighLowAntiWindupIntegrator(du,u,p,t)
  lowlimit  = IfElse.ifelse(u[1]<=p[2],IfElse.ifelse(p[1]<0.0,true,false),false)
  highlimit = IfElse.ifelse(u[1]>=p[3],IfElse.ifelse(p[1]>0.0,true,false),false)
  du[1] = IfElse.ifelse(lowlimit==true,0.0,IfElse.ifelse(highlimit==true,0.0,p[1]))
end

function myfun_switch(du,u,p,t)
  du[1] = p[1]
  du[2] = IfElse.ifelse(u[1]>0,u[1],IfElse.ifelse(u[1] >= -0.1,0,-1))
end

function AVR(du,u,p,t)
  #p[1] = ΔV0
  du[1] = (p[1]*70.0-u[1])/20.0
  du[2] = u[2] - (u[1]+du[1]*10.0) #algebraic

  e = 10.0*(u[2] - u[3])
  lowlimit  = IfElse.ifelse(u[3]<=0.0,IfElse.ifelse(e<0.0,true,false),false)
  highlimit = IfElse.ifelse(u[3]>=4.0,IfElse.ifelse(e>0.0,true,false),false)
  du[3] = IfElse.ifelse(lowlimit==true,0.0,IfElse.ifelse(highlimit==true,0.0,e))
end

function OEL(du,u,p,t)
  #p[1] = Ifd
  ε = p[1] - 1.8991
  array_out  = IfElse.ifelse(ε > 0.0,ε,IfElse.ifelse(ε >= -0.1,0.0,-1.0))
  du[1] = IfElse.ifelse(u[1]<=-11.0,IfElse.ifelse(array_out<0.0,0.0,array_out),array_out)
  du[2] = u[2] - IfElse.ifelse(u[1]<0.0,-8.0,-ε) # algebraic
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

#u0 = [70*0.02603,70*0.02603,70*0.02603]
u0 = [-11.0,-8.0]
#p  = [0.0,20.0,10.0,70.0]
#p = [-1.0,-3.0,5.0]
p = [1.7]
tspan = (0.0,60.0)
M = [1 0 0; 0 0 0;0 0 1]
M = [1 0; 0 0]
f = ODEFunction(OEL,mass_matrix=M)
prob = ODEProblem(f,u0,tspan,p) #
#prob = ODEProblem(myfun2,u0,tspan,p)
mtprob = modelingtoolkitize(prob)

jumptime = [3.0, 5.0,7,8.0]
condition(u,t,integrator) = t ∈ jumptime
function affect!(integrator)
   integrator.p[1] += 0.15
 end

cb = DiscreteCallback(condition,affect!)
sol = solve(prob,Rodas4(), callback = cb, tstops=jumptime, dtmax=1e-3)
#sol = solve(prob,Rodas4(), dtmax=1e-3)
plot(sol)

begin
    plot(sol, vars = (3))
    test2 = DataFrame(CSV.File("C:\\Users\\liemann\\Desktop\\Matlab.csv"; header=false, delim=';', type=Float64))
    plot!(test2.Column1,(test2.Column2))
    #xlims!((0.9,10.))
end
