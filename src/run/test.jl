using Plots
using DifferentialEquations

function test1(du,u,p,t)
    du[1] = p[1]*(u[2]-u[1])
    du[2] = u[1]*(p[2]-u[3]) - u[2]
    x = ifelse(p[3]>0.0,u[1]*u[2] - (8/3)*u[3],0.0)
    du[3] =x
end


u0 = [1.0;0.0;0.0]
tspan = (0.0,30.0)
p = [10.0, 28.0,8/3]
prob = ODEProblem(test1,u0,tspan,p)
mtk  = modelingtoolkitize(prob)

sol = solve(prob,dtmax=1e-4, callback = CallbackSet(cb1,cb2,cb3),tstops=[5.0,10.0,15.0])
plot(sol)
