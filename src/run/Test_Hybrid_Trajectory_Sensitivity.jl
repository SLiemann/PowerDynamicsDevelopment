using Plots
using DifferentialEquations
using ModelingToolkit
using LinearAlgebra


function example1(du,u,p,t)
    #display(t)
    du[1] = u[1] + u[4]*u[2]
    du[2] = u[5]*u[1] + u[2]
    if u[3] < 0.0
        du[3] = p[1]*u[1] -u[2] -u[6]*u[3]
    else
        du[3] = u[2] -0.36*u[1] -u[6]*u[3]
    end
    du[4] = 0.0
    du[5] = 0.0
    du[6] = 0.0
end

function example11(du,u,p,t)
    du[1] = u[1] + p[2]*u[2]
    du[2] = p[3]*u[1] + u[2]
end

condition_up(u,t,integrator) = integrator.p[1]*u[1]-u[2]
condition_down(u,t,integrator) = u[2]-0.36*u[1]

function affect(integrator)
    tmp = deepcopy(integrator.p[2])
    integrator.p[2] = integrator.p[3]
    integrator.p[3] = tmp
end

x0 = [0.0,1.0]
p  = [2.75,-100.0,10.0]
cb1 = ContinuousCallback(condition_down,affect,nothing)
cb2 = ContinuousCallback(condition_up,affect,nothing)

fun = ODEFunction(example11,mass_matrix=Diagonal([1,1]))
prob =  ODEProblem(fun,x0,(0.0,0.2),p)
sol = solve(prob,Rodas4(),callback=CallbackSet(cb1,cb2),dtmax=1e-3) #

plot(sol,vars=[(0,1),(0,2)])
