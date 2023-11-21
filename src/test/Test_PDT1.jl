using DifferentialEquations
using Plots
using LinearAlgebra



function ftmp(dx,x,p,t)
    input,T1,T2 = p
    dx[1] = (input-x[1])/1e-3 # möglichst kleine Zahl (am besten vermutlich so groß wie die Schrittweite)
    dx[2] = (input+dx[1]*T1 + -x[2])/T2
end

function ftmp2(dx,x,p,t)
    input,T1,T2 = p
    dx[1] = (input-x[1])/T1 
    dx[2] = x[2] - (x[1] +dx[1]*T2)
end

f1 = ODEFunction(ftmp)
f2 = ODEFunction(ftmp2,mass_matrix = Diagonal([true,false]))
prob = ODEProblem(f1, [0,0], (0.0,10),[0,2,1])
prob2 = ODEProblem(f2, [0,0], (0.0,10),[0,2,1])

affect!(integrator) = integrator.p[1] = 1.0
cb = PresetTimeCallback([1.0],affect!)

sol = solve(prob,Tsit5(),dtmax=1e-3,callback = cb,tstops=[1.0])
sol2 = solve(prob2,Rodas4(),dtmax=1e-3,callback = cb,tstops=[1.0])

plot(sol)
plot!(sol2)

