using DifferentialEquations


function ftmp(dx,x,p,t)
    input,T1,T2 = p
    dx[1] = (input-x[1])/1e-5 # möglichst kleine Zahl (am besten vermutlich so groß wie die Schrittweite)
    dx[2] = (input+dx[1]*T1 + -x[2])/T2
end

f1 = ODEFunction(ftmp)
prob = ODEProblem(f1, [0,0], (0.0, 10),[0,2,1])

affect!(integrator) = integrator.p[1] = 1.0
cb = PresetTimeCallback([1.0],affect!)

sol = solve(prob,Tsit5(),callback = cb,tstops=[1.0])

plot(sol)

