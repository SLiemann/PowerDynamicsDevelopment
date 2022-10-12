using DifferentialEquations


function ftmp(dx,x,p,t)
    input,T1,T2 = p
    dx[1] = (dx[2]*T2 + x[2] -input)/T1
    dx[2] = (dx[1]*T1 + x[1] -x[2])/T2
end

f1 = ODEFunction(ftmp)
prob = ODEProblem(f1, [0,0], (0.0, 10),[0,1,1])

affect!(integrator) = integrator.p[1] = 1.0
cb = PresetTimeCallback([1.0],affect!)

sol = solve(prob,Tsit5(),callback = cb,tstops=[1.0])

plot(sol)

