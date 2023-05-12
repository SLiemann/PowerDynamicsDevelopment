using DifferentialEquations
using Plots
using SciMLSensitivity

function f(dx,x,p,t)
    Ud,w,L,C,R = p
    dx[1] = (Ud*cos(w*t) - x[2])/L * x[3]
    dx[2] = (x[1] - x[2]/R)/C
    dx[3] = 0
end

s1(u, t, integrator) = u[1]
h1(integrator) = integrator.u[3] = 0.0
cb1 = ContinuousCallback(s1, h1)

function s2(u, t, integrator)
    Ud,w,L,C,R = integrator.p
    (Ud*cos(w*t) - integrator.u[2]) > 0.0
end
h2(integrator) = integrator.u[3] = 1.0
cb2 = DiscreteCallback(s2, h2)

u0 = [0.0,0.0,1.0]
p = [sqrt(2),100*pi,5e-3,1e-3,50]
prob = ODEProblem(f, u0, (0.0, 0.25), p)
sol = solve(prob, Tsit5(), callback = CallbackSet(cb1,cb2),dtmax=1e-3)
plot(sol)

#vars=(1)
