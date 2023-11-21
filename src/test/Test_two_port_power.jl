using DifferentialEquations
using Plots
using LinearAlgebra

function f(du,u,p,t)
    ω  = p[1]
    U1 = p[2]
    U2 = p[3]
    K = p[4]
    L = p[5]
    R = p[6]

    du[1] = u[1] - U1*sin(ω*t)
    du[2] = (u[1] - u[2]*R - U2*sin(ω*t+u[3])) / L
    du[3] = K
end

mm = Diagonal([false,true,true])
fun = ODEFunction(f,mass_matrix =  mm)

u0 = [0.0,0.0,0.0]
p = [100*pi,230*sqrt(2),230*sqrt(2),0.1,1/(100*pi),1]
prob = ODEProblem(fun,u0,(0,1),p)

sol = solve(prob,Rodas4(),dtmax=1e-6)

plot(sol,vars=(2))

i = sol[2,:]
u1 = sol[1,:]
δ = sol[3,:]
um = u1 .- i*1
plot(sol.t,230*sqrt(2).*sin.(100*pi.*sol.t+δ))
plot(sol.t,um)
xlims!((0.9,1))


using NLsolve


function VoltageCurve(P,Q)
    function f(F,u,p)
        Z  = p[1]
        Q = p[2]
        U1 = p[3]
        X = p[4]


        F[1] = U1*u[1]/X*sin(u[2]) - u[1]^2/Z
        F[2] = u[1]^2/X-U1*u[1]*cos(u[2])/X - Q
    end

    p = [1,0,1,0.1]
    p[1] = P
    p[2] = Q
    sol = nlsolve((F,x) ->f(F,x,p), [1.0,pi/8],ftol = 1e-12)
    return sol.zero
end

Zload = 1:-0.001:0.001
Qload = -0:0.1:1
V2 = VoltageCurve.(Zload,0)

U2 = zeros(length(V2))
winkel = zeros(length(V2))
for (ind,val) in enumerate(V2)
    U2[ind] = val[1]
    winkel[ind] = rem(val[2]*180/pi,180)
end
Pload = U2.^2.0 ./Zload

plot(Pload,U2)
plot!(Pload,winkel./180*pi)
