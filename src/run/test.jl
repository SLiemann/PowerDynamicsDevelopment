using Plots
using DifferentialEquations

function test1(du,u,p,t)
    du[1] = p[1]*(u[2]-u[1])
    du[2] = u[1]*(p[2]-u[3]) - u[2]
    du[3] = u[1]*u[2] -1im*p[3]*u[3]
end


u0 = [1.0;0.0;0.0]
tspan = (0.0,30.0)
p = [10.0, 28.0,8/3*1im]
prob = ODEProblem(test1,u0,tspan,p)
mtk  = modelingtoolkitize(prob)

sol = solve(prob,dtmax=1e-4, callback = CallbackSet(cb1,cb2,cb3),tstops=[5.0,10.0,15.0])
plot(sol)

@parameters z

function test(a::Union{Complex{Float64},Complex{Num}},b,c::Union{Complex{Float64},Complex{Num}})
    return a+b+c
end
r = test(z*im,2.0,2.0*im);

typeof(1.0+z)
