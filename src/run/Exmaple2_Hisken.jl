using DifferentialEquations
using Plots
using SciMLSensitivity
using ForwardDiff
include("C:/Users/liemann/github/PowerDynamicsDevelopment/src/sensitivity_analyses/Local_Sensitivity.jl")

function Example2()
    function f(dx,x,p,t)
        λ = p
        dx[1] = x[4]
        dx[2] = x[6]
        dx[3] = x[3] - x[1]
        dx[4] = x[4] - (x[2] - x[5])
        dx[5] = 0.0
        dx[6] = 0.0
    end
    M = [1. 0  0  0  0  0
         0  1. 0  0  0  0
         0  0  0  0  0  0
         0  0  0  0  0  0
         0  0  0  0  1. 0
         0  0  0  0  0  1.];


    evr = Array{Float64}(undef,0,4)
    s1(u, t, integrator) = u[3]

    function h1(integrator) 
    integrator.u[5] = integrator.u[5] + (1-integrator.p[1]) * integrator.u[4]
    integrator.u[6] = -integrator.u[6]
    evr = vcat(evr, [integrator.t 1 1 1])
    end
    cb1 = ContinuousCallback(s1,h1)

    u0 = [0.25,0,0.25, 0,0,-1]
    p = [0.8]
    ode_fun = ODEFunction(f,mass_matrix = M)
    prob = ODEProblem(ode_fun, u0, (0.0, 6), p)
    sol = solve(prob, Rodas4(), callback = CallbackSet(cb1),dtmax=1e-4)
    return sol, evr
end
sol,evr = Example2();
plot(sol[1,:])
plot!(sol[2,:])
plot!(sol[4,:])
plot(sol,idxs=(1,4))


#######

mtk = modelingtoolkitize(sol.prob)
sym_states = states(mtk)
symp = parameters(mtk)

hs = Vector{Vector{Num}}(undef,0)
ident = Num.(vcat(sym_states[[1,2,5,6]],symp))
ident[3] = sym_states[5] + (1-symp[1])*sym_states[4]
ident[4] = -sym_states[6]
push!(hs,ident)

@parameters t
s = Vector{Num}(undef,0)
push!(s,sym_states[3])

hybrid_sen = CalcHybridTrajectorySensitivity([mtk],sol,evr,s,hs);
plot(sol.t,hybrid_sen[5][1,:])
plot!(sol.t,hybrid_sen[5][7,:])