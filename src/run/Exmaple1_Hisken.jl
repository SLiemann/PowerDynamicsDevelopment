using DifferentialEquations
using Plots
using SciMLSensitivity
using ForwardDiff
include("C:/Users/liemann/github/PowerDynamicsDevelopment/src/sensitivity_analyses/Local_Sensitivity.jl")

function Example1(;p=[2.75])
    function f(dx,x,p,t)
        dx[1] = x[1] + x[3]*x[2] 
        dx[2] = x[4]*x[1] + x[2]
        dx[3] = 0.0
        dx[4] = 0.0
        dx[5] = 0.0
        dx[6] = x[6]
    end
    M = [1. 0  0  0  0  0
         0  1. 0  0  0  0
         0  0  1. 0  0  0
         0  0  0  1. 0  0
         0  0  0  0  1. 0
         0  0  0  0  0  0];
    evr = Array{Float64}(undef,0,4)
    s1(u, t, integrator) = (integrator.p[1]*u[1] - u[2])/u[5] + u[5] -1
    s2(u, t, integrator) = (0.36*u[1] - u[2])/u[5] + u[5] + 1

    function h1(integrator)
        tmp =  deepcopy(integrator.u[3])
        integrator.u[3] = integrator.u[4] 
        integrator.u[4] = tmp
        integrator.u[5] = -integrator.u[5]
        evr = vcat(evr, [integrator.t 1 1 1])
    end
    function h2(integrator) 
        tmp =  deepcopy(integrator.u[3])
        integrator.u[3] = integrator.u[4] 
        integrator.u[4] = tmp
        integrator.u[5] = -integrator.u[5]
        evr = vcat(evr, [integrator.t 1 2 1])
    end
    cb1 = ContinuousCallback(s1,h1)
    cb2 = ContinuousCallback(s2,h2)

    u0 = [0.0,1,-100,10,1,0]
    
    ode_fun = ODEFunction(f,mass_matrix = M)#
    prob = ODEProblem(ode_fun, u0, (0.0, 0.2), p)
    sol = solve(prob, Rodas4(), callback = CallbackSet(cb1,cb2),dtmax=1e-4)
    return sol, evr
end

sol,evr = Example1();
plot(sol.t,sol[1,:])
plot!(sol.t,sol[2,:])
# plot(sol.t,sol[5,:])
# plot(sol.t,sol[6,:])

# plot(sol,idxs=(1,2),xlim=(-4,2),ylim=(-2,2))
# x = collect(-4:0.1:2)
# y = 0.36.*x
# plot!(x,y)
# plot!(x,2.75/0.36*y)

#######
mtk = modelingtoolkitize(sol.prob)
sym_states = states(mtk)
symp = parameters(mtk)

hs = Vector{Vector{Num}}(undef,0)
ident = Num.(vcat(sym_states[1:5],symp))
ident[3] = sym_states[4]
ident[4] = sym_states[3]
ident[5] = -sym_states[5]
push!(hs,ident)

@parameters t
s = Vector{Num}(undef,0)
push!(s,(symp[1]*sym_states[1]-sym_states[2])/sym_states[5]+sym_states[5]-1)
push!(s,(0.36*sym_states[1]-sym_states[2])/sym_states[5]+sym_states[5]+1)

hybrid_sen,Δτ = CalcHybridTrajectorySensitivity([mtk],sol,evr,s,hs);
# plot(sol.t,hybrid_sen[6][1,:])
# plot!(sol.t,hybrid_sen[6][2,:])

############
Δp = 3.0 - 2.75
sol_appr = ApproximatedTrajectory(sol,hybrid_sen[6],Δp)
sol_per,evr = Example1(p=[3.0]);

scatter(sol.t,sol[1,:])
scatter!(sol.t,sol_appr[1,:],xlims=(0.11,0.12),ylims=(-0.55,-0.49))
scatter!(sol_per.t,sol_per[1,:])

scatter(sol.t,hybrid_sen[6][1,:],xlims=(0.11,0.12))

# ind1 = Anfang der Indizes für die kontinuierlichen Sensis
ind1 = vcat(1,setdiff(indexin(evr[:,1],sol.t),[nothing]).+1)
# ind2 = Ende der Indizes für die kontinuierlichen Sensis
ind2 = vcat(setdiff(indexin(evr[:,1],sol.t),[nothing]),length(sol.t))
# zwischen den Anfang und End Indizes liegt das Event mit Zustandsänderung
ind_sol = Int.(hcat(ind1,ind2))