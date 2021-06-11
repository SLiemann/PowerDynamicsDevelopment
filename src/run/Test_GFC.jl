using PowerDynamics
using DifferentialEquations
using Plots
#using ModelingToolkit
#using CSV
#using DataFrames

begin
    include("C:/Users/liemann/github/PowerDynamicsDevelopment/src/include_costum_nodes_lines_utilities.jl")
    include("C:/Users/liemann/github/PowerDynamicsDevelopment/src/grids/GFC_Test_Grid.jl")
end
pg = GFC_Test_Grid()
U,δ,ic0 = PowerFlowClassic(pg)
Ykk = NodalAdmittanceMatrice(pg)
Uc = U.*exp.(1im*δ/180*pi)
I_c = Ykk*Uc
S = conj(Ykk*Uc).*Uc

pg1 ,ic = InitializeInternalDynamics(pg,I_c,ic0)
prob = ODEProblem(rhs(pg1),ic,(0.0,30.0),[1])
pgsol = simGFC(prob)

plot(pgsol,["bus3"],:iabs,legend =false)
ylims!((0.45,1.15))
plot(pgsol,["bus3"],:v)
plot(pgsol,["bus3"],:umabs)
plot(pgsol,["bus3"],:umangle)

plot(pgsol,["bus3"],:p)
plot(pgsol,["bus3"],:q)
plot!(pgsol,["bus3"],:θ)
plot!(pgsol,["bus3"],:φ)
plot(pgsol,["bus3"],:ω)
plot!(pgsol,["bus3"],:Qm)
xlims!((4.9,5.50))

plot(pgsol,["bus3"],:e_uq)
plot(pgsol,["bus3"],:e_ud)

rhs(pg).syms .=> ic

uend = pgsol.dqsol.u[end]

[rhs(pg).syms uend ic]
ic = uend

function simGFC(prob)
    pg_new = GFC_Test_Grid(p_new = -1.3)
    tstep = 5.0
    function Loadstep(integrator)
        new_f = rhs(pg_new)
        op_prob = ODEProblem(new_f, integrator.sol[end], (0.0, 1e-6), initializealg = BrownFullBasicInit())
        ic_new = solve(op_prob,Rodas5())
        integrator.f = new_f
        integrator.cache.tf.f = integrator.f
        integrator.u = ic_new.u[end]
    end

    cb = DiscreteCallback(((u,t,integrator) -> t in tstep), Loadstep)
    sol = solve(prob, Rodas4(), callback = cb, tstops=[tstep], dtmax = 1e-3)
    #return sol
    return PowerGridSolution(sol,pg_new)
end
