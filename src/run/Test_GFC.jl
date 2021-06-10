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

prob = ODEProblem(rhs(pg1),uend,(0.0,100.0),[1])
sol = solve(prob,Rodas4())
pgsol = PowerGridSolution(sol,pg)
plot(pgsol,["bus3"],:iabs)
plot(pgsol,["bus3"],:v)
plot(pgsol,["bus3"],:umabs)
plot(pgsol,["bus3"],:umangle)

plot(pgsol,["bus3"],:p)
plot(pgsol,["bus3"],:θ)
plot(pgsol,["bus3"],:φ)
plot(pgsol,["bus3"],:ω)
plot(pgsol,["bus3"],:Q)

plot(pgsol,["bus3"],:e_uq)
plot(pgsol,["bus3"],:e_ud)

rhs(pg).syms .=> ic

uend = pgsol.dqsol.u[end]

[rhs(pg).syms uend ic]
ic = uend
