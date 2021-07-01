using PowerDynamics
using DifferentialEquations
using Plots
#using ModelingToolkit
#using CSV
#using DataFrames

begin
    include("C:/Users/liemann/github/PowerDynamicsDevelopment/src/include_costum_nodes_lines_utilities.jl")
    include("C:/Users/liemann/github/PowerDynamicsDevelopment/src/grids/GFC_Test_Grid.jl")
    include("C:/Users/liemann/github/PowerDynamicsDevelopment/src/sensitivity_analyses/Local_Sensitivity.jl")

    pg = GFC_Test_Grid()
    U,δ,ic0 = PowerFlowClassic(pg, iwamoto = true, max_tol = 1e-7)
    Ykk = NodalAdmittanceMatrice(pg)
    Uc = U.*exp.(1im*δ/180*pi)
    I_c = Ykk*Uc
    S = conj(Ykk*Uc).*Uc

    pg1 ,ic = InitializeInternalDynamics(pg,I_c,ic0)
    params = GFC_params()
    prob = ODEProblem(rhs(pg1),ic,(0.0,10.0),params)
end
@time pgsol,evr = simGFC(prob)

plot(pgsol,["bus3"],:i_abs, label = "Kvi = " * string(pg.nodes["bus3"].Kvi) * ", K_vq = " *string(pg.nodes["bus3"].K_vq))
ylims!((2.5,2.99))
xlims!((0.99,2.5))
plot(pgsol,collect(keys(pg.nodes))[2:end],:v)
plot(pgsol,["bus2"],:iabs)
plot(pgsol,["bus3"],:p)
plot(pgsol,["bus3"],:q)
plot(pgsol,["bus3"],:θ)
plot(pgsol,["bus3"],:φ)
plot(pgsol,["bus3"],:ω)
plot(pgsol,["bus3"],:Qm)
xlims!((4.9,5.50))

plot(pgsol,["bus3"],:e_uq)
plot(pgsol,["bus3"],:e_ud)

mtk = GetMTKSystem(pg,(0.0,10.0),params)
tmp = CalcEigenValues(pg1,params,output = true)
#Checking initialization
[rhs(pg).syms pgsol.dqsol.u[end] ic]


#return sol

Imax = 1.0
I_tresh = 0.95
ΔXR = 10.0
V0 = 1.0
a = (Imax - I_tresh)^2 * (1 + ΔXR^2)
b = 2.0*(Imax - I_tresh) * (ΔXR)
c = - V0^2/Imax^2

kp = (-b + sqrt(b^2 - 4.0*a*c))/(2*a)
