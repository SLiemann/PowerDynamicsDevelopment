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
end
pg = GFC_Test_Grid()
U,δ,ic0 = PowerFlowClassic(pg)
Ykk = NodalAdmittanceMatrice(pg)
Uc = U.*exp.(1im*δ/180*pi)
I_c = Ykk*Uc
S = conj(Ykk*Uc).*Uc

pg1 ,ic = InitializeInternalDynamics(pg,I_c,ic0)
params = GFC_params()
prob = ODEProblem(rhs(pg1),ic,(0.0,100.0),params)
@time pgsol = simGFC(prob)

plot(pgsol,["bus3"],:iabs,legend =false)
ylims!((0.94,1.01))
plot(pgsol,["bus3"],:v)
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

function simGFC(prob)
    #pg_new = GFC_Test_Grid(p_new = -1.3)
    pg_new = GFC_Test_Grid(y_new = 150.0)
    params = GFC_params()
    tstep = [1.0,1.15]
    function fault_state(integrator)
        new_f = rhs(pg_new)
        op_prob = ODEProblem(new_f, integrator.sol[end], (0.0, 1e-6),params, initializealg = BrownFullBasicInit())
        ic_new = solve(op_prob,Rodas5())
        integrator.f = new_f
        integrator.cache.tf.f = integrator.f
        integrator.u = ic_new.u[end]
    end

    function postfault_state(integrator)
        sol = integrator.sol
        ic_tmp = deepcopy(integrator.sol.u[indexin(tstep[1],integrator.sol.t)[1]])
        ic_tmp = getPreFaultVoltages(pg_new,ic_tmp,deepcopy(sol[end]))
        op_prob = ODEProblem(prob.f, ic_tmp, (0.0, 1e-6),params, initializealg = BrownFullBasicInit())
        ic_new = solve(op_prob,Rodas5())
        integrator.f = prob.f
        integrator.cache.tf.f = integrator.f
        integrator.u = ic_new.u[end]
    end

    cb = DiscreteCallback(((u,t,integrator) -> t in tstep[1]), fault_state)
    cb1 = DiscreteCallback(((u,t,integrator) -> t in tstep[2]), postfault_state)

    sol = solve(prob, Rodas4(), tstops= tstep, dtmax = 1e-2,progress=true)#,callback = CallbackSet(cb,cb1)
    return PowerGridSolution(sol,pg_new)
end

#return sol

Imax = 1.0
I_tresh = 0.95
ΔXR = 10.0
V0 = 1.0
a = (Imax - I_tresh)^2 * (1 + ΔXR^2)
b = 2.0*(Imax - I_tresh) * (ΔXR)
c = - V0^2/Imax^2

kp = (-b + sqrt(b^2 - 4.0*a*c))/(2*a)
