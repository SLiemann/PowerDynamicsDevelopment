using PowerDynamics
using OrderedCollections: OrderedDict
using Distributed
import DiffEqBase: initialize_dae!
@everywhere using IfElse
include("C:/Users/liemann/github/PowerDynamicsDevelopment/src/include_custom_nodes_lines_utilities.jl")


buses=OrderedDict(
    "bus_slack" => SlackAlgebraicParam(U = 1.0,p_ind=[1]),
    "bus_load" => WeccPeLoad(P=-0.5, Q=-0.1, Vd1=0.8, Vd2=0.4,p_ind=[2]))

branches=OrderedDict(
    "line"=> PiModelLine(from= "bus_slack", to = "bus_load",y=1.0/(0.01+1im*0.1), y_shunt_km=0.0, y_shunt_mk=0.0))
pg=  PowerGrid(buses, branches)

U1,δ1,ic0,cu = PowerFlowClassic(pg,iwamoto = true,max_tol = 1e-7,ind_sl=1)
pg, ic = InitializeInternalDynamics(pg,ic0)

index_U_load = PowerDynamics.variable_index(pg.nodes,"bus_load",:u_r)

function peload_decrease(u, t, integrator)
    pg.nodes["bus_load"].Vd2 <= sqrt(u[index_U_load]*u[index_U_load] + u[index_U_load+1]*u[index_U_load+1]) < pg.nodes["bus_load"].Vd1
end

function peload_disconnect(u, t, integrator)
    sqrt(u[index_U_load]*u[index_U_load] + u[index_U_load+1]*u[index_U_load+1]) < pg.nodes["bus_load"].Vd2 
end

function affect_peload_decrease(integrator)
    U = sqrt(integrator.u[index_U_load]*integrator.u[index_U_load] + integrator.u[index_U_load+1]*integrator.u[index_U_load+1])
    vd1 = pg.nodes["bus_load"].Vd1
    vd2 = pg.nodes["bus_load"].Vd2
    integrator.p[pg.nodes["bus_load"].p_ind[1]] =  (abs(U)-vd2)/(vd1 - vd2)
    initialize_dae!(integrator,BrownFullBasicInit())
    #auto_dt_reset!(integrator)
end

function affect_peload_disconnect(integrator)
    integrator.p[pg.nodes["bus_load"].p_ind[1]] =  0.0
    initialize_dae!(integrator,BrownFullBasicInit())
    #auto_dt_reset!(integrator)
end

function voltage_decrease(u, t, integrator)
    1.0 <= t < 2.0
end

function voltage_increase(u, t, integrator)
    2.0 <= t < 3.0
end

function affect_voltage_decrease(integrator)
    u0 = 2.0-integrator.t
    if u0 <=0.1 
         u0 = 0.1 
    end
    integrator.p[pg.nodes["bus_slack"].p_ind[1]] =  u0
end

function affect_voltage_increase(integrator)
    u0 = integrator.t - 1.9
    if u0 >=1.0 
         u0 = 1.0 
    end
    integrator.p[pg.nodes["bus_slack"].p_ind[1]] =  u0
end

cb1 = DiscreteCallback(peload_decrease, affect_peload_decrease)
cb2 = DiscreteCallback(peload_disconnect, affect_peload_disconnect)
cb3 = DiscreteCallback(voltage_decrease, affect_voltage_decrease)
cb4 = DiscreteCallback(voltage_increase, affect_voltage_increase)

params = [1.0,1.0]
problem= ODEProblem{true}(rhs(pg),ic,(0.0,4.0),params)
sol = solve(problem, Rodas4(), callback = CallbackSet(cb1,cb2,cb3,cb4), dtmax = 1e-3,force_dtmin=false,maxiters=1e4, initializealg = BrownFullBasicInit())
pgsol = PowerGridSolution(sol,pg);
plotallvoltages(pgsol);
myplot(pgsol,"bus_load",:P0)
myplot(pgsol,"bus_load",:Q0)

