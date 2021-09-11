using PowerDynamics
using OrderedCollections: OrderedDict
using Plots

include("C:/Users/liemann/github/PowerDynamicsDevelopment/src/include_costum_nodes_lines_utilities.jl")

function GetTestGrid(;y_new = 0.0)
    buses = OrderedDict(
        "bus1" => SlackAlgebraic(U = 1.0),
        "bus2" => VoltageDependentLoad(P=-0.50,Q=0.1,U=1.0,A=0.0,B=0.0,Y_n = 0.0),
        "bus3" => FourthOrderEq(
                    H = 3,
                    P = 1.0,
                    D = 0.15*0,
                    Ω = 2 * pi * 50,
                    E_f = 1.0,
                    T_d_dash = 0.1,
                    T_q_dash = 0.1,
                    X_q_dash = 0.2,
                    X_d_dash = 0.2,
                    X_d = 1.3,
                    X_q = 1.3,
                ),
        "bus4" => VoltageDependentLoad(P=-0.25,Q=0.0,U=1.0,A=0.0,B=0.0,Y_n = y_new),
    )
    branches = OrderedDict(
        "branch12" => PiModelLine(from="bus1",to="bus2",y=1.0/(1im*0.1),y_shunt_km=0.0,y_shunt_mk=0.0),
        "branch13" => PiModelLine(from="bus1",to="bus3",y=1.0/(1im*0.1),y_shunt_km=0.0,y_shunt_mk=0.0),
        "branch23" => PiModelLine(from="bus2",to="bus3",y=1.0/(1im*0.1),y_shunt_km=0.0,y_shunt_mk=0.0),
        "branch24" => PiModelLine(from="bus2",to="bus4",y=1.0/(1im*0.1),y_shunt_km=0.0,y_shunt_mk=0.0),
    )
    return pg = PowerGrid(buses, branches)
end

pg = GetTestGrid()
U,δ,ic0 = PowerFlowClassic(pg)
pg, ic = InitializeInternalDynamics(pg,ic0)
prob1 = ODEProblem(rhs(pg),ic,(0.0,5.0))
pgsol = simTestGrid(prob1)

plot(pgsol,"bus1",:v, legend = (0.9,0.5), label = "U1")
plot!(pgsol,"bus2",:v, legend = (0.9,0.5), label = "U2")
plot!(pgsol,"bus3",:v, legend = (0.9,0.5), label = "U3")
plot!(pgsol,"bus4",:v, legend = (0.9,0.5), label = "U4")
plot(pgsol,"bus4",:θ, legend = (0.9,0.5), label = "θ")
plot(pgsol,"bus4",:ω, legend = (0.9,0.5), label = "omega")


function simTestGrid(prob)
    pg_new = GetTestGrid(y_new = 10)
    tstep = [0.5,0.7]
    function fault_state(integrator)
        new_f = rhs(pg_new)
        op_prob = ODEProblem(new_f, integrator.sol[end], (0.0, 1e-6),[1], initializealg = BrownFullBasicInit())
        ic_new = solve(op_prob,Rodas5())
        integrator.f = new_f
        integrator.cache.tf.f = integrator.f
        integrator.u = ic_new.u[end]
    end

    function postfault_state(integrator)
        sol = integrator.sol
        #ic_tmp = deepcopy(integrator.sol.u[indexin(tstep[1],integrator.sol.t)[1]])
        #ic_tmp = getPreFaultVoltages(pg_new,ic_tmp,deepcopy(sol[end]))
        op_prob = ODEProblem(prob.f, sol.u[end], (0.0, 1e-6),[1], initializealg = BrownFullBasicInit())
        ic_new = solve(op_prob,Rodas5())
        integrator.f = prob.f
        integrator.cache.tf.f = integrator.f
        integrator.u = ic_new.u[end]
    end

    cb = DiscreteCallback(((u,t,integrator) -> t in tstep[1]), fault_state)
    cb1 = DiscreteCallback(((u,t,integrator) -> t in tstep[2]), postfault_state)

    sol = solve(prob, Rodas4(), tstops= tstep,callback = CallbackSet(cb,cb1), dtmax = 1e-2,progress=true)
    return PowerGridSolution(sol,pg_new)
end










Uc1 = U.*exp.(1im*δ/180*pi)
