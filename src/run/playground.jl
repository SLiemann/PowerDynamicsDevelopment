using DifferentialEquations
using PowerDynamics
using OrderedCollections: OrderedDict
using Plots
using DataFrames
using CSV

Ubase = 380e3
Sbase  = 100e6
Zbase = Ubase^2/Sbase
Z_EHV_Line = (9.6 + 1im*64)/Zbase
B_half     = 1im*1498.54*1e-6 / 2.0 *Zbase

include("C:/Users/liemann/github/PowerDynamicsDevelopment/src/utility/utility_functions.jl")
buses=OrderedDict(
    "bus1" => SlackAlgebraic(U=1.0),
    "busv" => VoltageDependentLoad(P=0.0, Q=0.0, U=1.0, A=0., B=0.,Y_n = complex(0.0)),
    "bus2" => VoltageDependentLoad(P=-1.0, Q=-0.5, U=1.0, A=0., B=0.,Y_n = complex(0.0)))

branches=OrderedDict(
    "Line_1-2"=> PiModelLine(from= "bus1", to = "bus2",y=1.0/Z_EHV_Line, y_shunt_km=B_half, y_shunt_mk=B_half),
    "Line_1-v"=> PiModelLine(from= "bus1", to = "busv",y=1.0/(Z_EHV_Line/2.0), y_shunt_km=B_half, y_shunt_mk=0.0),
    "Line_v-2"=> PiModelLine(from= "busv", to = "bus2",y=1.0/(Z_EHV_Line/2.0), y_shunt_km=0.0, y_shunt_mk=B_half))

buses_fault=OrderedDict(
    "bus1" => SlackAlgebraic(U=1.0),
    "busv" => VoltageDependentLoad(P=0.0, Q=0.0, U=1.0, A=0., B=0.,Y_n = complex(1.0/(1im*1/Zbase))),
    "bus2" => VoltageDependentLoad(P=-1.0, Q=-0.5, U=1.0, A=0., B=0.,Y_n = complex(0.0)))

pg       = PowerGrid(buses, branches)
pg_fault = PowerGrid(buses_fault, branches)

function error_state(integrator)
    sol1 = integrator.sol
    #x2 = find_valid_initial_condition(np_powergrid, sol1[end]) # Jump the state to be valid for the new system.
    ode = rhs(pg_fault)
    op_prob = ODEProblem(ode, sol1[end], (0.0, 1e-6), nothing, initializealg = BrownFullBasicInit())
    x2 = solve(op_prob,Rodas5())
    x2 = x2.u[end]
    integrator.f = rhs(pg_fault)
    integrator.u = x2#sol1[end]
end

function regularState(integrator)
    sol2 = integrator.sol
    #x3 = find_valid_initial_condition(powergrid, sol2[end]) # Jump the state to be valid for the new system.
    ode = rhs(pg)
    ic_tmp = deepcopy(integrator.sol.u[indexin(tfault[1],integrator.sol.t)[1]])
    ic_tmp = getPreFaultVoltages(pg,ic_tmp,deepcopy(sol2[end]))
    op_prob = ODEProblem(ode, ic_tmp, (0.0, 1e-6), nothing, initializealg = BrownFullBasicInit())
    x3 = solve(op_prob,Rodas5(),dtmax=1e-6)
    x3 = x3.u[end]

    integrator.f = rhs(pg)
    integrator.u = x3#sol2[end]
end

tfault = [1.0,1.15]
cb1 = DiscreteCallback(((u,t,integrator) -> t in tfault[1]), error_state)
cb2 = DiscreteCallback(((u,t,integrator) -> t in tfault[2]), regularState)
ic = find_operationpoint(pg)
prob = ODEProblem(rhs(pg),ic.vec,(0.0,5.0))
pgsol = solve(prob,Rodas4(),callback =CallbackSet(cb1,cb2),tstops=tfault)

pgsol = PowerGridSolution(pgsol, pg)
plot(pgsol,collect(keys(pg.nodes)), :v,size = (1000, 500),legend = (0.6, 0.75))

begin
    plot(pgsol,collect(keys(pg.nodes)), :v,size = (1000, 500),legend = (0.6, 0.75))
    test = DataFrame(CSV.File("C:\\Users\\liemann\\Desktop\\u_pf.csv"; header=false, delim=';', type=Float64))
    plot!(test.Column1,test.Column2,label = "PF-bus1")
    plot!(test.Column1,test.Column3,label = "PF-bus2")
    plot!(test.Column1,test.Column4,label = "PF-busv")
end


#plot!(pgsol)

#=
ic_new = deepcopy(pgsol.u[end])
deleteat!(ic_new,3:4)
prob_fault = ODEProblem{true}(rhs(pg_postfault),ic_new,(pgsol.t[end],5.0))
pgsol2 = solve(prob_fault,Rodas4())

u2 = deepcopy(pgsol2.u)
for (ind,val) in enumerate(u2)
    u2[ind] = [val[1:2];0.0;0.0;val[3:4]]
end

t = [pgsol.t;pgsol2.t]
u = [pgsol.u;u2]

u_new = [0 0 0 0 0 0]
for (ind,val) in enumerate(u)
    u_new = append!(u_new,val')
end
plot(t,u')
plot([]


integ = init(prob,Rodas4());
deleteat_non_user_cache!(integrator,3:4)
deleteat!(integ,3:4)

integ.f = rhs(pg_postfault)
ic_temp = find_operationpoint(pg_postfault)
integ.u = ic_temp.vec
step!(integ)
=#
