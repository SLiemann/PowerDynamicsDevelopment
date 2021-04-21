using DifferentialEquations
using PowerDynamics
using OrderedCollections: OrderedDict
using Plots

Ubase = 380e3
Sbase  = 100e6
Zbase = Ubase^2/Sbase
Z_EHV_Line = (9.6 + 1im*64)/Zbase
B_half     = 1im*1498.54*1e-6 / 2.0 *Zbase

buses=OrderedDict(
    "bus1" => SlackAlgebraic(U=1.0),
    "busv" => VoltageDependentLoad(P=0.0, Q=0.0, U=1.0, A=0., B=0.,Y_n = complex(0.0)),
    "bus2" => VoltageDependentLoad(P=-1.0, Q=-0.5, U=1.0, A=0., B=0.,Y_n = complex(0.0)))

branches=OrderedDict(
    "Line_1-2"=> PiModelLine(from= "bus1", to = "bus2",y=1.0/Z_EHV_Line, y_shunt_km=B_half, y_shunt_mk=B_half),
    "Line_1-v"=> PiModelLine(from= "bus1", to = "busv",y=1.0/(Z_EHV_Line/2.0), y_shunt_km=B_half, y_shunt_mk=0.0),
    "Line_v-2"=> PiModelLine(from= "busv", to = "bus2",y=1.0/(Z_EHV_Line/2.0), y_shunt_km=0.0, y_shunt_mk=B_half))


buses_postfault =OrderedDict(
    "bus1" => SlackAlgebraic(U=1.0),
    "bus2" => VoltageDependentLoad(P=-1.0, Q=-0.5, U=1.0, A=0., B=0.,Y_n = complex(0.0)))

branches_postfault=OrderedDict(
    "Line_1-2"=> PiModelLine(from= "bus1", to = "bus2",y=1.0/Z_EHV_Line, y_shunt_km=B_half, y_shunt_mk=B_half))

pg           = PowerGrid(buses, branches)
pg_postfault = PowerGrid(buses_postfault, branches_postfault)

function switch_off(integrator)
    deleteat!(integrator,3:4)
    #deleteat_non_user_cache!(integrator,3:4) #tried also, but same error
    #resize!(integrator,4) #tried also, but same error
    #resize_non_user_cache!(integrator,4) #tried also, but same error
    integrator.f = rhs(pg_postfault)
    integrator.cache.uf.f = integrator.f
    ic_temp = find_operationpoint(pg_postfault)
    integrator.u = ic_temp.vec
    # reinit!(integrator,ic_temp.vec)  #tried also, but same error
end

tfault =1.0
cb = DiscreteCallback(((u,t,integrator) -> t in tfault), switch_off)
ic = find_operationpoint(pg)
prob = ODEProblem(rhs(pg),ic.vec,(0.0,5.0))

pgsol = solve(prob,Rodas4(),callback =cb,tstops=[tfault])


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
