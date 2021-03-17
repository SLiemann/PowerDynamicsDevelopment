using PowerDynamics: SlackAlgebraic, VoltageDependentLoad
using PowerDynamics
using OrderedCollections: OrderedDict
using ModelingToolkit

@variables x #x is here a Num

iszero(0*x) #true
iszero(0.0*x) #false


buses=OrderedDict(
    "bus1"=> SlackAlgebraic(U=1),
    "bus2"=> VoltageDependentLoad(P=-0.3, Q=0.3, U=1.0, A=0., B=0.,Y_n = 0.0))

branches=OrderedDict(
    "branch1"=> PiModelLine(from= "bus1", to = "bus2",y=1.0/(0.05+1im*0.15), y_shunt_km=0., y_shunt_mk=0.))
pg = PowerGrid(buses, branches)
ic = rand(4,1)
ODEProb = ODEProblem(rhs(pg),ic,(0,1.0))
mtsys = modelingtoolkitize(ODEProb)
