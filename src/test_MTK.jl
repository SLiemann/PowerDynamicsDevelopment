using PowerDynamics: SlackAlgebraic, VoltageDependentLoad
using PowerDynamics
import PowerDynamics: PiModel
using OrderedCollections: OrderedDict
using DifferentialEquations
#using ModelingToolkit
using Plots
using DiffEqSensitivity

buses = OrderedDict(
    "bus1" => SlackAlgebraic(U = 1),
    "bus2" => VoltageDependentLoad(
        P = -0.3,
        Q = 0.3,
        U = 1.0,
        A = 0.0,
        B = 0.0,
        Y_n = 0.0,
    ),
)

branches = OrderedDict(
    "branch1" => PiModelLine(
        from = "bus1",
        to = "bus2",
        y = 1.0 / (0.05 + 1im * 0.15),
        y_shunt_km = 0.0,
        y_shunt_mk = 0.0,
    ),
)
pg = PowerGrid(buses, branches)

include("operationpoint/PowerFlow.jl")
U,δ1,ic = PowerFlowClassic(pg,iwamoto = true)
#prob = ODELocalSensitivityProblem(rhs(pg),ic,(0.0,10.0),1)
ODEProb = ODEProblem(rhs(pg), ic, (0, 1.0))
new_f = ODEFunction(ODEProb.f.f, syms = ODEProb.f.syms, mass_matrix = Int.(ODEProb.f.mass_matrix))
ODEProb = ODEProblem(new_f,ic,(0,1.0), 1) # 1 is a dummy parameter
mtsys = modelingtoolkitize(ODEProb)


#mtsys = modelingtoolkitize(ODEProb)
