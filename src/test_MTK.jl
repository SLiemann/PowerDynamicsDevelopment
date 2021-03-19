using PowerDynamics: SlackAlgebraic, VoltageDependentLoad, SlackAlgebraicParam
using PowerDynamics
import PowerDynamics: PiModel
using OrderedCollections: OrderedDict
using DifferentialEquations
#using ModelingToolkit
using Plots
#using DiffEqSensitivity

buses = OrderedDict(
    "bus1" => SlackAlgebraicParam(U = 0.98),
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
ODEProb = ODEProblem(rhs(pg), ic, (0, 1.0),1.0)
new_f = ODEFunction(ODEProb.f.f, syms = ODEProb.f.syms, mass_matrix = Int.(ODEProb.f.mass_matrix))
ODEProb = ODEProblem(new_f,ic,(0,1.5), 1.0,initializealg = BrownFullBasicInit()) # 1 is a dummy parameter
mtsys = modelingtoolkitize(ODEProb)


#mtsys = modelingtoolkitize(ODEProb)
copy_ode = deepcopy(ODEProb)
function condition(u,t,integrator)
    t == 1.0
end

function affect!(integrator)
    tmp_prob = ODEProblem(integrator.f,integrator.sol[end],(0.0,1e-6),0.8),initializealg = BrownFullBasicInit())
    sol = solve(tmp_prob, Rodas4())
    integrator.p = 0.8
    integrator.u = sol[end]
end

cb1 = DiscreteCallback(condition, affect!)

sol = solve(ODEProb, Rodas4(), callback = cb1, dt = 1e-4, adaptive = false, tstops = [1.0, 1.5])

plot(sol)
ylims!((0.5,1.05))
xlims!((0.999,1.001))
