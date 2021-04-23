using DifferentialEquations
using PowerDynamics
using OrderedCollections: OrderedDict

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

pg = PowerGrid(buses, branches)
pg_postfault = PowerGrid(buses_postfault, branches_postfault)

function switch_off(integrator)
    resize!(integrator,4)

    integrator.f = rhs(pg_postfault)
    ic_temp = find_operationpoint(pg_postfault)
    integrator.u = ic_temp.vec

    # descent into madness: the cache is a OrdinayDiffEq.Rodas4Cache
    # https://github.com/SciML/OrdinaryDiffEq.jl/blob/181dcf265351ed3c02437c89a8d2af3f6967fa85/src/caches/rosenbrock_caches.jl#L326
    # which holds an SciMLBase.TimeGradientWrapper (tf)
    # https://github.com/SciML/SciMLBase.jl/blob/9e5ba6f10a0d347d48ba4d2bae4a3c610a589bc1/src/function_wrappers.jl#L1
    # which holds the reference to the `ODEFunction` (f)
    # update this reference and it works.
    integrator.cache.tf.f = integrator.f
end

tfault = 1.0
cb = PresetTimeCallback(tfault, switch_off)
ic = find_operationpoint(pg)
prob = ODEProblem(rhs(pg), ic.vec, (0.0,5.0))

pgsol = solve(prob, Rodas4(), callback=cb)
