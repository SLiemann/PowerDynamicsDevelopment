using PowerDynamics
using OrderedCollections: OrderedDict
using Plots
using DifferentialEquations
using CSV #read PF DataFrames
using DataFrames #for CSV

include("C:/Users/liemann/github/PowerDynamicsDevelopment/src/utility/utility_functions.jl")
include("C:/Users/liemann/github/PowerDynamicsDevelopment/src/operationpoint/InitializeInternals.jl")
include("C:/Users/liemann/github/PowerDynamicsDevelopment/src/operationpoint/simulate_fault.jl")

begin
    include(
        "C:/Users/liemann/github/PowerDynamicsDevelopment/src/operationpoint/PowerFlow.jl",
    )
    Ubase = 380e3
    Sbase = 100e6
    Zbase = (Ubase^2) / (Sbase)
    buses = OrderedDict(
        "bus1" => SlackAlgebraic(U = 1.0),
        "bus2" => VoltageDependentLoad(
            P = -0.3,
            Q = -0.3,
            U = 1.0,
            A = 0.0,
            B = 0.0,
            Y_n = complex(0.0),
        ),
    )

    branches = OrderedDict(
        #"branch1"=> DynamicPowerTransformer(from="bus1",to="bus2",Sbase=Sbase,Srated=100e6,uk=0.1581138,XR_ratio=3, i0=6.35,Pv0=300e3,
        #                                   tap_side = "LV",tap_pos = 0,tap_inc = 1.0,tap_delay = 5.0,tap_max = 10,tap_min = -10,
        #                                   deadband_low = 0.99,deadband_high = 1.01, timer_start = -1,Δtap = 0.0,low_high = 0.0))
        "branch1" => StaticPowerTransformerTapParam(
            from = "bus1",
            to = "bus2",
            Sbase = Sbase,
            Srated = 100e6,
            uk = 0.01,
            XR_ratio = 3,
            i0 = 6.35,
            Pv0 = 300e3,
            tap_side = "LV",
            tap_pos = 0,
            tap_inc = 1.0,
            tap_max = 5,
            tap_min = -10,
        ),
    )
    pg = PowerGrid(buses, branches)
    U, δ1, ic = PowerFlowClassic(pg, iwamoto = true, Qlimit_iter_check = 2)
end

#prob = ODEProblem(rhs(pg),ic,(0.0,100.0),[0.0])
#sol = solve(prob, Rodas4(), dtmax = 1e-3)
#pgsol = PowerGridSolution(sol, pg)
pgsol = calcOLTC(pg,ic)
plot(pgsol,collect(keys(pg.nodes)),:v,legend = (0.8, 0.75))

function calcOLTC(pg, ic)
    timer_start = -1.0
    function TapState(integrator)
        timer_start = integrator.t
        integrator.p[1] += 1
    end

    function voltage_deadband(u, t, integrator)
        0.999 <= sqrt(u[3] * u[3] + u[4] * u[4]) <= 1.001
    end

    function timer_off(integrator)
        timer_start = -1
    end

    function voltage_outside(u, t, integrator)
        sqrt(u[3] * u[3] + u[4] * u[4]) < 0.999
    end

    function timer_on(integrator)
        if timer_start == -1
            timer_start = integrator.t
        end
    end

    function timer_hit(u, t, integrator)
        if timer_start == -1
            return false
        else
            return t - timer_start > 5
        end
    end

    cb1 = DiscreteCallback(voltage_deadband, timer_off)
    cb2 = DiscreteCallback(voltage_outside, timer_on)
    cb3 = DiscreteCallback(timer_hit, TapState)

    prob = ODEProblem(rhs(pg), ic, (0.0, 50), [0.0])
    sol = solve(
        prob,
        Rodas4(),
        callback = CallbackSet(cb1, cb2, cb3),
        dtmax = 1e-3,
    )
    return PowerGridSolution(sol, pg)
end
