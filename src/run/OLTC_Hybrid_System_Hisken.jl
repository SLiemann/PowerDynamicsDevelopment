using PowerDynamics
using DifferentialEquations
using Plots
using ModelingToolkit
#using Distributed
#@everywhere using IfElse: ifelse

function calcOLTCHIsken(pg, ic)
    timer_start = -1.0
    pg_fault = OLTC_Hybrid_Sensi(x_grid=0.65)
    function TapState(integrator)
        timer_start = integrator.t
        integrator.p[2] += 1.0
        p_tmp = deepcopy(integrator.p)
        ic_tmp = deepcopy(integrator.sol[end])
        op_prob = ODEProblem(integrator.f, ic_tmp, (0.0, 1e-6), p_tmp, initializealg = BrownFullBasicInit())
        ic_new = solve(op_prob,Rodas5())
        integrator.u = ic_new.u[end]
    end

    function voltage_deadband(u, t, integrator)
        1.04 <= sqrt(u[7] * u[7] + u[8] * u[8])
        #0.99 <= sqrt(u[7] * u[7] + u[8] * u[8]) <= 1.01
    end

    function timer_off(integrator)
        timer_start = -1
    end

    function voltage_outside(u, t, integrator)
        sqrt(u[7] * u[7] + u[8] * u[8]) < 1.04 && t>=10.0
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
            return t - timer_start > 20.0
        end
    end

    function newGridImpedanze(integrator)
        #f = rhs(pg_fault)
        integrator.p[1] = 0.65
        p_tmp = deepcopy(integrator.p)
        ic_tmp = deepcopy(integrator.sol[end])
        op_prob = ODEProblem(integrator.f, ic_tmp, (0.0, 1e-6), p_tmp, initializealg = BrownFullBasicInit())
        ic_new = solve(op_prob,Rodas5())
        #integrator.f = rhs(pg_fault)
        #integrator.cache.tf.f = integrator.f
        #integrator.cache.uf.f = integrator.f
        integrator.u = ic_new.u[end]
    end

    function check_voltage_low(u,t,integrator)
            sqrt(u[7]*u[7] + u[8]*u[8]) < 0.75
    end
    function check_voltage_high(u,t,integrator)
            sqrt(u[7]*u[7] + u[8]*u[8]) > 1.1
    end

    function stop_integration(integrator)
        println("Terminated at t = $(integrator.t)")
        terminate!(integrator)
        #necessary, otherwise PowerGridSolution throws error
        integrator.sol = DiffEqBase.solution_new_retcode(integrator.sol, :Success)
    end

    cb1 = DiscreteCallback(voltage_deadband, timer_off)
    cb2 = DiscreteCallback(voltage_outside, timer_on)
    cb3 = DiscreteCallback(timer_hit, TapState)
    cb4 = DiscreteCallback(((u,t,integrator) -> t in 10.0), newGridImpedanze)
    cb5 = DiscreteCallback(check_voltage_low, stop_integration)
    cb6 = DiscreteCallback(check_voltage_high, stop_integration)

    prob = ODEProblem(rhs(pg), ic, (0.0, 200.0), [0.25,0.0])
    sol = solve(
        prob,
        Rodas4(),
        callback = CallbackSet(cb1, cb2, cb3,cb4,cb5),
        dtmax = 1e-3,
        tstops=[10.0],
        maxiters = 1e5,
    )
    return PowerGridSolution(sol, pg)
end

include("C:/Users/liemann/github/PowerDynamicsDevelopment/src/grids/OLTC_Hybrid_Sensis.jl")

pg, ic0 = GetInitializedOLTCHisken()

pgsol = calcOLTCHIsken(pg,ic0)
plot(pgsol,["bus4"],:v);

mtk = GetMTKOLTCSystem((0.0, 200.0), [0.25,0.0])
eqs = equations(mtk)
st  = states(mtk)
par = parameters(mtk)

include("C:/Users/liemann/github/PowerDynamicsDevelopment/src/grids/LTVS_Test_System.jl")
mtk = GetMTKLTVSSystem((0.0,10.0),nothing)

IfElse.ifelse(1>0,-1,-2)
