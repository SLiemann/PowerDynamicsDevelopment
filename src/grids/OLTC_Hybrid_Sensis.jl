using PowerDynamics
using OrderedCollections: OrderedDict
using ModelingToolkit

function OLTC_Hybrid_Sensi(;x_grid = 0.25)
    Ubase = 380e3
    Sbase = 100e6
    Zbase = (Ubase^2) / (Sbase)
    buses = OrderedDict(
        "bus1" => SlackAlgebraic(U = 1.05),
        "bus2" => VoltageDependentLoad(P=0.0,Q=0.0,U=1.0,A=0.0,B=0.0,Y_n = 0.0),
        "bus3" => VoltageDependentLoad(P=0.0,Q=0.0,U=1.0,A=0.0,B=0.0,Y_n = complex(0.0)),
        "bus4" => SimpleRecoveryLoadParam(
            P0 = -0.4,
            Q0 = 0.0,
            Pt = -0.4,
            Qt = 0.0,
            p_ind = [3,4],
        ),
    )

    branches = OrderedDict(
        "branch1" => PiModelLineParam(from="bus1",to="bus2",y=1.0/(1im*x_grid),y_shunt_km=0.0,y_shunt_mk=0.0,p_ind=1),
        "branch2" => StaticPowerTransformerTapParam(
            from = "bus2",
            to = "bus3",
            Sbase = Sbase,
            Srated = 100e6,
            uk = 0.01,
            XR_ratio = Inf,
            i0 = 0.0,
            Pv0 = 0.0,
            tap_side = "LV",
            tap_pos = 3,
            tap_inc = 1.25,
            tap_max = 8,
            tap_min = -8,
            p_ind=2
        ),
        "branch3" => PiModelLine(from="bus3",to="bus4",y=1.0/(1im*0.80/1.0125),y_shunt_km=0.0,y_shunt_mk=0.0),
    )
    pg = PowerGrid(buses, branches)
end
tspanOLTCHisken() = (0.0,200.0)
GetParametersOLTCHisken(Tp) =  [0.25,0.0,Tp,5.0,20]

function GetInitializedOLTCHisken()
    pg = OLTC_Hybrid_Sensi()
    U,δ,ic0 = PowerFlowClassic(pg,iwamoto = false)
    return InitializeInternalDynamics(pg,ic0)
end

function GetMTKOLTCSystem(;Tp = 5.0)
    pg, ic = GetInitializedOLTCHisken()
    tspan = tspanOLTCHisken()
    params = GetParametersOLTCHisken(Tp)
    prob   = ODEProblem(rhs(pg),ic,tspan,params)
    new_f = ODEFunction(prob.f.f, syms = prob.f.syms, mass_matrix = Int.(prob.f.mass_matrix))
    ODEProb = ODEProblem(new_f,ic,tspan,params)
    return modelingtoolkitize(ODEProb)
end

function SimulateOLTCHIsken(;Tp = 5.0)
    pg, ic = GetInitializedOLTCHisken()
    timer_start = -1.0
    event_recorder = Array{Float64,2}(undef,0,3+length(GetParametersOLTCHisken(Tp)))
    function TapState(integrator)
        timer_start = integrator.t
        tap_max = pg.lines["branch2"].tap_max
        tap_min = pg.lines["branch2"].tap_min
        tap_pos = pg.lines["branch2"].tap_pos
        if integrator.p[2] + tap_pos < tap_max && integrator.p[2] + tap_pos > tap_min
            integrator.p[2] += 1.0
            p_tmp = deepcopy(integrator.p)
            ic_tmp = deepcopy(integrator.sol[end])
            op_prob = ODEProblem(integrator.f, ic_tmp, (0.0, 1e-6), p_tmp, initializealg = BrownFullBasicInit())
            ic_new = solve(op_prob,Rodas5())
            integrator.u = ic_new.u[end]
            event_recorder = vcat(event_recorder,[integrator.t integrator.p' 2 1])
        end
    end

    function voltage_deadband(u, t, integrator)
        1.04 <= sqrt(u[7] * u[7] + u[8] * u[8])
    end

    function timer_off(integrator)
        if timer_start != -1
            timer_start = -1
            event_recorder = vcat(event_recorder,[integrator.t integrator.p' 4 1])
        end
    end

    function voltage_outside(u, t, integrator)
        sqrt(u[7] * u[7] + u[8] * u[8]) < 1.04 && t>10.0
    end

    function timer_on(integrator)
        if timer_start == -1
            timer_start = integrator.t
            event_recorder = vcat(event_recorder,[integrator.t integrator.p' 4 1])
        end
    end

    function timer_hit(u, t, integrator)
        if timer_start == -1
            return false
        else
            return t - timer_start > integrator.p[5]
        end
    end

    function newGridImpedanze(integrator)
        integrator.p[1] = 0.65
        p_tmp = deepcopy(integrator.p)
        ic_tmp = deepcopy(integrator.sol[end])
        op_prob = ODEProblem(integrator.f, ic_tmp, (0.0, 1e-6), p_tmp, initializealg = BrownFullBasicInit())
        ic_new = solve(op_prob,Rodas5())
        integrator.u = ic_new.u[end]
        event_recorder = vcat(event_recorder,[integrator.t integrator.p' 1 1])
    end

    function check_voltage_low(u,t,integrator)
            sqrt(u[7]*u[7] + u[8]*u[8]) < 0.45
    end
    function check_voltage_high(u,t,integrator)
            sqrt(u[7]*u[7] + u[8]*u[8]) > 1.4
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

    params = GetParametersOLTCHisken(Tp)
    #prob = ODEProblem(rhs(pg), ic, tspanOLTCHisken(),params)
    prob = ODEForwardSensitivityProblem(rhs(pg),ic,tspanOLTCHisken(),params,ForwardDiffSensitivity())
    sol = solve(
        prob,
        Rodas4(),
        callback = CallbackSet(cb1, cb2, cb3,cb4,cb5),
        dt = 1e-2,
        adaptive = false,
        tstops=[10.0],
        maxiters = 1e5,
        progress = true,
    )
    #return PowerGridSolution(sol, pg), event_recorder
    return sol, event_recorder
end

function GetTriggCondsOLTCHisken(mtk::ODESystem)
    st = states(mtk)
    @variables t
    s = [
        0.0 ~ t - 10.0,
        0.0 ~ t - 20.0,
        0.0 ~ hypot(st[7], st[8]) - 1.04,
        0.0 ~ hypot(st[7], st[8]) - 1.04,
    ]
    return s
end

function GetStateResFunOLTCHisken(mtk::ODESystem)
    eqs, aeqs, D_states, A_states = GetSymbolicEquationsAndStates(mtk)
    return [zeros(length(D_states),1) .~ D_states]
end
