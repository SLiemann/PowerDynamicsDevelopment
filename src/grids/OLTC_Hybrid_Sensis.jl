using PowerDynamics
using OrderedCollections: OrderedDict
using ModelingToolkit

include("C:/Users/liemann/github/PowerDynamicsDevelopment/src/operationpoint/PowerFlow.jl")
include("C:/Users/liemann/github/PowerDynamicsDevelopment/src/utility/utility_functions.jl")
include("C:/Users/liemann/github/PowerDynamicsDevelopment/src/operationpoint/InitializeInternals.jl")

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
            XR_ratio = 999,
            i0 = 0.0,
            Pv0 = 0.0,
            tap_side = "LV",
            tap_pos = 3,
            tap_inc = 1.25,
            tap_max = 8,
            tap_min = -8,
            p_ind=2
        ),
        "branch3" => PiModelLine(from="bus3",to="bus4",y=1.0/(1im*0.8),y_shunt_km=0.0,y_shunt_mk=0.0),
    )
    pg = PowerGrid(buses, branches)
end
tspanOLTCHisken() = (0.0,200.0)
GetParametersOLTCHisken() =  [0.25,0.0,5.0,5.0]

function GetInitializedOLTCHisken()
    pg = OLTC_Hybrid_Sensi()
    U,δ,ic0 = PowerFlowClassic(pg,iwamoto = true)
    Ykk = NodalAdmittanceMatrice(pg)
    Uc = U.*exp.(1im*δ/180*pi)
    I_c = Ykk*Uc
    S = conj(Ykk*Uc).*Uc
    return InitializeInternalDynamics(pg,I_c,ic0)
end

function GetMTKOLTCSystem()
    pg, ic = GetInitializedOLTCHisken()
    tspan = tspanOLTCHisken()
    params = GetParametersOLTCHisken()
    prob   = ODEProblem(rhs(pg),ic,tspan,params)
    new_f = ODEFunction(prob.f.f, syms = prob.f.syms, mass_matrix = Int.(prob.f.mass_matrix))
    ODEProb = ODEProblem(new_f,ic,tspan,params)
    return modelingtoolkitize(ODEProb)
end

function SimulateOLTCHIsken()
    pg, ic = GetInitializedOLTCHisken()
    timer_start = -1.0
    event_recorder = Array{Float64,2}(undef,0,3+length(GetParametersOLTCHisken()))
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
            event_recorder = vcat(event_recorder,[integrator.t integrator.p' 3 1])
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
        integrator.p[1] = 0.65
        p_tmp = deepcopy(integrator.p)
        ic_tmp = deepcopy(integrator.sol[end])
        op_prob = ODEProblem(integrator.f, ic_tmp, (0.0, 1e-6), p_tmp, initializealg = BrownFullBasicInit())
        ic_new = solve(op_prob,Rodas5())
        integrator.u = ic_new.u[end]
        event_recorder = vcat(event_recorder,[integrator.t integrator.p' 1 1])
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

    params = GetParametersOLTCHisken()
    prob = ODEProblem(rhs(pg), ic, tspanOLTCHisken(),params)
    sol = solve(
        prob,
        Rodas4(),
        callback = CallbackSet(cb1, cb2, cb3,cb4,cb5),
        dtmax = 1e-3,
        tstops=[10.0],
        maxiters = 1e5,
    )
    return PowerGridSolution(sol, pg), event_recorder
end

GetTriggeringConditions() = GetTriggeringConditions(GetMTKOLTCSystem()) #--> ZU SPEZIFISCH
function GetTriggeringConditions(mtk::ODESystem) #--> ZU SPEZIFISCH
    st = states(mtk)
    @variables t
    h = [
        0.0 ~ t - 10.0,
        0.0 ~ t - 20.0,
        0.0 ~ hypot(st[7], st[8]) - 1.04,
        0.0 ~ hypot(st[7], st[8]) - 1.04,
    ]
    return h
end

GetStateResetFunctions() = GetStateResetFunctions(GetMTKOLTCSystem()) #--> ZU SPEZIFISCH
function GetStateResetFunctions(mtk::ODESystem) #--> ZU SPEZIFISCH
    st = states(mtk)
    return [zeros(length(st),1) .~ st]
end

function CalcTriggerAndStateResetJacobians(mtk::ODESystem)
    eqs, aeqs, x, y = GetSymbolicEquationsAndStates(mtk)
    h = GetTriggeringConditions(mtk) #--> zukünftig Übergabeparameter
    s = GetStateResetFunctions(mtk) #--> zukünftig Übergabeparameter
    hx = Array{Array{Num}}(undef,length(h),1)
    hy = similar(hx)
    sx = Array{Array{Num}}(undef,length(s),1)
    sy = similar(sx)
    for i=1:length(h)
        hx[i] = GetJacobian([h[i]],x)
        hy[i] = GetJacobian([h[i]],y)
    end
    for i=1:length(s)
        sx[i] = GetJacobian(s[i],x)
        sy[i] = GetJacobian(s[i],y)
    end
    return hx,hy,sx,sy
end

function CalcSensitivityAfterJump(
    states,
    params,
    xx0_pre,
    yx0_pre,
    x0_pre,
    x0_post,
    p_pre,
    p_post,
    f,
    g,
    J,
    hx,
    hy,
    sx,
    sy,
)
    f, g, x, y = sys
    fx, fy, gx, gy = J

    subs_pre = [states .=> x0_pre; params .=> p_pre]
    subs_post = [states .=> x0_post; params .=> p_post]

    f_pre = Substitute(f, subs_pre)
    f_post = Substitute(f, subs_post)

    fx_pre = Substitute(fx, subs_pre)
    fy_pre = Substitute(fy, subs_pre)
    gx_pre = Substitute(gx, subs_pre)
    gy_pre = Substitute(gy, subs_pre)

    gx_post = Substitute(gx, subs_post)
    gy_post = Substitute(gy, subs_post)

    hx_pre = Substitute(hx, subs_pre)
    hy_pre = Substitute(hy, subs_pre)
    sx_pre = Substitute(sx, subs_pre)
    sy_pre = Substitute(sy, subs_pre)

    gygx = inv(gy_pre) * gx_pre

    hx_star = hx_pre - hy_pre * gygx

    s_star = sx_pre - sy_pre * gygx

    τx0 = s_star * xx0_pre / (s_star * f_pre)

    xx0_post = xx0_pre * hx_star - (f_post - hx_star * f_pre) * τx0
    yx0_post = inv(gy_post) * gx_post * xx0_post

    return xx0_post, yx0_post
end
