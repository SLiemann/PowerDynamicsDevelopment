using PowerDynamics
using OrderedCollections: OrderedDict
using DiffEqSensitivity

Ubase = 380e3
Sbase = 100e6
Zbase = (Ubase^2) / (Sbase)

yfault() = 0.1*150.0
tfault_on() = 0.001
tfault_off() = 0.30
dt_max() = 1e-2

function GFC_Test_Grid(;p_new = 0.0,q_new = 0.0,y_new = 0.0)
    buses = OrderedDict(
        "bus1" => SlackAlgebraic(U = 1.00),
        "bus2" => VoltageDependentLoad(
            P = p_new,
            Q = q_new,
            U = 1.0,
            A = 0.0,
            B = 0.0,
            Y_n = y_new,
        ),
        "bus3" => GridFormingConverterCSAAntiWindup(
            Sbase = Sbase,
            Srated = 6*Sbase,
            p0set = 3.0, # based on Sbase!
            q0set = 0.01,
            u0set = 1.00,
            Kp_droop = 0.01,
            Kq_droop = 0.001,
            ωf_P = 10.0 * 2 * pi,
            ωf_Q = 5.0 * 2 * pi,
            xlf = 0.0177*(380/320)^2, #0.01257,     #0.15*(320/380)^2,  *(380/320)^2
            rf =  0.00059095*(380/320)^2, #0.00042,     #0.005*(320/380)^2,
            xcf = 1.7908*(380/320)^2,#1.26990*2,1.7908   # 15.51,# 15.51*(320/380)^2, #1.0/(2.0*pi*50.0*1.231e-6)/Zbase, #
            Kp_u = 0.52, #1.0
            Ki_u = 1.161022,
            Kp_i = 0.738891, # 0.73
            Ki_i = 1.19,
            imax = 1.0,
            Kvi = 0.5, #0.8272172037144201, # 0.677
            σXR = 3.0,
            K_vq = 0.01,
            imax_csa = 1.2,
            p_ind = collect(1:16),
        ),
    )

    branches = OrderedDict(
        "branch1" => PiModelLine(
            from = "bus1",
            to = "bus2",
            y = 1.0 / (1im * 0.1),
            y_shunt_km = 0.0,
            y_shunt_mk = 0.0,
        ),
        "branch2"=> StaticPowerTransformer(from="bus2",to="bus3",Sbase=Sbase,Srated=1200e6,uk=0.15,XR_ratio=Inf,
                                           i0=0.0,Pv0=0.0,tap_side = "LV",tap_pos = 0,tap_inc = 1.0))
    pg = PowerGrid(buses, branches)
end

function GFC_params()
    pg = GFC_Test_Grid()
    GFC = pg.nodes["bus3"]
    return [
        GFC.Kp_droop,
        GFC.Kq_droop,
        GFC.ωf_P,
        GFC.ωf_Q,
        GFC.xlf,
        GFC.rf,
        GFC.xcf,
        GFC.Kp_u,
        GFC.Ki_u,
        GFC.Kp_i,
        GFC.Ki_i,
        GFC.imax,
        GFC.Kvi,
        GFC.σXR,
        GFC.K_vq,
        GFC.imax_csa,
    ]
end

function simGFC(prob)
    #pg_new = GFC_Test_Grid(p_new = -1.3)
    pg_new = GFC_Test_Grid(y_new = yfault())
    tstep = [tfault_on(),tfault_off()]
    event_recorder = Array{Float64,2}(undef,0,4+length(GFC_params()))
    function fault_state(integrator)
        new_f = rhs(pg_new)
        op_prob = ODEProblem(new_f, integrator.sol[end], (0.0, 1e-6),integrator.p, initializealg = BrownFullBasicInit())
        ic_new = solve(op_prob,Rodas5())
        integrator.f = new_f
        integrator.cache.tf.f = integrator.f
        integrator.u = ic_new.u[end]
        event_recorder = vcat(event_recorder,[integrator.t 2 integrator.p' 1 1])
    end

    function postfault_state(integrator)
        sol = integrator.sol
        ic_tmp = deepcopy(integrator.sol.u[indexin(tstep[1],integrator.sol.t)[1]])
        ic_tmp = getPreFaultVoltages(pg_new,ic_tmp,deepcopy(sol[end]))
        #ic_tmp = getPreFaultAlgebraicStates(pg_new,ic_tmp,deepcopy(sol[end]))
        op_prob = ODEProblem(prob.f, ic_tmp, (0.0, 1e-6),integrator.p, initializealg = BrownFullBasicInit())
        ic_new = solve(op_prob,Rodas5())
        integrator.f = prob.f
        integrator.cache.tf.f = integrator.f
        integrator.u = ic_new.u[end]
        event_recorder = vcat(event_recorder,[integrator.t 1 integrator.p' 2 1])
    end

    cb = DiscreteCallback(((u,t,integrator) -> t in tstep[1]), fault_state)
    cb1 = DiscreteCallback(((u,t,integrator) -> t in tstep[2]), postfault_state)

    sol = solve(prob, Rodas4(), tstops= tstep,callback = CallbackSet(cb,cb1), dtmax = dt_max(),progress=true)
    return PowerGridSolution(sol,pg_new), event_recorder
end

function GetTriggCondsGFCTest(mtk::ODESystem)
    st = states(mtk)
    @variables t
    s = [
        0.0 ~ t - tfault_on(),
        0.0 ~ t - tfault_off(),
    ]
    return s
end

function GetStateResFunGFCTest(mtk::ODESystem)
    eqs, aeqs, D_states, A_states = GetSymbolicEquationsAndStates(mtk)
    return [zeros(length(D_states),1) .~ D_states]
end
