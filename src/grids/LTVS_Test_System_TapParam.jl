using PowerDynamics
using OrderedCollections: OrderedDict
using Distributed
@everywhere using IfElse
using ModelingToolkit

Sbase = 100e6
Ubase = 380e3
Ibase = Sbase/Ubase/sqrt(3)
Zbase = Ubase^2/Sbase

zfault() = 40.0/Zbase
tfault_on() = 1.0
tfault_off() = 1.1
dt_max() = 1e-2

function GFC_LTVS_Test_SystemTapParam(;nTap = 0.0)
    Sbase = 100e6
    Ubase = 380e3
    Ibase = Sbase/Ubase/sqrt(3)
    Zbase = Ubase^2/Sbase

    buses=OrderedDict(
        "bus1" => SlackAlgebraic(U=1.0),
        "bus2" => VoltageDependentLoad(P=0.0, Q=0.0, U=1.0, A=0., B=0.,Y_n = complex(0.0)),
        "bus3" => VoltageDependentLoad(P=-9.6, Q = -2.8, U=1.0, A=1.0, B=0.0,Y_n = complex(0.0)),
        "bus4" => GridFormingConverterCSAAntiWindup(
            Sbase = Sbase,
            Srated = 5.5*Sbase,
            p0set = 5.0, # based on Sbase!
            q0set = 0.001*0,
            u0set = 1.00,
            Kp_droop = 0.020,
            Kq_droop = 0.001,
            ωf_P = 62.8,
            ωf_Q = 62.8,
            xlf = 0.15, #0.01257,     #0.15*(320/380)^2,  *(380/320)^2
            rf =  0.005, #0.00042,     #0.005*(320/380)^2,
            xcf = 15.1515151515,#1.26990*2,1.7908   # 15.51,# 15.51*(320/380)^2, #1.0/(2.0*pi*50.0*1.231e-6)/Zbase, #
            Kp_u = 0.52, #1.0
            Ki_u = 1.161022,
            Kp_i = 0.738891, # 0.73
            Ki_i = 1.19,
            imax = 1.0,
            Kvi = 0.055, #0.8272172037144201, # 0.677
            σXR = 10.0,
            K_vq = 0.1,
            imax_csa = 1.10,
            #iprio = "none",
            p_ind = collect(1:16),
        ),
        "busv" => VoltageDependentLoad(P=0.0, Q=0.0, U=1.0, A=0., B=0.,Y_n = complex(0.0)))

    Z_EHV_Line = (9.6 + 1im*64)/Zbase
    B_half     = 1im*1498.54*1e-6 / 2.0 *Zbase #already an admittance: 1498.54 = 2*pi*50*4.77001*10e-6
    branches=OrderedDict(
        "Line_1-2"=> PiModelLine(from= "bus1", to = "bus2",y=1.0/Z_EHV_Line, y_shunt_km=B_half, y_shunt_mk=B_half),
        "Line_1-v"=> PiModelLine(from= "bus1", to = "busv",y=1.0/(Z_EHV_Line*0.50), y_shunt_km=B_half, y_shunt_mk=0.0),
        "Line_v-2"=> PiModelLine(from= "bus2", to = "busv",y=1.0/(Z_EHV_Line*0.50), y_shunt_km=B_half, y_shunt_mk=0.0),
        "branch3"=> StaticPowerTransformer(from="bus2",to="bus4",Sbase=Sbase,Srated=600e6,uk=0.15,XR_ratio=Inf,
                                           i0=0.0,Pv0=0.0,tap_side = "HV",tap_pos = 0,tap_inc = 1.0),
        "branch4"=> StaticPowerTransformerTapParam(from="bus2",to="bus3",Sbase=Sbase,Srated=1200e6,uk=0.15,XR_ratio=Inf,
                                           i0=0.0,Pv0=0.0,tap_side = "LV",tap_pos = nTap,tap_inc = 1.2,tap_max = 20,tap_min = 0,p_ind=17))
        return PowerGrid(buses, branches)
end

function GFC_LTVS_params_TapParam()
    pg = GFC_LTVS_Test_SystemTapParam()
    GFC = pg.nodes["bus4"]
    OLTC = pg.lines["branch4"]
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
        OLTC.tap_pos
    ]
end

function GetInitializedLTVSSystem(pg::PowerGrid)
    Qmax   = [Inf, Inf, Inf,Inf, Inf,Inf*sqrt(1-0.95^2)]
    Qmin   = -Qmax
    U,δ,ic0 = PowerFlowClassic(pg,iwamoto = true, Qmax = Qmax, Qmin = Qmin, Qlimit_iter_check = 2,max_tol = 1e-6)
    return InitializeInternalDynamics(pg,ic0)
end

function GetMTKLTVSSystemTapParam(;pg_state = "gfc_normal")
    pg0 = GFC_LTVS_Test_SystemTapParam()
    if pg_state == "gfc_fault"
        pg = GetFaultLTVSPG(pg0)
        ic = zeros(systemsize(pg))
    elseif pg_state == "gfc_postfault"
        pg = GetPostFaultLTVSPG(pg0)
        ic = zeros(systemsize(pg))
    else
        pg, ic = GetInitializedLTVSSystem(pg0)
    end
    p = GFC_LTVS_params_TapParam()
    tspan = (0.0,1.0)
    prob   = ODEProblem(rhs(pg),ic,tspan,p)
    new_f = ODEFunction(prob.f.f, syms = prob.f.syms, mass_matrix = Int.(prob.f.mass_matrix))
    ODEProb = ODEProblem(new_f,ic,tspan,p)
    return modelingtoolkitize(ODEProb)
end

function run_LTVS_simulationTapParam(pg::PowerGrid,ic1::Array{Float64,1},tspan::Tuple{Float64,Float64};t_stop_droop = Inf)
    tfault = [tfault_on(), tfault_off()]

    pg_fault = GetFaultLTVSPG(pg)
    pg_postfault = GetPostFaultLTVSPG(pg)

    params = GFC_LTVS_params_TapParam()
    #problem = ODEProblem{true}(rhs(pg),ic1,tspan,params)
    problem = ODEForwardSensitivityProblem(rhs(pg),ic1,tspan,params,ForwardDiffSensitivity())
    timer_start = -1.0
    timer_now   = 0.0
    branch_oltc = "branch4"
    tap = pg.lines[branch_oltc].tap_pos
    OLTC = deepcopy(pg.lines[branch_oltc])
    tap_max = 20.0
    postfault_state = false
    fault_state = false
    #index_U_oltc = getNodeVoltageSymbolPosition(pg,pg.lines[branch_oltc].to)
    index_U_oltc = PowerDynamics.variable_index(pg.nodes,pg.lines[branch_oltc].to,:u_r)
    index_U_load = PowerDynamics.variable_index(pg.nodes,"bus3",:u_r)
    event_recorder = Array{Float64,2}(undef,0,4+length(params))
    function TapState(integrator)
        timer_start = integrator.t
        tap_max = pg.lines["branch4"].tap_max
        tap_min = pg.lines["branch4"].tap_min
        tap_pos = pg.lines["branch4"].tap_pos
        if integrator.p[17] + tap_pos < tap_max && integrator.p[17] + tap_pos > tap_min
            integrator.p[17] += 1.0
            p_tmp = deepcopy(integrator.p)
            ic_tmp = deepcopy(integrator.sol[end])
            op_prob = ODEProblem(integrator.f, ic_tmp, (0.0, 1e-6), p_tmp, initializealg = BrownFullBasicInit())
            ic_new = solve(op_prob,Rodas5())
            integrator.u = ic_new.u[end]

            active_pg = GetActivePG(fault_state,postfault_state)
            event_recorder = vcat(event_recorder,[integrator.t active_pg integrator.p' 3 1])
        end
    end

    function voltage_deadband(u,t,integrator)
         0.99 <= sqrt(u[index_U_oltc]*u[index_U_oltc] + u[index_U_oltc+1]*u[index_U_oltc+1]) <= 1.01
    end

    function timer_off(integrator)
        if timer_start != -1
            timer_start = -1
            active_pg = GetActivePG(fault_state,postfault_state)
            event_recorder = vcat(event_recorder,[integrator.t active_pg integrator.p' 4 1])
        end
    end

    function voltage_outside(u,t,integrator)
         sqrt(u[index_U_oltc]*u[index_U_oltc] + u[index_U_oltc+1]*u[index_U_oltc+1]) < 0.99
    end

    function timer_on(integrator)
        if timer_start == -1
            timer_start = integrator.t
            active_pg = GetActivePG(fault_state,postfault_state)
            event_recorder = vcat(event_recorder,[integrator.t active_pg integrator.p' 4 1])
        end
    end

    function timer_hit(u,t,integrator)
        if timer_start == -1
            return false
        else
            return t-timer_start > 5.0
        end
    end

    function errorState(integrator)
        sol1 = integrator.sol
        ode = rhs(pg_fault)
        op_prob = ODEProblem(ode, sol1[end], (0.0, 1e-6), params, initializealg = BrownFullBasicInit())
        x2 = solve(op_prob,Rodas5())
        x2 = x2.u[end]

        integrator.f = ode
        integrator.cache.tf.f = integrator.f
        integrator.u = x2
        fault_state = true
        active_pg = GetActivePG(fault_state,postfault_state)
        event_recorder = vcat(event_recorder,[integrator.t active_pg integrator.p' 1 1])
    end

    function regularState(integrator)
        sol = integrator.sol
        ode   = rhs(pg_postfault)
        index = PowerDynamics.variable_index(pg.nodes,"busv",:u_r)

        ic_tmp = deepcopy(integrator.sol.u[indexin(tfault[1],integrator.sol.t)[1]])
        #ic_tmp = sol[end]
        ic_tmp = getPreFaultVoltages(pg,ic_tmp,deepcopy(sol[end]))
        #deleteat!(ic_tmp,index:index+1)
        op_prob = ODEProblem(ode, ic_tmp, (0.0, 1e-6), params, initializealg = BrownFullBasicInit())
        x3 = solve(op_prob,Rodas5())
        x3 = x3.u[end]

        #resize!(integrator,length(ic0)-2)
        integrator.f = rhs(pg_postfault)
        integrator.cache.tf.f = integrator.f
        integrator.u = x3#sol2[end]

        index_U_oltc = PowerDynamics.variable_index(pg_postfault.nodes,pg_postfault.lines[branch_oltc].to,:u_r)
        index_U_load = PowerDynamics.variable_index(pg_postfault.nodes,"bus3",:u_r)
        postfault_state = true
        fault_state = false
        active_pg = GetActivePG(fault_state,postfault_state)
        event_recorder = vcat(event_recorder,[integrator.t active_pg integrator.p' 2 1])
    end

    function check_voltage(u,t,integrator)
            sqrt(u[index_U_load]*u[index_U_load] + u[index_U_load+1]*u[index_U_load+1]) < 0.65
    end

    function stop_integration(integrator)
        println("Terminated at $(integrator.t)")
        terminate!(integrator)
        #necessary, otherwise PowerGridSolution throws error
        integrator.sol = DiffEqBase.solution_new_retcode(integrator.sol, :Success)
    end

    function deactivate_droop(integrator)
        integrator.p[3] = 0.0
        integrator.p[15] = 0.0
        integrator.u[10] = 0.0
    end

    cb1 = DiscreteCallback(voltage_deadband, timer_off)
    cb2 = DiscreteCallback(voltage_outside, timer_on)
    cb3 = DiscreteCallback(timer_hit, TapState)
    cb4 = DiscreteCallback(((u,t,integrator) -> t in tfault[1]), errorState)
    cb5 = DiscreteCallback(((u,t,integrator) -> t in tfault[2]), regularState)
    cb6 = DiscreteCallback(check_voltage, stop_integration)
    cb7 = DiscreteCallback(((u,t,integrator) -> t in t_stop_droop), deactivate_droop)

    sol = solve(problem, Rodas4(), callback = CallbackSet(cb1,cb2,cb3,cb4,cb5,cb6,cb7), tstops=[tfault[1],tfault[2],t_stop_droop], dtmax = dt_max(),progress =true) #
    #sol = AddNaNsIntoSolution(pg,pg_postfault,deepcopy(sol))

    #return PowerGridSolution(sol, pg), event_recorder
    return sol, event_recorder
end

function GetActivePG(fault_state::Bool, postfault_state::Bool)
    if fault_state
        active_pg = 2
    elseif postfault_state
        active_pg = 3
    else
        active_pg = 1
    end
    return active_pg
end

function GetFaultLTVSPG(pg::PowerGrid)
    pg_fault = deepcopy(pg)
    pg_fault.nodes["busv"] = VoltageDependentLoad(P=0.0, Q=0.0, U=1.0, A=0., B=0.,Y_n = complex(1.0/(zfault())))
    return pg_fault
end
function GetPostFaultLTVSPG(pg::PowerGrid)
    nodes_postfault = deepcopy(pg.nodes)
    branches_postfault = deepcopy(pg.lines)
    #delete!(nodes_postfault,"busv")
    #delete!(branches_postfault,"Line_1-v")
    delete!(branches_postfault,"Line_v-2")
    return PowerGrid(nodes_postfault,branches_postfault)
end

function GetTriggCondsLTVS(mtk::ODESystem)
    st = states(mtk)
    @variables t
    s = [
        0.0 ~ t - tfault_on(),
        0.0 ~ t - tfault_off(),
        0.0 ~ t - 5.0,
        0.0 ~ hypot(st[5], st[6]) - 0.99,
    ]
    return s
end

function GetStateResFunLTVS(mtk::ODESystem)
    eqs, aeqs, D_states, A_states = GetSymbolicEquationsAndStates(mtk)
    return [zeros(length(D_states),1) .~ D_states]
end
