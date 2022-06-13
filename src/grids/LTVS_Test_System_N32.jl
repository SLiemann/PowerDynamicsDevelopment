using PowerDynamics
using OrderedCollections: OrderedDict
using Distributed
@everywhere using IfElse
using ModelingToolkit

zfault() = 0.01
tfault_on() = 1.0
tfault_off() = 1.1
dt_max() = 1e-2

function LTVS_Test_System_N32()
    Sbase = 100e6
    Ubase = 400e3
    Ibase = Sbase/Ubase/sqrt(3)
    Zbase = Ubase^2/Sbase
    Ybase = 1.0/Zbase
    #Yshunt = (400e3)^2*1im*100*pi*273.1061e-6*(400/130)^2/Ybase
    Q_Shunt = 1450e6/Sbase
    cos_phi_load = 0.95 #ind
    Pload = -7580e6 /Sbase
    QLoad = Pload * sqrt(1.0/cos_phi_load^2 -1.0)

    buses=OrderedDict(
        "bus1" => SlackAlgebraic(U=1.0),
        "bus2" => VoltageDependentLoad(P=0.0, Q=0.0, U=1.0, A=0., B=0.,Y_n = complex(0.0)),
        "bus3" => VoltageDependentLoad(P=0.0, Q=Q_Shunt, U=1.0, A=1.0, B=0.,Y_n = complex(0.0)),
        "bus4" => GeneralVoltageDependentLoad(P=Pload, Q = QLoad, U=1.0, Ap=0.0, Bp=1.0,Aq = 1.0, Bq= 0.0,Y_n = complex(0.0)),
        "bus5" => SixOrderMarconatoMachineAVROEL(Sbase=Sbase,Srated=5300e6,H = 3, P=4440e6/Sbase, D=0., Ω=50, R_a = 0.0,
                                             T_ds=0.9545455,T_qs=0.3,T_dss=0.0333333,T_qss=0.025,
                                             X_d=2.2,X_q=2.0, X_ds=0.3,X_qs=0.4, X_dss=0.2,
                                             X_qss=0.2,T_AA=0.,V0 = 1.0, Ifdlim = 3.0618,
                                             L1 = -18.0, G1 = 120.0, Ta = 5.0, Tb = 12.5,
                                             G2 = 10.0, L2 = 5.0),
        "busv" => VoltageDependentLoad(P=0.0, Q=0.0, U=1.0, A=0., B=0.,Y_n = complex(0.0)))
    #Lines
    Z_SumLine = (3.140255 + 1im*17.48548)/Zbase
    B_half_SumLine = 1im*100*pi*19.49005*1e-6/2.0*Zbase
    Z_4032_4044 = (9.6 + 1im*80.0)/Zbase
    B_half_4032_4044 = 1im*100*pi*4.77*1e-6/2.0*Zbase
    #Transformers
    Z_Trafo_Netz = 0.13*(400e3)^2/7580e6/Zbase #same base power
    Z_OLTC= 0.1*(130e3)^2/7580e6 * (400/130)^2 / Zbase
    Z_Trafo_SM= 0.15*(130e3)^2/5300e6 * (400/130)^2 / Zbase
    branches=OrderedDict(
        "Line_1-2"=> PiModelLine(from= "bus1", to = "bus2",y=1.0/Z_SumLine, y_shunt_km=B_half_SumLine, y_shunt_mk=B_half_SumLine),
        "Line_1-v"=> PiModelLine(from= "bus1", to = "busv",y=1.0/(Z_4032_4044/2.0), y_shunt_km=B_half_4032_4044, y_shunt_mk=0.0),
        "Line_v-2"=> PiModelLine(from= "bus2", to = "busv",y=1.0/(Z_4032_4044/2.0), y_shunt_km=0.0, y_shunt_mk=B_half_4032_4044),
        "Trafo_Netz"=> StaticPowerTransformer(from="bus2",to="bus3",Sbase=Sbase,Srated=7580e6,uk=0.13,XR_ratio=Inf,
                                           i0=0.0,Pv0=0.0,tap_side = "HV",tap_pos = 7,tap_inc = 1.0),
        "OLTC"=> StaticPowerTransformer(from="bus3",to="bus4",Sbase=Sbase,Srated=7580e6,uk=0.1,XR_ratio=Inf,
                                           i0=0.0,Pv0=0.0,tap_side = "HV",tap_pos = -1,tap_inc = 1.0),
        "Trafo_SM"=> StaticPowerTransformer(from="bus3",to="bus5",Sbase=Sbase,Srated=5300e6,uk=0.15,XR_ratio=Inf,
                                          i0=0.0,Pv0=0.0,tap_side = "HV",tap_pos = 5,tap_inc = 1.0))
        return PowerGrid(buses, branches)
end

function run_LTVS_simulation(pg::PowerGrid,ic1::Array{Float64,1},tspan::Tuple{Float64,Float64};t_stop_droop = Inf)
    tfault = [tfault_on(), tfault_off()]

    pg_fault = GetFaultLTVSPG(pg)
    pg_postfault = GetPostFaultLTVSPG(pg)

    params = GFC_LTVS_params()
    problem = ODEProblem{true}(rhs(pg),ic1,tspan,params)
    timer_start = -1.0
    timer_now   = 0.0
    branch_oltc = "OLTC"
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
        sol1 = integrator.sol
        if tap < tap_max
            tap += 1
            node = StaticPowerTransformer(from=OLTC.from,to=OLTC.to,Srated=OLTC.Srated,
                                          uk=OLTC.uk,XR_ratio=OLTC.XR_ratio,i0=OLTC.i0,
                                          Pv0=OLTC.Pv0,Sbase=OLTC.Sbase,
                                          tap_side = OLTC.tap_side, tap_pos = tap,tap_inc = OLTC.tap_inc)
            if postfault_state
                np_pg = deepcopy(pg_postfault)
            elseif fault_state
                np_pg = deepcopy(pg_fault)
            else
                np_pg = deepcopy(pg)
            end
            np_pg.lines[branch_oltc] = node
            ode = rhs(np_pg)
            op_prob = ODEProblem(ode, sol1[end], (0.0, 1e-6), params, initializealg = BrownFullBasicInit())
            x2 = solve(op_prob,Rodas4())
            x2 = x2.u[end]
            integrator.f = ode
            integrator.cache.tf.f = integrator.f
            integrator.u = x2#sol1[end]
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

    return PowerGridSolution(sol, pg), event_recorder
end
