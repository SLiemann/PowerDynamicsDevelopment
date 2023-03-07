using PowerDynamics
using OrderedCollections: OrderedDict
using Distributed
@everywhere using IfElse
using ModelingToolkit

Sbase = 8000e6
Ubase = 400e3
Ibase = Sbase/Ubase/sqrt(3)
Zbase = Ubase^2/Sbase

zfault() = (20+1im*20)/Zbase
tfault_on() = 1.0
tfault_off() = 1.1
dt_max() = 1e-3

function LTVS_Test_System_N32(;tap=6)
    Ybase = 1.0/Zbase
    #Yshunt = (400e3)^2*1im*100*pi*273.1061e-6*(400/130)^2/Ybase
    Q_Shunt_EHV = 600e6/Sbase
    Q_Shunt_HV = 850e6/Sbase
    cos_phi_load = 0.95 #ind
    Pload = -7580e6 /Sbase
    QLoad = -2243.7e6/Sbase #Pload * sqrt(1.0/cos_phi_load^2 -1.0)
    position_fault = 0.9 #0 at slack 1.0 at bus 2

    buses=OrderedDict(
        "bus0" => SlackAlgebraic(U=1.054080675),# 1.05717
        "bus1" => VoltageDependentLoad(P=0.0, Q=0.0, U=1.0, A=1.0, B=0.0,Y_n = complex(0.0)),
        "bus_ehv" => VoltageDependentLoad(P=0.0, Q=Q_Shunt_EHV, U=1.0, A=1.0, B=0.0,Y_n = complex(0.0)),
        "bus_hv" => VoltageDependentLoad(P=0.0, Q=Q_Shunt_HV,  U=1.0, A=1.0, B=0.0,Y_n = complex(0.0)),
        "bus_load" => GeneralVoltageDependentLoad(P=Pload, Q = QLoad, U=1.0, Ap=0.0, Bp=1.0,Aq = 1.0, Bq= 0.0,Y_n = complex(0.0)),
 

        
        "bus_sm" => gentpjAVROEL(Sbase=Sbase,Srated=5300e6, H=6.0, P=4440e6/Sbase, D=0.0, Ω=50, R_a=0, T_d0s=7.0, T_q0s=1.5, T_d0ss=0.05,
                                 T_q0ss=0.05, X_d=2.2, X_q=2.0, X_ds=0.3, X_qs=0.4, X_dss=0.2, X_qss=0.2, X_l=0.15, S_10=0.1, S_12=0.3,K_is=0.0,
                                 V0 = 1.0, Ifdlim = 3.0618, L1 = -20.0, G1 = 120.0, Ta = 5.0, Tb = 12.5, G2 = 10.0, L2 = 5.0),
       "busv" => VoltageDependentLoad(P=0.0, Q=0.0, U=1.0, A=0., B=0.,Y_n = complex(0.0)))
    #Lines
    Z_SumLine = (3.140255 + 1im*17.48548)/Zbase
    B_half_SumLine = 1im*100*pi*19.49005*1e-6/2.0*Zbase
    Z_4032_4044 = (9.6 + 1im*80.0)/Zbase
    B_half_4032_4044 = 1im*100*pi*4.770001*1e-6/2.0*Zbase
    #Transformers
    Z_Trafo_Netz = 0.13*(400e3)^2/7580e6/Zbase #same base power
    Z_OLTC= 0.1*(130e3)^2/7580e6 * (400/130)^2 / Zbase
    Z_Trafo_SM= 0.15*(130e3)^2/5300e6 * (400/130)^2 / Zbase
    R1 = 1.514082/Zbase
    X1 = 17.24593/Zbase
    branches=OrderedDict(
        "Line_0-1"=> PiModelLine(from= "bus0", to = "bus1",y=1.0/(R1+1im*X1), y_shunt_km=0.0, y_shunt_mk=0.0),
        "Line_1-2"=> PiModelLine(from= "bus1", to = "bus_ehv",y=1.0/Z_SumLine, y_shunt_km=B_half_SumLine, y_shunt_mk=B_half_SumLine),
        "Line_1-v"=> PiModelLine(from= "bus1", to = "busv",y=1.0/(Z_4032_4044*position_fault), y_shunt_km=B_half_4032_4044, y_shunt_mk=0.0),
        "Line_v-2"=> PiModelLine(from= "bus_ehv", to = "busv",y=1.0/(Z_4032_4044*(1.0-position_fault)), y_shunt_km=B_half_4032_4044, y_shunt_mk=0),
        "Trafo_Netz"=> StaticPowerTransformer(from="bus_ehv",to="bus_hv",Sbase=Sbase,Srated=8000e6,uk=0.12,XR_ratio=Inf,
                                           i0=0.0,Pv0=0.0,tap_side = "HV",tap_pos = 5,tap_inc = 1.0),
        "OLTC"=> StaticPowerTransformer(from="bus_hv",to="bus_load",Sbase=Sbase,Srated=8000e6,uk=0.11,XR_ratio=Inf,
                                           i0=0.0,Pv0=0.0,tap_side = "LV",tap_pos = tap,tap_inc = 1.0),
        "Trafo_SM"=> StaticPowerTransformer(from="bus_hv",to="bus_sm",Sbase=Sbase,Srated=5300e6,uk=0.15,XR_ratio=Inf,
                                          i0=0.0,Pv0=0.0,tap_side = "HV",tap_pos = 5,tap_inc = 1.0))
        return PowerGrid(buses, branches)
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
    branches_postfault["Line_1-v"] = PiModelLine(from= "bus_ehv", to = "busv",y=1.0/(1*(1.0-0.1)), y_shunt_km=0.0, y_shunt_mk=0)
    return PowerGrid(nodes_postfault,branches_postfault)
end

function run_LTVS_N32_simulation(pg::PowerGrid,ic1::Array{Float64,1},tspan::Tuple{Float64,Float64})
    tfault = [tfault_on(), tfault_off()]

    pg_fault = GetFaultLTVSPG(pg)
    pg_postfault = GetPostFaultLTVSPG(pg)

    params = []#GFC_LTVS_params()
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
    index_U_load = PowerDynamics.variable_index(pg.nodes,"bus_load",:u_r)
    #event_recorder = Array{Float64,2}(undef,0,4+length(params))
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
            #active_pg = GetActivePG(fault_state,postfault_state)
            #event_recorder = vcat(event_recorder,[integrator.t active_pg integrator.p' 3 1])
        end
    end

    function voltage_deadband(u,t,integrator)
         0.99 <= sqrt(u[index_U_oltc]*u[index_U_oltc] + u[index_U_oltc+1]*u[index_U_oltc+1]) <= 1.01
    end

    function timer_off(integrator)
        if timer_start != -1
            timer_start = -1
            #active_pg = GetActivePG(fault_state,postfault_state)
            #event_recorder = vcat(event_recorder,[integrator.t active_pg integrator.p' 4 1])
        end
    end

    function voltage_outside(u,t,integrator)
         sqrt(u[index_U_oltc]*u[index_U_oltc] + u[index_U_oltc+1]*u[index_U_oltc+1]) < 0.99
    end

    function timer_on(integrator)
        if timer_start == -1
            timer_start = integrator.t
            #active_pg = GetActivePG(fault_state,postfault_state)
            #event_recorder = vcat(event_recorder,[integrator.t active_pg integrator.p' 4 1])
        end
    end

    function timer_hit(u,t,integrator)
        if timer_start == -1
            return false
        else
            return t-timer_start > 10.0
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
        #active_pg = GetActivePG(fault_state,postfault_state)
        #event_recorder = vcat(event_recorder,[integrator.t active_pg integrator.p' 1 1])
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
        index_U_load = PowerDynamics.variable_index(pg_postfault.nodes,"bus_load",:u_r)
        postfault_state = true
        fault_state = false
        #active_pg = GetActivePG(fault_state,postfault_state)
        #event_recorder = vcat(event_recorder,[integrator.t active_pg integrator.p' 2 1])
    end

    function check_voltage(u,t,integrator)
            sqrt(u[index_U_load]*u[index_U_load] + u[index_U_load+1]*u[index_U_load+1]) < 0.3
    end

    function stop_integration(integrator)
        println("Terminated at $(integrator.t)")
        terminate!(integrator)
        #necessary, otherwise PowerGridSolution throws error
        integrator.sol = DiffEqBase.solution_new_retcode(integrator.sol, :Success)
    end

    cb1 = DiscreteCallback(voltage_deadband, timer_off)
    cb2 = DiscreteCallback(voltage_outside, timer_on)
    cb3 = DiscreteCallback(timer_hit, TapState)
    cb4 = DiscreteCallback(((u,t,integrator) -> t in tfault[1]), errorState)
    cb5 = DiscreteCallback(((u,t,integrator) -> t in tfault[2]), regularState)
    cb6 = DiscreteCallback(check_voltage, stop_integration)

    sol = solve(problem, Rodas4(), callback = CallbackSet(cb1,cb2,cb3,cb4,cb5,cb6), tstops=[tfault[1],tfault[2]], dtmax = dt_max(),progress =true)
    #sol = AddNaNsIntoSolution(pg,pg_postfault,deepcopy(sol))

    return PowerGridSolution(sol, pg)#, event_recorder
end
