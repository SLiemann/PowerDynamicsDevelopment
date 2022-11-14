using PowerDynamics
using OrderedCollections: OrderedDict
using Distributed
@everywhere using IfElse

Sbase = 8000e6
Ubase = 400e3
Ibase = Sbase/Ubase/sqrt(3)
Zbase = Ubase^2/Sbase

zfault() = (20+1im*20)/Zbase
tfault_on() = 0.5
tfault_off() = 0.6
dt_max() = 1e-2

function LTVS_Test_System_N32_GFM(;gfm=1,awu=1.0) #1 = droop, 2 = matching, 3 = dVOC, 4 = VSM
    Q_Shunt_EHV = 600e6/Sbase
    Q_Shunt_HV = 850e6/Sbase
    Pload = -7580e6 /Sbase
    QLoad = -2243.7e6/Sbase 
    position_fault = 0.9 #0 at slack 1.0 at bus 2

    Srated = 5300e6
    pref = 4440e6/Sbase
    imax_csa = 1.0
    imax_dc = 1.2
    anti_windup = awu

    Ubase_gfm = 1e3
    Sbase_gfm = 100e6
    Zbase_gfm = Ubase_gfm^2/Sbase_gfm

    Zbase_dc_gfm = (3*1e3*sqrt(2/3))^2/Sbase

    Rdc_pu = 1.2/Zbase_dc_gfm
    Gdc = 1.0/Rdc_pu
    Xcdc = 1.0/(100*pi*0.008*200) /Zbase_dc_gfm
    Cdc = 1.0/Xcdc/100/pi

    R_f = 0.001/200/Zbase_gfm
    L_f = 1*10^-6;
    Xlf = L_f * 100*pi /Zbase_gfm
    C_f = 200*300*10^-6;
    Xcf = 1.0/(100*pi*C_f) /Zbase_gfm

    buses=OrderedDict(
        "bus0" => SlackAlgebraic(U=1.0532), # 
        "bus1" => VoltageDependentLoad(P=0.0, Q=0.0, U=1.0, A=1.0, B=0.0,Y_n = complex(0.0)),
        "bus_ehv" => VoltageDependentLoad(P=0.0, Q=Q_Shunt_EHV, U=1.0, A=1.0, B=0.0,Y_n = complex(0.0)),
        "bus_hv" => VoltageDependentLoad(P=0.0, Q=Q_Shunt_HV,  U=1.0, A=1.0, B=0.0,Y_n = complex(0.0)),
        "bus_load" => GeneralVoltageDependentLoad(P=Pload, Q = QLoad, U=1.0, Ap=0.0, Bp=1.0,Aq = 1.0, Bq= 0.0,Y_n = complex(0.0)),
        "busv" => ThreePhaseFault(Yfault = 1.0/zfault(),p_ind=1))

    if gfm == 1
        buses["bus_gfm"] = droop(Sbase = Sbase,Srated = Srated, p0set = pref, u0set = 1.00,Kp_droop = pi,Kp_uset = 0.001, Ki_uset = 0.5,
                                 Kdc = 100.0, gdc = Gdc, cdc = Cdc, xlf = Xlf, rf = R_f, xcf =  Xcf, Tdc = 0.05, Kp_u = 0.52,
                                 Ki_u = 1.161022, Kp_i = 0.738891, Ki_i = 1.19, imax_csa = imax_csa, imax_dc = imax_dc, p_red = anti_windup, 
                                 p_ind = collect(5:21))
    elseif gfm == 2 
        buses["bus_gfm"] = MatchingControlRed(Sbase = Sbase,Srated = Srated,p0set = pref,u0set = 1.00,Kp_uset = 0.001, Ki_uset = 0.5,
                                              Kdc = 100.0, gdc = Gdc,cdc = Cdc,xlf = Xlf,rf = R_f, xcf =  Xcf,Tdc = 0.05,Kp_u = 0.52,
                                              Ki_u = 1.161022,Kp_i = 0.738891,Ki_i = 1.19,imax_csa = imax_csa,imax_dc = imax_dc,p_red = anti_windup,
                                                p_ind = collect(5:20))
    elseif gfm == 3
        buses["bus_gfm"] = dVOC(Sbase = Sbase,Srated = Srated,p0set = pref, q0set = 0.0,u0set = 1.00,eta = pi,alpha = 0.1*2/3*1000^2,
                                Kdc = 100.0, gdc = Gdc,cdc = Cdc,xlf = Xlf,rf = R_f, xcf =  Xcf,Tdc = 0.05,Kp_u = 0.52,Ki_u = 1.161022,
                                Kp_i = 0.738891,Ki_i = 1.19,imax_csa = imax_csa,imax_dc = imax_dc,p_red = anti_windup,ϵ = 1e-9*0, p_ind = collect(5:21))
    elseif gfm == 4
        buses["bus_gfm"] = VSM(Sbase = Sbase, Srated = Srated, p0set = pref, u0set = 1.00, J = 2, Dp = 100, Kp_uset = 0.001, Ki_uset = 0.5,
                               Kdc = 100.0, gdc = Gdc, cdc = Cdc, xlf = Xlf, rf = R_f, xcf =  Xcf, Tdc = 0.05, Kp_u = 0.52, Ki_u = 1.161022,
                               Kp_i = 0.738891, Ki_i = 1.19, imax_csa = imax_csa, imax_dc = imax_dc, p_red = anti_windup, p_ind = collect(5:22))
    else
        error("wrong number, gfm should be between 1-4")
    end

    #Lines
    Z_SumLine = (3.140255 + 1im*17.48548)/Zbase
    B_half_SumLine = 1im*100*pi*19.49005*1e-6/2.0*Zbase
    Z_4032_4044 = (9.6 + 1im*80.0)/Zbase
    B_half_4032_4044 = 1im*100*pi*4.770001*1e-6/2.0*Zbase
    #Slack internal resistance  
    R1 = 1.514082/Zbase
    X1 = 17.24593/Zbase
    branches=OrderedDict(
        "Line_0-1"=> PiModelLine(from= "bus0", to = "bus1",y=1.0/(R1+1im*X1), y_shunt_km=0.0, y_shunt_mk=0.0),
        "Line_1-2"=> PiModelLine(from= "bus1", to = "bus_ehv",y=1.0/Z_SumLine, y_shunt_km=B_half_SumLine, y_shunt_mk=B_half_SumLine),
        "Line_1-v"=> PiModelLineParam(from= "bus1", to = "busv",y=1.0/(Z_4032_4044*position_fault), y_shunt_km=B_half_4032_4044, y_shunt_mk=0.0,p_ind=2),
        "Line_v-2"=> PiModelLineParam(from= "bus_ehv", to = "busv",y=1.0/(Z_4032_4044*(1.0-position_fault)), y_shunt_km=B_half_4032_4044, y_shunt_mk=0,p_ind=3),
        "Trafo_Netz"=> StaticPowerTransformer(from="bus_ehv",to="bus_hv",Sbase=Sbase,Srated=8000e6,uk=0.12,XR_ratio=Inf,
                                           i0=0.0,Pv0=0.0,tap_side = "HV",tap_pos = 5,tap_inc = 1.0),
        "OLTC"=> StaticPowerTransformerTapParam(from="bus_hv",to="bus_load",Sbase=Sbase,Srated=8000e6,uk=0.11,XR_ratio=Inf,
                                           i0=0.0,Pv0=0.0,tap_side = "LV",tap_pos = 6,tap_inc = 1.0,tap_max=20,tap_min=-20,p_ind=4),
        "Trafo_SM"=> StaticPowerTransformer(from="bus_hv",to="bus_gfm",Sbase=Sbase,Srated=5300e6,uk=0.15,XR_ratio=Inf,
                                          i0=0.0,Pv0=0.0,tap_side = "HV",tap_pos = 5,tap_inc = 1.0))
        return PowerGrid(buses, branches)
end

function GetParamsGFM(pg::PowerGrid)
    node = pg.nodes["bus_gfm"]
    tap = pg.lines["OLTC"].tap_pos
    println(tap)
    params = Vector{Float64}()
    push!(params,0.0) #for 3ph fault, start without fault
    push!(params,1) # for 1st PiModelLineParam
    push!(params,1) # for 2nd PiModelLineParam
    push!(params,0) # for OLTC
    if typeof(node) == droop
        return vcat(params,getallParameters(node)[5:21])
    elseif typeof(node) == MatchingControlRed
        return vcat(params,getallParameters(node)[5:20])
    elseif typeof(node) == dVOC
        return vcat(params,getallParameters(node)[6:22])
    elseif typeof(node) == VSM
        return vcat(params,getallParameters(node)[5:22])
    else
        error("bus_gfm is not a valid GFM node")
    end
end

function run_LTVS_N32_simulation(gfm_choice,awu_choice,tspan::Tuple{Float64,Float64})
    pg = LTVS_Test_System_N32_GFM(gfm=gfm_choice,awu=awu_choice)
    Qmax   = [Inf,Inf, Inf, Inf,Inf,5300/8000*sqrt(1-0.85^2),Inf]
    Qmin   = -Qmax
    U1,δ1,ic0,cu = PowerFlowClassic(pg,iwamoto = false,max_tol = 1e-4,iter_max = 100,Qmax = Qmax, Qmin = Qmin,Qlimit_iter_check=80)
    #println.(keys(pg.nodes) .=> U1)
    #println.(keys(pg.nodes) .=> δ1)
    pg, ic1 = InitializeInternalDynamics(pg,ic0)

    tfault = [tfault_on(), tfault_off()]
    params = GetParamsGFM(pg)
    problem = ODEProblem{true}(rhs(pg),ic1,tspan,params)
    timer_start = -1.0

    branch_oltc = "OLTC"
    index_U_oltc = PowerDynamics.variable_index(pg.nodes,pg.lines[branch_oltc].to,:u_r)
    index_U_load = PowerDynamics.variable_index(pg.nodes,"bus_load",:u_r)

    function TapState(integrator)
        timer_start = integrator.t
        sol1 = integrator.sol
        integrator.p[4] += 1 
        op_prob = ODEProblem(integrator.f, sol1[end], (0.0, 1e-6), integrator.p, initializealg = BrownFullBasicInit())
        x2 = solve(op_prob,Rodas4())
        x2 = x2.u[end]
        integrator.u = deepcopy(x2)
        auto_dt_reset!(integrator)
    end

    function voltage_deadband(u,t,integrator)
         0.99 <= sqrt(u[index_U_oltc]*u[index_U_oltc] + u[index_U_oltc+1]*u[index_U_oltc+1]) <= 1.01
    end

    function timer_off(integrator)
        if timer_start != -1
            timer_start = -1
        end
    end

    function voltage_outside(u,t,integrator)
         sqrt(u[index_U_oltc]*u[index_U_oltc] + u[index_U_oltc+1]*u[index_U_oltc+1]) < 0.99
    end

    function timer_on(integrator)
        if timer_start == -1
            timer_start = integrator.t
        end
    end

    function timer_hit(u,t,integrator)
        if timer_start == -1
            return false
        else
            return t-timer_start > 4.0
        end
    end

    function errorState(integrator)
        integrator.p[1] = 1
        sol1 = integrator.sol
        ode =integrator.f
        op_prob = ODEProblem(ode, sol1[end], (0.0, 1e-6), integrator.p, initializealg = BrownFullBasicInit())
        x2 = solve(op_prob,Rodas4())
        x2 = x2.u[end]
        integrator.u = deepcopy(x2)
        auto_dt_reset!(integrator)
    end

    function regularState(integrator)
        integrator.p[1] = 0.0
        integrator.p[2] = 1e-10
        integrator.p[3] = 1e-10

        sol = integrator.sol
        #ode   = rhs(pg_postfault)
        ic_tmp = deepcopy(integrator.sol.u[indexin(tfault[1],integrator.sol.t)[1]])
        ic_tmp = getPreFaultVoltages(pg,ic_tmp,deepcopy(sol[end]))
        op_prob = ODEProblem(integrator.f, ic_tmp, (0.0, 1e-6), integrator.p, initializealg = BrownFullBasicInit())
        x3 = solve(op_prob,Rodas4())
        x3 = x3.u[end]
        integrator.u = deepcopy(x3)
    
        auto_dt_reset!(integrator)
    end

    function check_voltage(u,t,integrator)
            sqrt(u[index_U_load]*u[index_U_load] + u[index_U_load+1]*u[index_U_load+1]) < 0.3 && t > 5.0
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

    sol = solve(problem, Rodas4(), callback = CallbackSet(cb1,cb2,cb3,cb4,cb5,cb6), tstops=[tfault[1],tfault[2]], dtmax = dt_max(),force_dtmin=false,maxiters=1e5, initializealg = BrownFullBasicInit())
    #sol = AddNaNsIntoSolution(pg,pg_postfault,deepcopy(sol))
    if sol.retcode != :Success
        sol = DiffEqBase.solution_new_retcode(sol, :Success)
    end
    return PowerGridSolution(sol, pg)
end
