using PowerDynamics
using OrderedCollections: OrderedDict
using Distributed
import DiffEqBase: initialize_dae!
@everywhere using IfElse

Sbase = 8000e6
Ubase = 400e3
Ibase = Sbase/Ubase/sqrt(3)
Zbase = Ubase^2/Sbase

zfault() = (20+1im*20)/Zbase
tfault_on() = 0.1
tfault_off() =  tfault_on() + 0.10
dt_max() = 1e-2

function LTVS_Test_System_N32_GENTPJ() #1 = droop, 2 = matching, 3 = dVOC, 4 = VSM
    Q_Shunt_EHV = 600e6/Sbase
    Q_Shunt_HV = 850e6/Sbase
    Pload = -7580e6 /Sbase
    QLoad = -2243.7e6/Sbase 
    position_fault = 0.9 #0 at slack, 1.0 at bus 2

    Srated = 5300e6 #5150e6 for LTVS, 5125e6
    pref = 4440e6/Sbase
 

    buses=OrderedDict(
        "bus0" => SlackAlgebraic(U=1.054275078250000), # 1.0532 for 100% I , 1.0495 for 100% Z 1.054080675
        "bus1" => VoltageDependentLoad(P=0.0, Q=0.0, U=1.0, A=1.0, B=0.0,Y_n = complex(0.0)),
        "bus_ehv" => VoltageDependentLoad(P=0.0, Q=Q_Shunt_EHV, U=1.0, A=1.0, B=0.0,Y_n = complex(0.0)),
        "bus_hv" => VoltageDependentLoad(P=0.0, Q=Q_Shunt_HV,  U=1.0, A=1.0, B=0.0,Y_n = complex(0.0)),
        "bus_load" => GeneralVoltageDependentLoad(P=Pload, Q = QLoad, U=1.0, Ap=0.01, Bp=0.99,Aq = 1.0, Bq= 0.00,Y_n = complex(0.0)),
        "busv" => ThreePhaseFault(rfault=8e3,xfault=8e3,p_ind=collect(1:2)),
        "bus_sm" => gentpjAVROEL(Sbase=Sbase,Srated=Srated, H=6.0, P=pref, D=0.0, Ω=50, R_a=0, T_d0s=7.0, T_q0s=1.5, T_d0ss=0.05,
        T_q0ss=0.05, X_d=2.2, X_q=2.0, X_ds=0.3, X_qs=0.4, X_dss=0.2, X_qss=0.2, X_l=0.15, S_10=0.1, S_12=0.3,K_is=0.0,
        V0 = 1.0, Ifdlim = 3.0618, L1 = -20.0, G1 = 120.0, Ta = 5.0, Tb = 12.5, G2 = 10.0, L2 = 5.0))

        #buses["bus_gfm"] = droop(Sbase = Sbase,Srated = Srated, p0set = pref, u0set = 1.00,Kp_droop = pi,Kp_uset = 0.001, Ki_uset = 0.5,
        #                         Kdc = 100.0, gdc = Gdc, cdc = Cdc, xlf = Xlf, rf = R_f, xcf =  Xcf, Tdc = 0.05, Kp_u = 0.52,
        #                         Ki_u = 1.161022, Kp_i = 0.738891, Ki_i = 1.19, imax_csa = imax_csa, imax_dc = imax_dc, p_red = anti_windup, LVRT_on = 0.0,
        #                         p_ind = collect(6:25))

    #Lines
    Z_SumLine = (3.140255 + 1im*17.48548)/Zbase
    B_half_SumLine = 1im*100*pi*19.49005*1e-6/2.0*Zbase
    Z_4032_4044 = (9.6 + 1im*80.0)/Zbase
    B_half_4032_4044 = 1im*100*pi*4.770001*1e-6/2.0*Zbase
    #Slack internal resistance  
    R1 = 1.514081970099058/Zbase
    X1 = 17.245934327062923/Zbase
    branches=OrderedDict(
        "Line_0-1"=> PiModelLine(from= "bus0", to = "bus1",y=1.0/(R1+1im*X1), y_shunt_km=0.0, y_shunt_mk=0.0),
        "Line_1-2"=> PiModelLine(from= "bus1", to = "bus_ehv",y=1.0/Z_SumLine, y_shunt_km=B_half_SumLine, y_shunt_mk=B_half_SumLine),
        "Line_1-v"=> PiModelLineParam(from= "bus1", to = "busv",y=1.0/(Z_4032_4044*position_fault), y_shunt_km=B_half_4032_4044, y_shunt_mk=0.0,p_ind=3),
        "Line_v-2"=> PiModelLineParam(from= "bus_ehv", to = "busv",y=1.0/(Z_4032_4044*(1.0-position_fault)), y_shunt_km=B_half_4032_4044, y_shunt_mk=0,p_ind=4),
        "Trafo_Netz"=> StaticPowerTransformer(from="bus_ehv",to="bus_hv",Sbase=Sbase,Srated=8000e6,uk=0.12,XR_ratio=Inf,
                                           i0=0.0,Pv0=0.0,tap_side = "HV",tap_pos = 5,tap_inc = 1.0),
        "OLTC"=> StaticPowerTransformerTapParam(from="bus_hv",to="bus_load",Sbase=Sbase,Srated=8000e6,uk=0.11,XR_ratio=Inf,
                                           i0=0.0,Pv0=0.0,tap_side = "LV",tap_pos = 6,tap_inc = 1.0,tap_max=20,tap_min=-20,p_ind=5),
        "Trafo_SM"=> StaticPowerTransformer(from="bus_hv",to="bus_sm",Sbase=Sbase,Srated=5300e6,uk=0.15,XR_ratio=Inf,
                                          i0=0.0,Pv0=0.0,tap_side = "HV",tap_pos = 5,tap_inc = 1.0))
        return PowerGrid(buses, branches)
end

function GetParamsGENTPJ(pg::PowerGrid)
    tap = pg.lines["OLTC"].tap_pos
    params = Vector{Float64}()
    push!(params,8e3) #for 3ph fault, start without fault
    push!(params,8e3) #for 3ph fault, start without fault
    push!(params,1) # for 1st PiModelLineParam
    push!(params,1) # for 2nd PiModelLineParam
    push!(params,0) # for OLTC, Δtap_pos
    return params
end

function GetFaultLTVSPG(pg::PowerGrid)
    pg_fault = deepcopy(pg)
    pg_fault.nodes["busv"] = VoltageDependentLoad(P=0.0, Q=0.0, U=1.0, A=0., B=0.,Y_n = complex(1.0/(zfault())))
    return pg_fault
end

function GetContFaultPG(pg::PowerGrid)
    pg_fault = deepcopy(pg)
    pg_fault.nodes["busv"] = ThreePhaseFaultContinouos(rfault=8e3,xfault=8e3,Tf=1e-3/2,p_ind=collect(1:2))
    return pg_fault
end

function GetPostFaultLTVSPG(pg::PowerGrid)
    nodes_postfault = deepcopy(pg.nodes)
    branches_postfault = deepcopy(pg.lines)
    delete!(branches_postfault,"Line_v-2")
    branches_postfault["Line_1-v"] = PiModelLineParam(from= "bus_ehv", to = "busv",y=1.0/(1*(1.0-0.1)), y_shunt_km=0.0, y_shunt_mk=0,p_ind=3)
    return PowerGrid(nodes_postfault,branches_postfault)
end

function Initialize_N32_GENTPJ()
    pg = LTVS_Test_System_N32_GENTPJ()
    Qmax   = [Inf,Inf, Inf, Inf,Inf,5300/8000*sqrt(1-0.8377^2),Inf]
    Qmin   = -Qmax
    U1,δ1,ic0,cu = PowerFlowClassic(pg,iwamoto = false,max_tol = 1e-4,iter_max = 100,Qmax = Qmax, Qmin = Qmin,Qlimit_iter_check=80)
    pg, ic = InitializeInternalDynamics(pg,ic0)
    display(U1.=>δ1)
    return pg,ic
end

function MakeODEProblemN32(pg::PowerGrid,ic::Vector{Float64},params::Vector{Float64},tspan::Tuple{Float64,Float64})
    return ODEProblem{true}(rhs(pg),ic,tspan,params)
end

function run_LTVS_N32_simulation(tspan::Tuple{Float64,Float64})
    pg,ic = Initialize_N32_GENTPJ()
    return simulate_LTVS_N32_simulation_GENTPJ(pg,ic,tspan,zfault())
end

function simulate_LTVS_N32_simulation_GENTPJ(pg::PowerGrid,ic::Vector{Float64},tspan::Tuple{Float64,Float64},zfault::Union{Float64,Complex{Float64},Complex{Int64}})
    pg_postfault = GetPostFaultLTVSPG(pg)
    params = GetParamsGENTPJ(pg)
    problem= ODEProblem{true}(rhs(pg),ic,tspan,params)
    tfault = [tfault_on(), tfault_off()]
    timer_start = -1.0
    FRT = 1.0 # -1 for LVRT, 0 for HVRT, 1.0 for stable
    tap_dir = 1
    rfault = real(zfault) <= 0.0 ? 1e-5 : real(zfault)
    xfault = imag(zfault) <= 0.0 ? 1e-5 : imag(zfault)

    branch_oltc = "OLTC"
    index_U_oltc = PowerDynamics.variable_index(pg.nodes,pg.lines[branch_oltc].to,:u_r)
    index_U_sm = PowerDynamics.variable_index(pg.nodes,"bus_sm",:u_r)
    index_U_load = PowerDynamics.variable_index(pg.nodes,"bus_load",:u_r)

    function TapState(integrator)
        timer_start = integrator.t
        sol1 = integrator.sol
        integrator.p[5] += 1*tap_dir 
        initialize_dae!(integrator,BrownFullBasicInit())
        #op_prob = ODEProblem(integrator.f, sol1[end], (0.0, 1e-6), integrator.p, initializealg = BrownFullBasicInit())
        #x2 = solve(op_prob,Rodas4())
        #x2 = x2.u[end]
        #integrator.u = deepcopy(x2)
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

    function voltage_outside_low(u,t,integrator)
         sqrt(u[index_U_oltc]*u[index_U_oltc] + u[index_U_oltc+1]*u[index_U_oltc+1]) < 0.99
    end

    function voltage_outside_high(u,t,integrator)
        sqrt(u[index_U_oltc]*u[index_U_oltc] + u[index_U_oltc+1]*u[index_U_oltc+1]) > 1.01
   end

    function timer_on_low(integrator)
        tap_dir = 1
        if timer_start == -1
            timer_start = integrator.t
        end
    end

    function timer_on_high(integrator)
        tap_dir = -1
        if timer_start == -1
            timer_start = integrator.t
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
        integrator.p[1] = rfault
        integrator.p[2] = xfault
        
        pg_cfault = GetContFaultPG(pg);
        ic_init= deepcopy(integrator.sol[end])
        len = length(symbolsof(pg.nodes["bus_sm"]))
        # insert two extra states for continouos fault
        ic_tmp = vcat(ic_init[1:end-len],[pg_cfault.nodes["busv"].rfault,pg_cfault.nodes["busv"].xfault],ic_init[end-len+1:end])
        # create problem and simulate for 10ms
        op_prob = ODEProblem(rhs(pg_cfault), ic_tmp, (0.0, 0.010), integrator.p)
        x2 = solve(op_prob,Rodas4(),dtmax=1e-4,initializealg = BrownFullBasicInit(),alg_hints=:stiff,verbose=false,abstol=1e-8,reltol=1e-8)

        ic_end = x2.u[end]
        # delete states of continouos fault
        ic_end = vcat(ic_end[1:end-len-2],ic_end[end-len+1:end])
        # change only algebraic states of original problem
        ind_as = findall(x-> iszero(x),diag(integrator.f.mass_matrix))
        #ind_as = getVoltageSymbolPositions(pg)
        for i in ind_as
            ic_init[i] = ic_end[i]
        end
        
        if x2.retcode == :Success
            integrator.u = deepcopy(ic_init)
            initialize_dae!(integrator,BrownFullBasicInit())
            auto_dt_reset!(integrator)
        elseif x2.retcode == :Unstable
            println("Terminated at $(integrator.t) due to detected instability at fault initialisation.")
            terminate!(integrator)
        else
            ind = getVoltageSymbolPositions(pg)
            ic_tmp = deepcopy(integrator.sol[end])
            for i in ind
                ic_tmp[i] = ic_tmp[i]/4.0
            end
            initialize_dae!(integrator,BrownFullBasicInit())
            auto_dt_reset!(integrator)
        end
    end

    function regularState(integrator)
        integrator.p[1] = 8e3 #fault is zero again
        integrator.p[2] = 8e3 #fault is zero again

        #First create continouos fault and then post-fault grid
        pg_pcfault = GetContFaultPG(pg);
        #pg_pcfault = GetPostFaultLTVSPG(pg_pcfault);
        ic_init= deepcopy(integrator.sol[end])
        len = length(symbolsof(pg.nodes["bus_sm"]))
        # insert two extra states for continouos fault
        ic_tmp = vcat(ic_init[1:end-len],[rfault,xfault],ic_init[end-len+1:end])
        # create problem and simulate for 2 ms
        op_prob = ODEProblem(rhs(pg_pcfault), ic_tmp, (0.0, 0.002), integrator.p)
        x2 = solve(op_prob,Rodas4(),dtmax=1e-4,initializealg = BrownFullBasicInit(),alg_hints=:stiff,verbose=false,abstol=1e-8,reltol=1e-8)
        
        tmp_retcode = deepcopy(x2.retcode)

        ode   = rhs(pg_postfault)
        integrator.f = ode
        integrator.cache.tf.f = integrator.f
        

        if tmp_retcode == :Success
            ic_end = x2.u[end]
            # delete states of continouos fault
            ic_end = vcat(ic_end[1:end-len-2],ic_end[end-len+1:end])
            # change only algebraic states of original problem
            ind_as = findall(x-> iszero(x),diag(integrator.f.mass_matrix))
            for i in ind_as
                ic_init[i] = ic_end[i]
            end
            integrator.u = deepcopy(ic_init)
            initialize_dae!(integrator,BrownFullBasicInit())
            auto_dt_reset!(integrator)
        elseif tmp_retcode == :Unstable
            println("Terminated at $(integrator.t) due to detected instability at post-fault initialisation.")
            terminate!(integrator)
        else
            ic_tmp = deepcopy(integrator.sol.u[indexin(tfault[1],integrator.sol.t)[1]])
            #integrator.u  = getPreFaultVoltages(pg,ic_tmp,deepcopy(sol[end]))
            integrator.u = getPreFaultAlgebraicStates(pg,ic_tmp,deepcopy(integrator.sol[end]))
            initialize_dae!(integrator,BrownFullBasicInit())
            auto_dt_reset!(integrator)
        end

        #=sol = integrator.sol
        ode   = rhs(pg_postfault)
        ic_tmp = deepcopy(integrator.sol.u[indexin(tfault[1],integrator.sol.t)[1]])
        #ic_tmp = getPreFaultVoltages(pg,ic_tmp,deepcopy(sol[end]))
        ic_tmp = getPreFaultAlgebraicStates(pg,ic_tmp,deepcopy(sol[end]))
        op_prob = ODEProblem(ode, ic_tmp, (0.0, 1e-6), integrator.p, initializealg = BrownFullBasicInit())
        x3 = solve(op_prob,Rodas5(),initializealg = BrownFullBasicInit(),alg_hints=:stiff,abstol=1e-8,reltol=1e-8)
        x3 = x3.u[end]
        integrator.u = deepcopy(x3)
        integrator.f = ode
        integrator.cache.tf.f = integrator.f
        auto_dt_reset!(integrator)=#
    end

    function check_OLTC_voltage(u,t,integrator)
            sqrt(u[index_U_load]*u[index_U_load] + u[index_U_load+1]*u[index_U_load+1]) < 0.3 && t > 5.0
    end

    function stop_integration(integrator)
        println("Terminated at $(integrator.t)")
        terminate!(integrator)
        #necessary, otherwise PowerGridSolution throws error
        integrator.sol = DiffEqBase.solution_new_retcode(integrator.sol, :Success)
    end

    function stop_integration_LVRT(integrator)
        FRT = -1.0
        println("Terminated at $(integrator.t) due to LVRT")
        terminate!(integrator)
        integrator.sol = DiffEqBase.solution_new_retcode(integrator.sol, :Success)
    end

    function stop_integration_HVRT(integrator)
        FRT = 0.0
        println("Terminated at $(integrator.t) due to HVRT")
        terminate!(integrator)
        integrator.sol = DiffEqBase.solution_new_retcode(integrator.sol, :Success)
    end

    function start_LVRT(integrator)
        integrator.p[end] = 0.85/2.85
    end
    
    function end_LVRT(integrator)
        integrator.p[end] = 0.0
    end

    #function check_GENTPJ_voltage_LVRT(u,t,integrator)
    #    sqrt(u[index_U_sm]*u[index_U_sm] + u[index_U_sm+1]*u[index_U_sm+1]) < integrator.u[end]
    #end

    #function check_GENTPJ_voltage_HVRT13(u,t,integrator)
    #    sqrt(u[index_U_sm]*u[index_U_sm] + u[index_U_sm+1]*u[index_U_sm+1]) > 1.3  && t < tfault[1]+0.1
    #end

    #function check_GENTPJ_voltage_HVRT12(u,t,integrator)
    #    sqrt(u[index_U_sm]*u[index_U_sm] + u[index_U_sm+1]*u[index_U_sm+1]) > 1.2  && t > tfault[1]+0.1
    #end

    cb1 = DiscreteCallback(voltage_deadband, timer_off)
    cb2 = DiscreteCallback(voltage_outside_low, timer_on_low)
    cb21 = DiscreteCallback(voltage_outside_high, timer_on_high)
    cb3 = DiscreteCallback(timer_hit, TapState)
    cb4 = PresetTimeCallback(tfault[1], errorState)
    cb5 = PresetTimeCallback(tfault[2], regularState)
    #cb6 = DiscreteCallback(check_OLTC_voltage, stop_integration)
    #cb7 = PresetTimeCallback(tfault[1]+0.15, start_LVRT)
    #cb8 = PresetTimeCallback(tfault[1]+3.0,end_LVRT)
    #cb9 = DiscreteCallback(check_GENTPJ_voltage_LVRT, stop_integration_LVRT)
    #cb10 = DiscreteCallback(check_GENTPJ_voltage_HVRT13, stop_integration_HVRT)
    #cb11 = DiscreteCallback(check_GENTPJ_voltage_HVRT12, stop_integration_HVRT)

    stiff  = repeat([:stiff],length(ic))
    sol = solve(problem, Rodas5(autodiff=true), callback = CallbackSet(cb1,cb2,cb21,cb3,cb4,cb5), tstops=[tfault[1],tfault[2]], dtmax = dt_max(),force_dtmin=false,maxiters=1e6, initializealg = BrownFullBasicInit(),alg_hints=:stiff,abstol=1e-8,reltol=1e-8) #
    # good values abstol=1e-8,reltol=1e-8 and Rodas5(autodiff=true) for droop
    success = deepcopy(sol.retcode)
    if sol.retcode != :Success
        sol = DiffEqBase.solution_new_retcode(sol, :Success)
    end
    return PowerGridSolution(sol, pg), success, FRT
end
