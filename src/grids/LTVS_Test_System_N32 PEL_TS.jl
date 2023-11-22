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
dt_max() = 1e-3

function LTVS_Test_System_N32_PEL_TS()
    Q_Shunt_EHV = 600e6/Sbase
    Q_Shunt_HV = 850e6/Sbase
    Pload = -7580e6 /Sbase
    QLoad = -2243.7e6/Sbase 
    position_fault = 0.9 #0 at slack, 1.0 at bus 2

    Srated = 5300e6 #5150e6 for LTVS, 5125e6
    pref = 4440e6/Sbase

    vbase = 230.0
    sbase = 1000.0
    zbase = vbase^2/sbase
    C = 700e-6;
    ω0 = 100*pi
    #xcpu = (1/(ω0*C))/zbase
    xcpu = 0.036
    Cpu = 1/(xcpu*ω0)
 
    buses=OrderedDict(
        "bus0" => SlackAlgebraic(U=1.054275078250000), # 1.0532 for 100% I , 1.0495 for 100% Z 1.054080675
        "bus1" => VoltageDependentLoad(P=0.0, Q=0.0, U=1.0, A=1.0, B=0.0,Y_n = complex(0.0)),
        "bus_ehv" => VoltageDependentLoad(P=0.0, Q=Q_Shunt_EHV, U=1.0, A=1.0, B=0.0,Y_n = complex(0.0)),
        "bus_hv" => VoltageDependentLoad(P=0.0, Q=Q_Shunt_HV,  U=1.0, A=1.0, B=0.0,Y_n = complex(0.0)),
        "bus_load" => nPFC(Cd=Cpu, Pdc = -Pload, p_ind=collect(6:7)),
        "busv" => ThreePhaseFault(rfault=8e3,xfault=8e3,p_ind=collect(1:2)),
        "bus_sm" => gentpjAVROEL(Sbase=Sbase,Srated=Srated, H=6.0, P=pref, D=0.0, Ω=50, R_a=0, T_d0s=7.0, T_q0s=1.5, T_d0ss=0.05,
        T_q0ss=0.05, X_d=2.2, X_q=2.0, X_ds=0.3, X_qs=0.4, X_dss=0.2, X_qss=0.2, X_l=0.15, S_10=0.1, S_12=0.3,K_is=0.0,
        V0 = 1.0, Ifdlim = 3.0618, L1 = -20.0, G1 = 120.0, Ta = 5.0, Tb = 12.5, G2 = 10.0, L2 = 5.0))

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
                                           i0=0.0,Pv0=0.0,tap_side = "HV",tap_pos = 5*0-2,tap_inc = 1.0),
        "OLTC"=> StaticPowerTransformerTapParam(from="bus_hv",to="bus_load",Sbase=Sbase,Srated=8000e6,uk=0.11,XR_ratio=Inf,
                                           i0=0.0,Pv0=0.0,tap_side = "LV",tap_pos = 6*0-5,tap_inc = 1.0,tap_max=20,tap_min=-20,p_ind=5),
        "Trafo_SM"=> StaticPowerTransformer(from="bus_hv",to="bus_sm",Sbase=Sbase,Srated=5300e6,uk=0.15,XR_ratio=Inf,
                                          i0=0.0,Pv0=0.0,tap_side = "HV",tap_pos = 5*0,tap_inc = 1.0))
        return PowerGrid(buses, branches)
end


function GetParams_PEL_TS(pg::PowerGrid)
    Cd = pg.nodes["bus_load"].Cd
    Pdc = pg.nodes["bus_load"].Pdc
    tap = pg.lines["OLTC"].tap_pos
    params = Vector{Float64}()
    push!(params,8e3) #for 3ph fault, start without fault
    push!(params,8e3) #for 3ph fault, start without fault
    push!(params,1) # for 1st PiModelLineParam
    push!(params,1) # for 2nd PiModelLineParam
    push!(params,0) # for OLTC, Δtap_pos
    push!(params,Cd) # for nPFC, 
    push!(params,Pdc) # for nPFC, 
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

function GetPostFaultLTVSPG_TS(pg::PowerGrid)
    nodes_postfault = deepcopy(pg.nodes)
    branches_postfault = deepcopy(pg.lines)
    delete!(branches_postfault,"Line_v-2")
    branches_postfault["Line_1-v"] = PiModelLineParam(from= "bus_ehv", to = "busv",y=1.0/(1*(1.0-0.1)), y_shunt_km=0.0, y_shunt_mk=0,p_ind=3)
    return PowerGrid(nodes_postfault,branches_postfault)
end

function Initialize_N32_PEL_TS()
    pg = LTVS_Test_System_N32_PEL_TS()
    Qmax   = [Inf,Inf, Inf, Inf,Inf,5300/8000*sqrt(1-0.8377^2),Inf]
    Qmin   = -Qmax
    U1,δ1,ic0,cu = PowerFlowClassic(pg,iwamoto = true,max_tol = 1e-5,iter_max = 100,Qmax = Qmax, Qmin = Qmin,Qlimit_iter_check=80)
    pg, ic = InitializeInternalDynamics(pg,ic0)
    display(U1.=>δ1)
    return pg,ic
end

function MakeODEProblemN32(pg::PowerGrid,ic::Vector{Float64},params::Vector{Float64},tspan::Tuple{Float64,Float64})
    return ODEProblem{true}(rhs(pg),ic,tspan,params)
end

function run_LTVS_N32_simulation_PEL_TS(gfm_choice,awu_choice,tspan::Tuple{Float64,Float64})
    pg,ic = Initialize_N32_PEL_TS()
    return simulate_LTVS_N32_simulation(pg,ic,tspan,zfault())
end

function simulate_LTVS_N32_simulation_PEL_TS(pg::PowerGrid,ic::Vector{Float64},tspan::Tuple{Float64,Float64},zfault::Union{Float64,Complex{Float64},Complex{Int64}})
    pg_postfault = GetPostFaultLTVSPG_TS(pg)
    params = GetParams_PEL_TS(pg)
    problem= ODEProblem{true}(rhs(pg),ic,tspan,params)
    tfault = [tfault_on(), tfault_off()]
    timer_start = -1.0
    FRT = 1.0 # -1 for LVRT, 0 for HVRT, 1.0 for stable
    tap_dir = 1
    rfault = real(zfault) <= 0.0 ? 1e-5 : real(zfault)
    xfault = imag(zfault) <= 0.0 ? 1e-5 : imag(zfault)

    # ODESystems = 1=prefault and fault, 2=postfault
    # wenn bei h und s "null" steht, steht dies für ein Ereignis welches zwar Parameter ändert,
    # aber eigentlich kein richtiges event ist! 
    evr = Array{Float64}(undef,0,4+length(params))
    ind_odesys = 1;

    branch_oltc = "OLTC"
    index_U_oltc = PowerDynamics.variable_index(pg.nodes,pg.lines[branch_oltc].to,:u_r)
    #index_LVRT_gfm = PowerDynamics.variable_index(pg.nodes,"bus_gfm",:LVRT)
    index_U_load = PowerDynamics.variable_index(pg.nodes,"bus_load",:u_r)
    index_vofft2_load = PowerDynamics.variable_index(pg.nodes,"bus_load",:vofft2)
    index_tsum_load = PowerDynamics.variable_index(pg.nodes,"bus_load",:tsum)
    index_ton_load = PowerDynamics.variable_index(pg.nodes,"bus_load",:ton)
    index_toff_load = PowerDynamics.variable_index(pg.nodes,"bus_load",:toff)
    index_p_load = PowerDynamics.variable_index(pg.nodes,"bus_load",:p1)
    index_q_load = PowerDynamics.variable_index(pg.nodes,"bus_load",:q1)
    p_pel_ind = pg.nodes["bus_load"].p_ind

    function TapState(integrator)
        timer_start = integrator.t
        sol1 = integrator.sol
        integrator.p[5] += 1*tap_dir 
        initialize_dae!(integrator,BrownFullBasicInit())
        auto_dt_reset!(integrator)
        evr = vcat(evr, [integrator.t ind_odesys 0 0 integrator.p'])
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

    ### PEL Callbacks START ###
    function f_tsum(u,t,integrator)
        true
    end
    
    function affect_tsum(integrator)
        integrator.u[index_tsum_load]  = integrator.u[index_tsum_load] + integrator.t - integrator.tprev  
    end

    function new_half_wave(u,t,integrator)
        mod(integrator.t/0.01,1.0)
    end

    function affect_new_half_wave(integrator)
        display("NEW HALF WAVE")
        Vabs = hypot(integrator.u[index_U_load],integrator.u[index_U_load+1])*sqrt(2)
        Cd = deepcopy(integrator.p[p_pel_ind[1]])
        Pdc = deepcopy(integrator.p[p_pel_ind[1]])
        T = 0.02
        ton = deepcopy(integrator.u[index_ton_l])

        if ton < 0.0 # Location A
            integrator.u[index_tsum_load]  = 0.0
            vofft2 = deepcopy(integrator.u[index_vofft2_load])
            voffT2_new = CalfnPFCVoffT2(vofft2,Pdc,Cd,(T/2-0.0))
            integrator.u[index_vofft2_load]  = voffT2_new
        else # Location A
            integrator.u[index_tsum_load]  = 0.0
            voff = Vabs*sin(ω0*toff)
            voffT2_new = CalfnPFCVoffT2(voff,Pdc,Cd,(T/2-toff))
            integrator.u[index_vofft2_load]  = voffT2_new
        end
    end

    function ton_pos(u,t,integrator)
        integrator[index_ton_load] >= 0.0
    end

    function ton_neg(u,t,integrator)
        integrator[index_ton_load] < 0.0
    end

    function affect_ton_neg(integrator)
        Vabs = hypot(integrator.u[index_U_load],integrator.u[index_U_load+1])*sqrt(2)
        Cd = deepcopy(integrator.p[p_pel_ind[1]])
        Pdc = deepcopy(integrator.p[p_pel_ind[2]])
        T = 0.02

        vofft2 = deepcopy(integrator.u[index_vofft2_load])
        t = deepcopy(integrator.t)

        if iszero(mod(t/0.02,1.0)) #Location A
            #display("A ton neg")
            #nothing
        elseif  tsum < toff && tsum < ton #Location B
           # display("Location B neg")
            toff_new = CalcnPFCtoff(Vabs,Pdc,Cd)
            ton_new = CalfnPFCton(Vabs,Pdc,Cd,vofft2)  

            integrator.u[index_toff_load]  = toff_new
            integrator.u[index_ton_load]  = ton_new
        elseif  tsum < toff && tsum >= ton
            #display("Location C neg")
            toff_new = CalcnPFCtoff(Vabs,Pdc,Cd)
            integrator.u[index_toff_load]  = toff_new
        elseif  tsum >= toff && tsum < ton
            #display("Location D  neg")
            #nothing
        else
            error("state tsum is not within length of half-cycle: $(tsum)")
        end 
        integrator.u[index_p_load]  = 0.0
        integrator.u[index_q_load]  = 0.0
        initialize_dae!(integrator,BrownFullBasicInit())
        auto_dt_reset!(integrator)       
    end

    function affect_ton_pos(integrator)
        Vabs = hypot(integrator.u[index_U_load],integrator.u[index_U_load+1])*sqrt(2)
        Vabs_prev = hypot(integrator.uprev[index_U_load],integrator.uprev[index_U_load+1])*sqrt(2)
        Cd = deepcopy(integrator.p[p_pel_ind[1]])
        Pdc = deepcopy(integrator.p[p_pel_ind[2]])
        T = 0.02

        vofft2 = deepcopy(integrator.u[index_vofft2_load])
        tsum = deepcopy(integrator.u[index_tsum_load])  
        ton = deepcopy(integrator.u[index_ton_load])
        toff = deepcopy(integrator.u[index_toff_load]) 
        t = deepcopy(integrator.t)

        #display(integrator.u[index_vofft2_load:index_q_load])

        if abs(Vabs - Vabs_prev) > 1e-3 
            if iszero(mod(t/0.02,1.0)) #Location A
                #display("A ton pos")
            elseif  tsum < toff && tsum < ton #Location B
                #display("B ton pos")
                p1_new, q1_new = CalcnPFCP1Q1(Vabs,Pdc,Cd,ton,toff)
                toff_new = CalcnPFCtoff(Vabs,Pdc,Cd)
                ton_new = CalfnPFCton(Vabs,Pdc,Cd,vofft2)  
                

                integrator.u[index_toff_load]  = toff_new
                integrator.u[index_ton_load]  = ton_new
                integrator.u[index_p_load]  = p1_new
                integrator.u[index_q_load]  = q1_new

                display(integrator.u[index_vofft2_load:index_q_load])

            elseif  tsum < toff && tsum >= ton
                #display("C ton pos")
                p1_new, q1_new = CalcnPFCP1Q1(Vabs,Pdc,Cd,ton,toff)
                toff_new = CalcnPFCtoff(Vabs,Pdc,Cd)
                
                integrator.u[index_toff_load]  = toff_new
                integrator.u[index_p_load]  = p1_new
                integrator.u[index_q_load]  = q1_new
            elseif  tsum >= toff && tsum < ton
                #display("D ton pos")
                p1_new, q1_new = CalcnPFCP1Q1(Vabs,Pdc,Cd,ton,toff)

                integrator.u[index_p_load]  = p1_new
                integrator.u[index_q_load]  = q1_new
            else
                error("state tsum is not within length of half-cycle: $(tsum)")
            end
        end
        initialize_dae!(integrator,BrownFullBasicInit())
        auto_dt_reset!(integrator)       
    end

    ### PEL Callbacks END ###

    function errorState(integrator)
        integrator.p[1] = rfault
        integrator.p[2] = xfault
        ind_odesys = 1;
        
        # Den evr würde man nur benötigen, wenn qimax und error NICHT  zusammenfallsen 
        #evr = vcat(evr, [integrator.t ind_odesys 0 0 integrator.p'])
        
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
        ind_odesys = 2;

        #First create continouos fault and then post-fault grid
        pg_pcfault = GetContFaultPG(pg);
        #pg_pcfault = GetPostFaultLTVSPG_TS(pg_pcfault);
        ic_init= deepcopy(integrator.sol[end])
        len = length(symbolsof(pg.nodes["bus_sm"]))
        # insert two extra states for continouos fault
        ic_tmp = vcat(ic_init[1:end-len],[rfault,xfault],ic_init[end-len+1:end])
        # create problem and simulate for 2 ms
        #ic_tmp[34] = 0.0
        #ic_tmp[35] = 0.0
        
        op_prob = ODEProblem(rhs(pg_pcfault), ic_tmp, (0.0, 0.02), integrator.p)
        x2 = solve(op_prob,Rodas4(),dtmax=1e-4,initializealg = BrownFullBasicInit(),alg_hints=:stiff,verbose=false,abstol=1e-8,reltol=1e-8)
        # pgsol_tmp = PowerGridSolution(x2,pg_pcfault)
        # plotallvoltages(pgsol_tmp);
        # plot(myplot(pgsol_tmp,"bus_gfm",:q_imax))

        tmp_retcode = deepcopy(x2.retcode)

        ode   = rhs(pg_postfault)
        integrator.f = ode
        integrator.cache.tf.f = integrator.f
        
        # integrator.u = deepcopy(ic_init)
        # initialize_dae!(integrator,BrownFullBasicInit())
        # auto_dt_reset!(integrator)

        if tmp_retcode == ReturnCode.Success
            ic_end = x2.u[end]
            # delete states of continouos fault
            ic_end = vcat(ic_end[1:end-len-2],ic_end[end-len+1:end])
            # change only algebraic states of original problem
            ind_as = findall(x-> iszero(x),diag(integrator.f.mass_matrix))
            for i in ind_as
                ic_init[i] = ic_end[i]
            end

            #ic_init[32] = 0.0
            #ic_init[33] = 0.0
            #evr = vcat(evr, [integrator.t ind_odesys 2 2 integrator.p'])

            integrator.u = deepcopy(ic_init)
            initialize_dae!(integrator,BrownFullBasicInit())
            auto_dt_reset!(integrator)
        elseif tmp_retcode == ReturnCode.Unstable
            println("Terminated at $(integrator.t) due to detected instability at post-fault initialisation.")
            terminate!(integrator)
        else
            ic_tmp = deepcopy(integrator.sol.u[indexin(tfault[1],integrator.sol.t)[1]])
            #integrator.u  = getPreFaultVoltages(pg,ic_tmp,deepcopy(sol[end]))
            integrator.u = getPreFaultAlgebraicStates(pg,ic_tmp,deepcopy(integrator.sol[end]))
            initialize_dae!(integrator,BrownFullBasicInit())
            auto_dt_reset!(integrator)
        end
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

    # function check_GFM_voltage_LVRT(u,t,integrator)
    #     sqrt(u[index_U_gfm]*u[index_U_gfm] + u[index_U_gfm+1]*u[index_U_gfm+1]) < integrator.u[index_LVRT_gfm]
    # end

    # function check_GFM_voltage_HVRT13(u,t,integrator)
    #     sqrt(u[index_U_gfm]*u[index_U_gfm] + u[index_U_gfm+1]*u[index_U_gfm+1]) > 1.3  && t < tfault[1]+0.1
    # end

    # function check_GFM_voltage_HVRT12(u,t,integrator)
    #     sqrt(u[index_U_gfm]*u[index_U_gfm] + u[index_U_gfm+1]*u[index_U_gfm+1]) > 1.2  && t > tfault[1]+0.1
    # end


    cb1 = DiscreteCallback(voltage_deadband, timer_off) 
    cb2 = DiscreteCallback(voltage_outside_low, timer_on_low)
    cb21 = DiscreteCallback(voltage_outside_high, timer_on_high)
    cb3 = DiscreteCallback(timer_hit, TapState)
    cb4 = PresetTimeCallback(tfault[1], errorState)
    cb5 = PresetTimeCallback(tfault[2], regularState)
    #cb7 = PresetTimeCallback(tfault[1]+0.15, start_LVRT)
    #cb8 = PresetTimeCallback(tfault[1]+3.0,end_LVRT)
    #cb9 = DiscreteCallback(check_GFM_voltage_LVRT, stop_integration_LVRT)
    #cb10 = DiscreteCallback(check_GFM_voltage_HVRT13, stop_integration_HVRT)
    #cb11 = DiscreteCallback(check_GFM_voltage_HVRT12, stop_integration_HVRT)
    #cb12 = DiscreteCallback(imax_on, affect_imax_on, save_positions=(false,true))

    cb_nhw =ContinuousCallback(new_half_wave, affect_new_half_wave, save_positions=(true,true))
    cb_tpos =DiscreteCallback(ton_pos, affect_ton_pos, save_positions=(true,true))
    cb_tneg =DiscreteCallback(ton_neg, affect_ton_neg, save_positions=(true,true))
    cb_tsum =DiscreteCallback(f_tsum, affect_tsum, save_positions=(true,true))



    stiff  = repeat([:stiff],length(ic)) #, callback = CallbackSet(cb1,cb2,cb21,cb3,cb4,cb5)
    #cb7,cb8,cb9,cb10,cb11,
    tstops_sim =[tfault[1];tfault[2];collect(tspan[1]:0.02:tspan[2])]
    sort!(tstops_sim)

    sol = solve(problem, Rodas4(autodiff=true),callback = CallbackSet(cb4,cb5,cb_nhw,cb_tpos,cb_tneg,cb_tsum), tstops=tstops_sim, dtmax = dt_max(),force_dtmin=false,maxiters=1e6, initializealg = BrownFullBasicInit(),alg_hints=:stiff,abstol=1e-8,reltol=1e-8) #
    # good values abstol=1e-8,reltol=1e-8 and Rodas5(autodiff=true) for droop
    #success = deepcopy(sol.retcode)
    if sol.retcode != :Success
        sol = DiffEqBase.solution_new_retcode(sol, :Success)
    end
    return PowerGridSolution(sol, pg), evr
end