using PowerDynamics
using OrderedCollections: OrderedDict
using Distributed
import DiffEqBase: initialize_dae!
using SciMLSensitivity
using ForwardDiff
@everywhere using IfElse

Sbase = 8000e6
Ubase = 400e3
Ibase = Sbase/Ubase/sqrt(3)
Zbase = Ubase^2/Sbase

zfault() = (20+1im*20)/Zbase
tfault_on() = 0.1
tfault_off() =  tfault_on() + 0.1
dt_max() = 1e-3

function LTVS_Test_System_N32_GFM_PEL_TS(;gfm=1,awu=1.0,share_pe = 0.3) #1 = droop, 2 = matching, 3 = dVOC, 4 = VSM
    Q_Shunt_EHV = 600e6/Sbase
    Q_Shunt_HV = 850e6/Sbase
    Pload = -7580e6 /Sbase
    QLoad = -2243.7e6/Sbase 
    position_fault = 0.9 #0 at slack, 1.0 at bus 2

    Srated = 5300e6 #5150e6 for LTVS 
    pref = 4440e6/Sbase
    imax_csa = 1.0
    imax_dc = 1.2
    anti_windup = awu

    Ubase_gfm = 1e3
    Sbase_gfm = 100e6
    Zbase_gfm = Ubase_gfm^2/Sbase_gfm

    Zbase_dc_gfm = (3*1e3*sqrt(2/3))^2/Sbase_gfm

    Rdc_pu = 1.2/Zbase_dc_gfm
    Gdc = 1.0/Rdc_pu
    Xcdc = 1.0/(100*pi*0.008*200) /Zbase_dc_gfm
    Cdc = 1.0/Xcdc/100/pi

    R_f = 0.001/200/Zbase_gfm
    L_f = 1*10^-6;
    Xlf = L_f * 100*pi /Zbase_gfm
    C_f = 200*300*10^-6;
    Xcf = 1.0/(100*pi*C_f) /Zbase_gfm

    vbase = 230.0
    sbase = 1000.0
    zbase = vbase^2/sbase
    C = 700e-6;
    ω0 = 100*pi
    #xcpu = (1/(ω0*C))/zbase
    xcpu = 0.036 
    Cpu = 1/(xcpu*ω0) * share_pe

    p_static, q_static = CalcnPFCPower(0.99034*sqrt(2),-Pload*share_pe,Cpu) #0.9904
    qoff = QLoad - q_static
    poff = Pload - p_static

    buses=OrderedDict(
        "bus0" => SlackAlgebraic(U=1.054275078250000), # 1.0532 for 100% I , 1.0495 for 100% Z 1.054080675
        "bus1" => VoltageDependentLoad(P=0.0, Q=0.0, U=1.0, A=1.0, B=0.0,Y_n = complex(0.0)),
        "bus_ehv" => VoltageDependentLoad(P=0.0, Q=Q_Shunt_EHV, U=1.0, A=1.0, B=0.0,Y_n = complex(0.0)),
        "bus_hv" => VoltageDependentLoad(P=0.0, Q=Q_Shunt_HV,  U=1.0, A=1.0, B=0.0,Y_n = complex(0.0)),
        "bus_load" => nPFC(Cd=Cpu, Pdc = -p_static, p_offset = poff, q_offset=qoff, p_ind=collect(25:26)),
        "busv" => ThreePhaseFault(rfault=8e3,xfault=8e3,p_ind=collect(1:2)))

    if gfm == 1
        buses["bus_gfm"] = droopTS(Sbase = Sbase,Srated = Srated, p0set = pref, u0set = 1.00,Kp_droop = pi,Kp_uset = 0.001, Ki_uset = 0.5,
                                 Kdc = 100.0, gdc = Gdc, cdc = Cdc, xlf = Xlf, rf = R_f, xcf =  Xcf, Tdc = 0.05, Kp_u = 0.52,
                                 Ki_u = 1.161022, Kp_i = 0.738891, Ki_i = 1.19, imax_csa = imax_csa, imax_dc = imax_dc, LVRT_on = 0.0,
                                 p_ind = collect(6:24))
    elseif gfm == 2 
        buses["bus_gfm"] = MatchingControlRed(Sbase = Sbase,Srated = Srated,p0set = pref,u0set = 1.00,Kp_uset = 0.001, Ki_uset = 0.5,
                                              Kdc = 100.0, gdc = Gdc,cdc = Cdc,xlf = Xlf,rf = R_f, xcf =  Xcf,Tdc = 0.05,Kp_u = 0.52,
                                              Ki_u = 1.161022,Kp_i = 0.738891,Ki_i = 1.19,imax_csa = imax_csa,imax_dc = imax_dc,p_red = anti_windup,LVRT_on = 0.0, p_ind = collect(6:22))
    elseif gfm == 3
        buses["bus_gfm"] = dVOC(Sbase = Sbase,Srated = Srated,p0set = pref, q0set = 0.0,u0set = 1.00,eta = pi,alpha = 0.1*2/3*1000^2,
                                Kdc = 100.0, gdc = Gdc,cdc = Cdc,xlf = Xlf,rf = R_f, xcf =  Xcf,Tdc = 0.05,Kp_u = 0.52,Ki_u = 1.161022,
                                Kp_i = 0.738891,Ki_i = 1.19,imax_csa = imax_csa,imax_dc = imax_dc,p_red = anti_windup,ϵ = 1e-9*0,LVRT_on = 0.0, p_ind = collect(6:23))
    elseif gfm == 4
        buses["bus_gfm"] = VSM(Sbase = Sbase, Srated = Srated, p0set = pref, u0set = 1.00, J = 2, Dp = 100, Kp_uset = 0.001, Ki_uset = 0.5,
                               Kdc = 100.0, gdc = Gdc, cdc = Cdc, xlf = Xlf, rf = R_f, xcf =  Xcf, Tdc = 0.05, Kp_u = 0.52, Ki_u = 1.161022,
                               Kp_i = 0.738891, Ki_i = 1.19, imax_csa = imax_csa, imax_dc = imax_dc, p_red = anti_windup,LVRT_on = 0.0, p_ind = collect(6:24))
    else
        error("wrong number, gfm should be between 1-4")
    end

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
        "Trafo_SM"=> StaticPowerTransformer(from="bus_hv",to="bus_gfm",Sbase=Sbase,Srated=5300e6,uk=0.15,XR_ratio=Inf,
                                          i0=0.0,Pv0=0.0,tap_side = "HV",tap_pos = 5,tap_inc = 1.0))
        return PowerGrid(buses, branches)
end

function GetParamsGFM_GFM_PEL_TS(pg::PowerGrid)
    Cd = pg.nodes["bus_load"].Cd
    Pdc = pg.nodes["bus_load"].Pdc
    node = pg.nodes["bus_gfm"]
    tap = pg.lines["OLTC"].tap_pos
    params = Vector{Float64}()
    push!(params,8e3) #for 3ph fault, start without fault
    push!(params,8e3) #for 3ph fault, start without fault
    push!(params,1) # for 1st PiModelLineParam
    push!(params,1) # for 2nd PiModelLineParam
    push!(params,0) # for OLTC, Δtap_pos
    #push!(params,0) # for OLTC, for timer
    if typeof(node) == droopTS
        params = vcat(params,getallParameters(node)[3:21])
    elseif typeof(node) == MatchingControlRed
        params = vcat(params,getallParameters(node)[5:21])
    elseif typeof(node) == dVOC
        params = vcat(params,getallParameters(node)[6:23])
    elseif typeof(node) == VSM
        params = vcat(params,getallParameters(node)[5:23])
    else
        error("bus_gfm is not a valid GFM node")
    end
    push!(params,Cd) # for nPFC, 
    push!(params,Pdc) # for nPFC, 
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

function GetPostFaultLTVSPG_TS(pg::PowerGrid)
    nodes_postfault = deepcopy(pg.nodes)
    branches_postfault = deepcopy(pg.lines)
    delete!(branches_postfault,"Line_v-2")
    branches_postfault["Line_1-v"] = PiModelLineParam(from= "bus_ehv", to = "busv",y=1.0/(1*(1.0-0.1)), y_shunt_km=0.0, y_shunt_mk=0,p_ind=3)
    return PowerGrid(nodes_postfault,branches_postfault)
end

function Initialize_N32_GFM_PEL_TS()
    pg = LTVS_Test_System_N32_GFM_PEL_TS()
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

function simulate_LTVS_N32_simulation_N32_GFM_PEL_TS(pg::PowerGrid,ic::Vector{Float64},tspan::Tuple{Float64,Float64},zfault::Union{Float64,Complex{Float64},Complex{Int64}})
    pg_postfault = GetPostFaultLTVSPG_TS(pg)
    params = GetParamsGFM_GFM_PEL_TS(pg)
    #params[8] += 0.5*params[8] 
    #problem= ODEProblem{true}(rhs(pg),ic,tspan,params)
    sens_prob = ODEForwardSensitivityProblem(rhs(pg), ic, tspan, params,ForwardDiffSensitivity(convert_tspan=true);)
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
    index_U_gfm = PowerDynamics.variable_index(pg.nodes,"bus_gfm",:u_r)
    index_LVRT_gfm = PowerDynamics.variable_index(pg.nodes,"bus_gfm",:LVRT)
    index_iset_abs_gfm = PowerDynamics.variable_index(pg.nodes,"bus_gfm",:iset_abs)
    index_qimax_gfm = PowerDynamics.variable_index(pg.nodes,"bus_gfm",:q_imax)
    index_qidcmax_gfm = PowerDynamics.variable_index(pg.nodes,"bus_gfm",:q_idcmax)
    index_U_load = PowerDynamics.variable_index(pg.nodes,"bus_load",:u_r)
    index_vofft2_load = PowerDynamics.variable_index(pg.nodes,"bus_load",:vofft2)
    index_ton_load = PowerDynamics.variable_index(pg.nodes,"bus_load",:ton)
    index_toff_load = PowerDynamics.variable_index(pg.nodes,"bus_load",:toff)
    index_voff_load = PowerDynamics.variable_index(pg.nodes,"bus_load",:voff)
    index_qon_load = PowerDynamics.variable_index(pg.nodes,"bus_load",:q_on) 
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
        Vabs = hypot(integrator.u[index_U_load],integrator.u[index_U_load+1])*sqrt(2)
        Cd = deepcopy(integrator.p[p_pel_ind[1]])
        Pdc = deepcopy(integrator.p[p_pel_ind[2]])
        T = 0.02

        voffT2 = deepcopy(integrator.u[index_vofft2_load])
        ton = deepcopy(integrator.u[index_ton_load])
        toff = deepcopy(integrator.u[index_toff_load]) 
        voff = deepcopy(integrator.u[index_voff_load])
        tsum = mod(Float32(integrator.t/0.01),1.0)/100.0

        if iszero(tsum) &&  !(tfault[1]<= integrator.t < tfault[1]+0.01) && !(tfault[2]<= integrator.t < tfault[2]+0.01)
            #display("pel: $(integrator.t)")
            #display("pel: $(ton)")
            if ton >= 0.0
                #voff = Vabs*sin(100*pi*toff)
                voffT2 = CalfnPFCVoffT2(voff,Pdc,Cd,(T/2-toff))
            else
                voff = deepcopy(integrator.u[index_vofft2_load])
                voffT2 = CalfnPFCVoffT2(voff,Pdc,Cd,(T/2-0.0))
            end
            #display(integrator.t)
            ton = CalfnPFCton(Vabs,Pdc,Cd,voffT2)
            if ton >= 0
                toff = CalcnPFCtoff(Vabs,Pdc,Cd)
                voff = Vabs*sin(100*pi*toff)
            end 
            #display(voff)
            #display("pel: $(ton)")
            integrator.u[index_vofft2_load]  = voffT2
        elseif tsum < toff && tsum < ton
            #display("B ton pos")
            ton = CalfnPFCton(Vabs,Pdc,Cd,voffT2)
            if ton >= 0
                toff = CalcnPFCtoff(Vabs,Pdc,Cd)
                voff = Vabs*sin(100*pi*toff)
            end 
        elseif  tsum < toff && tsum >= ton
            #display("C ton pos")
            if ton >= 0
                toff = CalcnPFCtoff(Vabs,Pdc,Cd)
                voff = Vabs*sin(100*pi*toff)
            end 
        end

        if ton >= 0.0
            integrator.u[index_qon_load] = 1.0
        else
            integrator.u[index_qon_load] = 0.0           
        end 

        integrator.u[index_voff_load]  = voff
        integrator.u[index_toff_load]  = toff
        integrator.u[index_ton_load]  = ton
        initialize_dae!(integrator,BrownFullBasicInit())
    end

    ### PEL Callbacks END ###

    function errorState(integrator)
        integrator.p[1] = rfault
        integrator.p[2] = xfault
        ind_odesys = 1;
        # Den evr würde man nur benötigen, wenn qimax und error NICHT  zusammenfallsen 
        #evr = vcat(evr, [integrator.t ind_odesys 0 0 integrator.p'])
        
        pg_cfault = GetContFaultPG(pg);
        integrator.u[index_qon_load] = 0.0
        initialize_dae!(integrator,BrownFullBasicInit())  
        ic_init= deepcopy(integrator.sol[end])
        len = length(symbolsof(pg.nodes["bus_gfm"]))
        # insert two extra states for continouos fault
        ic_tmp = vcat(ic_init[1:end-len],[pg_cfault.nodes["busv"].rfault,pg_cfault.nodes["busv"].xfault],ic_init[end-len+1:end])
        # create problem and simulate for 10ms
        #ic_tmp[41] = 1.0
        #ic_tmp[42] = 0.0
        op_prob = ODEProblem(rhs(pg_cfault), ic_tmp, (0.0, 0.1), integrator.p)
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

        if ic_end[index_iset_abs_gfm] > 1.0
            display(ic_end[index_iset_abs_gfm+2])
            integrator.u = ic_end
            integrator.u[index_qimax_gfm] = 1.0
            initialize_dae!(integrator,BrownFullBasicInit())
        else
            integrator.u = ic_init
        end

         ## Init PEL Model
         Vabs = hypot(ic_init[index_U_load],ic_init[index_U_load+1])*sqrt(2)
         voff = integrator[index_voff_load]
         Cd = deepcopy(integrator.p[p_pel_ind[1]])
         Pdc = deepcopy(integrator.p[p_pel_ind[2]])
         toff = integrator[index_toff_load]
         T = 0.02
         ω0 = 100*pi
 
         VoffT2 = deepcopy(integrator.u[index_vofft2_load])
         #voff_new = Vabs*sin(ω0*toff)
         #VoffT2_new = CalfnPFCVoffT2(voff_new,Pdc,Cd,(T/2-toff))
         #toff_new = CalcnPFCtoff(Vabs,Pdc,Cd)
         #voff_new = voff*sin(100*pi*toff_new)
         ton_new = CalfnPFCton(Vabs,Pdc,Cd,VoffT2)  
 
         if ton_new >= 0.0
             integrator.u[index_qon_load] = 1.0
         else
             integrator.u[index_qon_load] = 0.0           
         end 
         #display(ton_new)
         #integrator.u[index_vofft2_load]  = VoffT2_new
         #integrator.u[index_voff_load]  = voff_new
         #integrator.u[index_toff_load]  = toff_new
         integrator.u[index_ton_load]  = ton_new
         initialize_dae!(integrator,BrownFullBasicInit())
        # if x2.retcode == :Success
        #     integrator.u = deepcopy(ic_init)
        #     initialize_dae!(integrator,BrownFullBasicInit())
        #     auto_dt_reset!(integrator)
        # elseif x2.retcode == :Unstable
        #     println("Terminated at $(integrator.t) due to detected instability at fault initialisation.")
        #     terminate!(integrator)
        # else
        #     ind = getVoltageSymbolPositions(pg)
        #     ic_tmp = deepcopy(integrator.sol[end])
        #     for i in ind
        #         ic_tmp[i] = ic_tmp[i]/4.0
        #     end
        #     initialize_dae!(integrator,BrownFullBasicInit())
        #     auto_dt_reset!(integrator)
        # end
    end

    function regularState(integrator)
        integrator.p[1] = 8e3 #fault is zero again
        integrator.p[2] = 8e3 #fault is zero again
        ind_odesys = 2;
        display("regular")
        #First create continouos fault and then post-fault grid
        pg_pcfault = GetContFaultPG(pg);
        #pg_pcfault = GetPostFaultLTVSPG_TS(pg_pcfault);
        ic_init= deepcopy(integrator.sol[end])
        len = length(symbolsof(pg.nodes["bus_gfm"]))
        # insert two extra states for continouos fault
        ic_tmp = vcat(ic_init[1:end-len],[rfault,xfault],ic_init[end-len+1:end])
        # create problem and simulate for 2 ms
        ic_tmp[index_qimax_gfm] = 0.0
        ic_tmp[index_qidcmax_gfm] = 0.0
        op_prob = ODEProblem(rhs(pg_pcfault), ic_tmp, (0.0, 0.02), integrator.p)
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

        if ic_end[index_iset_abs_gfm] > 1.0
            display(ic_end[index_iset_abs_gfm+2])
            integrator.u = ic_end
            integrator.u[index_qimax_gfm] = 1.0
            initialize_dae!(integrator,BrownFullBasicInit())
        else
            integrator.u = ic_init
        end

        ode   = rhs(pg_postfault)
        integrator.f = ode
        integrator.cache.tf.f = integrator.f
        
        ## Init PEL Model
        ic_end = x2.u[end]
        # delete states of continouos fault
        ic_end = vcat(ic_end[1:end-len-2],ic_end[end-len+1:end])
        # change only algebraic states of original problem
        ind_as = findall(x-> iszero(x),diag(integrator.f.mass_matrix))
        for i in ind_as
            ic_init[i] = ic_end[i]
        end

        Vabs = hypot(ic_init[index_U_load],ic_init[index_U_load+1])*sqrt(2)
        voff = integrator[index_voff_load]
        Cd = deepcopy(integrator.p[p_pel_ind[1]])
        Pdc = deepcopy(integrator.p[p_pel_ind[2]])
        toff = integrator[index_toff_load]
        T = 0.02
        ω0 = 100*pi

        VoffT2 = deepcopy(integrator.u[index_vofft2_load])

        #voff_new = Vabs*sin(ω0*toff)
        #VoffT2_new = CalfnPFCVoffT2(voff_new,Pdc,Cd,(T/2-toff))
        #toff_new = CalcnPFCtoff(Vabs,Pdc,Cd)
        #voff_new = Vabs*sin(100*pi*toff_new)
        ton_new = CalfnPFCton(Vabs,Pdc,Cd,VoffT2)  
        if ton_new >= 0.0
            integrator.u[index_qon_load] = 1.0
        else
            integrator.u[index_qon_load] = 0.0           
        end 

        #integrator.u[index_vofft2_load]  = VoffT2_new
        #integrator.u[index_voff_load]  = voff_new
        #integrator.u[index_toff_load]  = toff_new
        integrator.u[index_ton_load]  = ton_new
        ictmp = deepcopy(integrator.u)
        initialize_dae!(integrator,BrownFullBasicInit())
        #for (ind,val) in enumerate(ictmp)
        #    display([rhs(pg).syms[ind]  ictmp[ind] integrator.u[ind]])
        #end
        # if tmp_retcode == ReturnCode.Success
        #     ic_end = x2.u[end]
        #     # delete states of continouos fault
        #     ic_end = vcat(ic_end[1:end-len-2],ic_end[end-len+1:end])
        #     # change only algebraic states of original problem
        #     ind_as = findall(x-> iszero(x),diag(integrator.f.mass_matrix))
        #     for i in ind_as
        #         ic_init[i] = ic_end[i]
        #     end

        #     ic_init[32] = 0.0
        #     ic_init[33] = 0.0
        #     evr = vcat(evr, [integrator.t ind_odesys 2 2 integrator.p'])

        #     integrator.u = deepcopy(ic_init)
        #     initialize_dae!(integrator,BrownFullBasicInit())
        #     auto_dt_reset!(integrator)
        # elseif tmp_retcode == ReturnCode.Unstable
        #     println("Terminated at $(integrator.t) due to detected instability at post-fault initialisation.")
        #     terminate!(integrator)
        # else
        #     ic_tmp = deepcopy(integrator.sol.u[indexin(tfault[1],integrator.sol.t)[1]])
        #     #integrator.u  = getPreFaultVoltages(pg,ic_tmp,deepcopy(sol[end]))
        #     integrator.u = getPreFaultAlgebraicStates(pg,ic_tmp,deepcopy(integrator.sol[end]))
        #     initialize_dae!(integrator,BrownFullBasicInit())
        #     auto_dt_reset!(integrator)
        # end
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

    function check_GFM_voltage_LVRT(u,t,integrator)
        sqrt(u[index_U_gfm]*u[index_U_gfm] + u[index_U_gfm+1]*u[index_U_gfm+1]) < integrator.u[index_LVRT_gfm]
    end

    function check_GFM_voltage_HVRT13(u,t,integrator)
        sqrt(u[index_U_gfm]*u[index_U_gfm] + u[index_U_gfm+1]*u[index_U_gfm+1]) > 1.3  && t < tfault[1]+0.1
    end

    function check_GFM_voltage_HVRT12(u,t,integrator)
        sqrt(u[index_U_gfm]*u[index_U_gfm] + u[index_U_gfm+1]*u[index_U_gfm+1]) > 1.2  && t > tfault[1]+0.1
    end

    ind_iset_abs = PowerDynamics.variable_index(pg.nodes,"bus_gfm",:iset_abs)
    ind_idc0 = PowerDynamics.variable_index(pg.nodes,"bus_gfm",:idc0)
    ind_qimax = PowerDynamics.variable_index(pg.nodes,"bus_gfm",:q_imax)
    ind_qidcmax = PowerDynamics.variable_index(pg.nodes,"bus_gfm",:q_idcmax)
    imax_csa_tmp = pg.nodes["bus_gfm"].imax_csa
    imax_dc_tmp = pg.nodes["bus_gfm"].imax_dc

    function imax_on(u,t,integrator)
        u[ind_iset_abs] - imax_csa_tmp - 10 * u[ind_qimax] >=0
    end
    function imax_off(u,t,integrator)
        (u[ind_iset_abs] - imax_csa_tmp) * u[ind_qimax] < 0
    end

    function idcmax_on(u,t,integrator)
        abs(u[ind_idc0]) - imax_dc_tmp - 10 * u[ind_qidcmax] >=0
    end

    function idcmax_off(u,t,integrator)
        (abs(u[ind_idc0]) - imax_dc_tmp) * u[ind_qidcmax] < 0
    end

    function affect_imax_on(integrator)
        evr = vcat(evr, [integrator.t ind_odesys 1 1 integrator.p'])
        integrator.u[ind_qimax] = 1.0
        initialize_dae!(integrator,BrownFullBasicInit())
        auto_dt_reset!(integrator)
    end

    function affect_imax_off(integrator)
        evr = vcat(evr, [integrator.t ind_odesys 2 2 integrator.p'])
        integrator.u[ind_qimax] = 0.0
        initialize_dae!(integrator,BrownFullBasicInit())
        auto_dt_reset!(integrator)
    end

    function affect_idcmax_on(integrator)
        evr = vcat(evr, [integrator.t ind_odesys 3 3 integrator.p'])
        integrator.u[ind_qidcmax] = 1.0
        initialize_dae!(integrator,BrownFullBasicInit())
        auto_dt_reset!(integrator)
    end

    function affect_idcmax_off(integrator)
        evr = vcat(evr, [integrator.t ind_odesys 4 4 integrator.p'])
        integrator.u[ind_qidcmax] = 0.0
        initialize_dae!(integrator,BrownFullBasicInit())
        auto_dt_reset!(integrator)
    end

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
    cb12 = DiscreteCallback(imax_on, affect_imax_on, save_positions=(false,true))
    cb13 = DiscreteCallback(imax_off, affect_imax_off, save_positions=(false,true))
    cb14 = DiscreteCallback(idcmax_on, affect_idcmax_on, save_positions=(false,true))
    cb15 = DiscreteCallback(idcmax_off, affect_idcmax_off, save_positions=(false,true))
    cb_tsum = DiscreteCallback(f_tsum, affect_tsum, save_positions=(false,false))

    tstops_sim =collect(tspan[1]:0.01:tspan[2]);
    sort!(tstops_sim)
    sol = solve(sens_prob, Rodas4(autodiff=true), callback = CallbackSet(cb1,cb2,cb21,cb3,cb4,cb5,cb12,cb13,cb14,cb15,cb_tsum), tstops=tstops_sim, dtmax = dt_max(),force_dtmin=false,maxiters=1e6, initializealg = BrownFullBasicInit(),alg_hints=:stiff,abstol=1e-8,reltol=1e-8) #
    # good values abstol=1e-8,reltol=1e-8 and Rodas5(autodiff=true) for droop
    #success = deepcopy(sol.retcode)
    #if sol.retcode != :Success
    #    sol = DiffEqBase.solution_new_retcode(sol, :Success)
    #end
    #return PowerGridSolution(sol, pg), evr
    return sol, evr
end