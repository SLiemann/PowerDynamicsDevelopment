using PowerDynamics, ModelingToolkit, DifferentialEquations
using OrderedCollections: OrderedDict
using Distributed
import DiffEqBase: initialize_dae!
@everywhere using IfElse
include("C:/Users/liemann/github/PowerDynamicsDevelopment/src/sensitivity_analyses/Local_Sensitivity.jl")
include("C:/Users/liemann/github/PowerDynamicsDevelopment/src/include_custom_nodes_lines_utilities.jl")

begin
    RX = 0.1; 
    Zline = 0.05; 
    Xl = Zline/sqrt(1+RX^2)
    Rl = RX*Xl

    buses=OrderedDict(
        "bus_grid" => SlackAlgebraicParam(U=1.0,p_ind=[1]),
        "bus_load" => nPFC(Cd=0.036, Pdc = 1.0, p_offset = 0.0, q_offset=0.0, p_ind=collect(2:3)))

    # HIER EXTREM KURZE LEITUNG UM DIE LAST DIREKT AM CONVERTER ZU SIMULIEREN
    branches=OrderedDict(
        "line"=> PiModelLine(from= "bus_grid", to = "bus_load",y=1.0/(Rl+1im*Xl), y_shunt_km=0.0, y_shunt_mk=0.0))
    pg=  PowerGrid(buses, branches)

    U1,δ1,ic0,cu = PowerFlowClassic(pg,iwamoto = true,max_tol = 1e-7,ind_sl=1)
    pg, ic = InitializeInternalDynamics(pg,ic0);

    index_U_load = PowerDynamics.variable_index(pg.nodes,"bus_load",:u_r)
    index_vofft2_load = PowerDynamics.variable_index(pg.nodes,"bus_load",:vofft2)
    index_tsum_load = PowerDynamics.variable_index(pg.nodes,"bus_load",:tsum)
    index_ton_load = PowerDynamics.variable_index(pg.nodes,"bus_load",:ton)
    index_toff_load = PowerDynamics.variable_index(pg.nodes,"bus_load",:toff)
    index_Vabstoff_load = PowerDynamics.variable_index(pg.nodes,"bus_load",:Vabstoff)
    index_qon_load = PowerDynamics.variable_index(pg.nodes,"bus_load",:q_on)
    p_pel_ind = pg.nodes["bus_load"].p_ind

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
        tsum = deepcopy(integrator.u[index_tsum_load])  
        ton = deepcopy(integrator.u[index_ton_load])
        toff = deepcopy(integrator.u[index_toff_load]) 
        tsum = deepcopy(integrator.u[index_tsum_load])
        Vabstoff = deepcopy(integrator.u[index_Vabstoff_load])

        if tsum >= 0.01         
            integrator.u[index_tsum_load]  = 0.0
            if ton >= 0.0
                voff = Vabstoff*sin(100*pi*toff)
                voffT2 = CalfnPFCVoffT2(voff,Pdc,Cd,(T/2-toff))
            else
                voff = deepcopy(integrator.u[index_vofft2_load])
                voffT2 = CalfnPFCVoffT2(voffT2,Pdc,Cd,(T/2-0.0))
            end
            toff = CalcnPFCtoff(Vabs,Pdc,Cd)
            Vabstoff = Vabs*sin(100*pi*toff)
            ton = CalfnPFCton(Vabs,Pdc,Cd,voffT2) 
            integrator.u[index_vofft2_load]  = voffT2
        elseif  tsum < toff && tsum < ton #Location B
            #display("B ton pos")
            toff = CalcnPFCtoff(Vabs,Pdc,Cd)
            Vabstoff = Vabs*sin(100*pi*toff)
            ton = CalfnPFCton(Vabs,Pdc,Cd,voffT2)
        elseif  tsum < toff && tsum >= ton
            # display("C ton pos")
            toff = CalcnPFCtoff(Vabs,Pdc,Cd)
            Vabstoff = Vabs*sin(100*pi*toff)
        end

        if ton >= 0.0
            integrator.u[index_qon_load] = 1.0
        else
            integrator.u[index_qon_load] = 0.0           
        end 

        integrator.u[index_Vabstoff_load]  = Vabstoff
        integrator.u[index_toff_load]  = toff
        integrator.u[index_ton_load]  = ton
        integrator.u[index_tsum_load]  = integrator.u[index_tsum_load] + integrator.t - integrator.tprev 
        initialize_dae!(integrator,BrownFullBasicInit())
    end

    cb_tsum = DiscreteCallback(f_tsum, affect_tsum, save_positions=(false,false))
    cb1 = PresetTimeCallback([0.1], error, save_positions=(false,false))
    cb2 = PresetTimeCallback([0.2], regular, save_positions=(false,false))

    tspan = (0.0,0.3)
    params = [1.0, 0.036,1.0]
    problem= ODEProblem{true}(rhs(pg),ic,tspan,params)
    tstops_sim =collect(tspan[1]:0.01:tspan[2]);
    sort!(tstops_sim);
    nothing
end

sol = solve(problem, Rodas4(autodiff=true),callback = CallbackSet(cb1,cb2,cb_tsum), tstops=tstops_sim, dtmax = 1e-4,force_dtmin=true,maxiters=1e6, initializealg = BrownFullBasicInit(),alg_hints=:stiff,abstol=1e-8,reltol=1e-8);
pgsol0 = PowerGridSolution(sol,pg);

plotallvoltages(pgsol0);

myplot(pgsol0,"bus_load",:q1)
myplot(pgsol0,"bus_load",:p1);
myplot(pgsol0,"bus_load",:tsum);
myplot(pgsol0,"bus_load",:ton);
myplot(pgsol0,"bus_load",:toff);
myplot(pgsol0,"bus_load",:vofft2);
myplot(pgsol0,"bus_load",:Vabstoff);
myplot(pgsol0,"bus_load",:q_on);


### PEL Callbacks END ###
function error(integrator) 
    integrator.p[1] = 0.5
    initialize_dae!(integrator,BrownFullBasicInit())
    ## Init PEL Model
    
    # init = false
    # while init == false
    #     Vabs = hypot(integrator[index_U_load],integrator[index_U_load+1])*sqrt(2)
    #     Vabstoff = integrator[index_Vabstoff_load]
    #     Cd = deepcopy(integrator.p[p_pel_ind[1]])
    #     Pdc = deepcopy(integrator.p[p_pel_ind[2]])
    #     toff = integrator[index_toff_load]
    #     T = 0.02
    #     ω0 = 100*pi

    #     # Hit with voltage?
    #     voff_new = Vabstoff*sin(ω0*toff)
    #     VoffT2_new = CalfnPFCVoffT2(voff_new,Pdc,Cd,(T/2-toff))
    #     ton_new = CalfnPFCton(Vabs,Pdc,Cd,VoffT2_new)  

    #     if ton_new >= 0.0
    #         integrator.u[index_qon_load] = 1.0
    #     else
    #         integrator.u[index_qon_load] = 0.0           
    #     end 

    #     toff_new = CalcnPFCtoff(Vabs,Pdc,Cd)
    #     Vabstoff_new = Vabs*sin(100*pi*toff_new)
    #     # display(toff_new)
    #     # display(Vabstoff_new)
    #     # display(VoffT2_new)
    #     # display(ton_new)
    #     integrator.u[index_vofft2_load]  = VoffT2_new
    #     integrator.u[index_Vabstoff_load]  = Vabstoff_new
    #     integrator.u[index_toff_load]  = toff_new
    #     integrator.u[index_ton_load]  = ton_new
    #     integrator.u[index_tsum_load]  = 0.0
    #     initialize_dae!(integrator,BrownFullBasicInit())
    #     Vabs_after = hypot(integrator[index_U_load],integrator[index_U_load+1])*sqrt(2)
    #     display(abs(Vabs_after-Vabs))
    #     if abs(Vabs_after - Vabs) < 0.01
    #         init = true
    #     end
    # end
end
function regular(integrator) 
    integrator.p[1] = 1.0

    initialize_dae!(integrator,BrownFullBasicInit())

    # init = false

    # while init == false
    #     Vabs = hypot(integrator[index_U_load],integrator[index_U_load+1])*sqrt(2)
    #     Vabstoff = integrator[index_Vabstoff_load]
    #     Cd = deepcopy(integrator.p[p_pel_ind[1]])
    #     Pdc = deepcopy(integrator.p[p_pel_ind[2]])
    #     toff = integrator[index_toff_load]
    #     T = 0.02
    #     ω0 = 100*pi

    #     # Hit with voltage?
    #     voff_new = Vabstoff*sin(ω0*toff)
    #     VoffT2_new = CalfnPFCVoffT2(voff_new,Pdc,Cd,(T/2-toff))
    #     ton_new = CalfnPFCton(Vabs,Pdc,Cd,VoffT2_new)  

    #     if ton_new >= 0.0
    #         integrator.u[index_qon_load] = 1.0
    #     else
    #         integrator.u[index_qon_load] = 0.0           
    #     end 

    #     toff_new = CalcnPFCtoff(Vabs,Pdc,Cd)
    #     Vabstoff_new = Vabs*sin(100*pi*toff_new)
    #     # display(toff_new)
    #     # display(Vabstoff_new)
    #     # display(VoffT2_new)
    #     # display(ton_new)
    #     integrator.u[index_vofft2_load]  = VoffT2_new
    #     integrator.u[index_Vabstoff_load]  = Vabstoff_new
    #     integrator.u[index_toff_load]  = toff_new
    #     integrator.u[index_ton_load]  = ton_new
    #     integrator.u[index_tsum_load]  = 0.0
    #     initialize_dae!(integrator,BrownFullBasicInit())
    #     Vabs_after = hypot(integrator[index_U_load],integrator[index_U_load+1])*sqrt(2)
    #     display(abs(Vabs_after-Vabs))
    #     if abs(Vabs_after - Vabs) < 0.01
    #         init = true
    #     end
    # end
end