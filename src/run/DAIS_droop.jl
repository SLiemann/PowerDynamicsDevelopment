using ModelingToolkit
using DifferentialEquations
using Plots, LinearAlgebra
import DiffEqBase: initialize_dae!
# independent variable and its differential 1.0008309834581939  0.015337314442955285
@variables t
D = Differential(t)
# Differential states
diff_states = @variables x_id(t)=0.0002499998300626934 x_iq(t)= 4.424770956641677e-5 x_vd(t)=0.0 x_vq(t)=0.0 x_vdroop(t)=1.0 p_f(t)=0.5 Δθ(t)=0.0 Δv_dc(t)=0.0 Δp_cf(t)=0.00012891554966654173 i_dcrefτ(t)=0.5501289155496666
# algebraic states
alg_states = @variables v_d(t)=1.0 v_q(t)=0.0 v_cd(t)=1.00 v_cq(t)=0.0  i_cdlim(t)=0.5 i_cqlim(t)=0.08849555921538757 p_ref(t)=0.5 i_dcref(t)=0.5501289155496666

tmp_params = @parameters rf=0.0005 xlf=0.03141592653589793 xcf=5.305164769729845 rload=2.0 xload=10.0 cdc=0.0956 rdc=20.0 Tdc=0.05 kii=1.19 kip=0.738891 kiv=1.161022 kvp=0.52 pref0=0.5 vdcref=1.0 icmax=1.0 kvidroop=0.5 kvpdroop=0.001 ωf=10*pi kd=pi vref=1.0 v_qref = 0.0 idcmax=1.2 kdc=100.0
# discrete states
discrete_states = @variables q_icmax(t)=0.0 q_idcrefτ(t)=0.0


# Auxillary variables and equations
v_dθ = v_d*cos(Δθ) + v_q*sin(Δθ)
v_qθ = v_q*cos(Δθ) - v_d*sin(Δθ)
v_cdθ = v_cd*cos(Δθ) + v_cq*sin(Δθ)
v_cqθ = v_cq*cos(Δθ) - v_cd*sin(Δθ)
i_dθ = v_dθ/rload + v_qθ/xload
i_qθ = v_qθ/rload - v_dθ/xload
i_cfdθ = -v_qθ/xcf
i_cfqθ = v_dθ/xcf
i_cdθ = i_dθ + i_cfdθ
i_cqθ = i_qθ + i_cfqθ
vcdm = v_dθ + x_id + (i_cdlim - i_cdθ)*kip - i_cqθ*xlf
vcqm = v_qθ + x_iq + (i_cqlim - i_cqθ)*kip + i_cdθ*xlf
vcdm0 = vcdm*cos(Δθ) - vcqm*sin(Δθ)
vcqm0 = vcqm*cos(Δθ) + vcdm*sin(Δθ)
v_dref = x_vdroop + kvpdroop*(vref - sqrt(v_dθ^2 +v_qθ^2))
i_cdref = i_dθ + x_vd + (v_dref - v_dθ)*kvp - v_qθ/xcf
i_cqref = i_qθ + x_vq + (v_qref - v_qθ)*kvp + v_dθ/xcf
i_absref = sqrt(i_cdref^2 + i_cqref^2)
p_measθ = v_dθ*i_dθ + v_qθ*i_qθ
i_cd = v_d/rload + v_q/xload - v_q/xcf
i_cq = v_q/rload - v_d/xload + v_d/xcf
i_dcref0 = p_ref/vdcref +  -Δv_dc*kdc + (1.0 + Δv_dc)/rdc + Δp_cf/vdcref
Δp_c = v_cdθ*i_cdθ + v_cqθ*i_cqθ - p_measθ

#i_d = i_dθ*cos(Δθ) - i_qθ*sin(Δθ)
#i_q = i_qθ*cos(Δθ) + i_dθ*sin(Δθ)
#i_abs = hypot(i_d,i_q)

eqs_gfm = [ 
    0.0 ~ v_d - v_cd +  rf*(v_d/rload + v_q/xload - v_q/xcf) - xlf*(v_q/rload - v_d/xload + v_d/xcf)
    0.0 ~ v_q - v_cq + xlf*(v_d/rload + v_q/xload - v_q/xcf) +  rf*(v_q/rload - v_d/xload + v_d/xcf)
    0.0 ~ v_cd - vcdm0*(1.0 + Δv_dc)/vdcref
    0.0 ~ v_cq - vcqm0*(1.0 + Δv_dc)/vdcref  
    D(x_id) ~ (i_cdlim - i_cdθ)*kii
    D(x_iq) ~ (i_cqlim - i_cqθ)*kii
    0.0 ~ i_cdlim  - (1.0 - q_icmax)*i_cdref - q_icmax*i_cdref*icmax/i_absref
    0.0 ~ i_cqlim  - (1.0 - q_icmax)*i_cqref - q_icmax*i_cqref*icmax/i_absref
    D(x_vd) ~ (v_dref - v_dθ)*kiv*(1.0 - q_icmax)
    D(x_vq) ~ (v_qref - v_qθ)*kiv*(1.0 - q_icmax)
    D(x_vdroop) ~ (vref - sqrt(v_dθ^2 + v_qθ^2))*kvidroop*(1.0 - q_icmax)
    D(p_f) ~ (p_measθ - p_f)*ωf
    D(Δθ) ~ (p_ref - p_f)*kd
    0.0 ~ p_ref - (1.0 - q_icmax)*pref0 - q_icmax*(i_cdlim*v_dθ + i_cqlim*v_qθ)
    D(Δv_dc) ~ (i_dcref -(1.0 +Δv_dc)/rdc - (v_cd*i_cd + v_cq*i_cq)/(1.0 +Δv_dc))/cdc
    0.0 ~ i_dcref - (1.0 - q_idcrefτ)*i_dcrefτ - q_idcrefτ*(sign(i_dcrefτ)*idcmax)
    D(i_dcrefτ) ~ (i_dcref0 - i_dcrefτ)/Tdc
    D(Δp_cf) ~ (Δp_c - Δp_cf)*ωf
    D(q_icmax) ~ 0.0
    D(q_idcrefτ) ~ 0.0
]

all_states = [diff_states; alg_states;discrete_states]
step_load = [1.0,2.0,3.0] => [rload ~ rload*5.0/(rload + 5.0)]

function affect_imax_on(integ,u,p,ctx)
    integ.u[u.q_icmax] = 1.0
    initialize_dae!(integ,BrownFullBasicInit())
end  
imax_on = (i_absref >= icmax) => (affect_imax_on,[q_icmax],[],nothing)

function affect_imax_off(integ,u,p,ctx)
    integ.u[u.q_icmax] = 0.0
    initialize_dae!(integ,BrownFullBasicInit())
end  
imax_off = (i_absref < icmax) => (affect_imax_off,[q_icmax],[],[])

function affect_idcmax_on(integ,u,p,ctx)
    integ.u[u.q_idcrefτ] = 1.0
    initialize_dae!(integ,BrownFullBasicInit())
end  
idcmax_on = (i_dcrefτ >= idcmax) => (affect_idcmax_on,[q_idcrefτ],[],[])

function affect_idcmax_off(integ,u,p,ctx)
    integ.u[u.q_idcrefτ] = 0.0
    initialize_dae!(integ,BrownFullBasicInit())
end  
idcmax_off = (i_dcrefτ < idcmax) => (affect_idcmax_off,[q_idcrefτ],[],[])

all_disc_events = [step_load,imax_on,imax_off,idcmax_on,idcmax_off]

ordered_states = [v_d,v_q,v_cd,v_cq,x_id,x_iq,i_cdlim,i_cqlim,x_vd,x_vq,x_vdroop,p_f,Δθ,p_ref,Δv_dc,i_dcref,i_dcrefτ,Δp_cf,q_icmax,q_idcrefτ]
loose_states = [diff_states; alg_states;discrete_states]

@named odesys = ODESystem(eqs_gfm,t,ordered_states,tmp_params; discrete_events = all_disc_events)
tmp = ModelingToolkit.get_defaults(odesys)

ind_states = indexin(states(odesys),collect(keys(tmp)))
u0 = states(odesys) .=> collect(values(tmp))[ind_states] 
ind_params = setdiff(indexin(parameters(odesys),collect(keys(tmp))),[nothing])
params = parameters(odesys) .=> collect(values(tmp))[ind_params] 

prob  = ODEProblem(odesys, u0,(0.0,5.0),params)
sol = solve(prob,Rodas4(),dtmax=1e-3,; tstops = [1.0,2.0,3.0], initializealg = BrownFullBasicInit(),alg_hints=:stiff,abstol=1e-8,reltol=1e-8);
plot(sol,idxs=[p_f])

plot!(pgsol,"bus_gfm",:Pf,ylims=[0.88,0.92])

plot(sol,idxs=[Δθ])
plot!(pgsol,"bus_gfm",:θ)

plot(sol,idxs=[q_cmax])

plot(sol,idxs=[Δv_dc])
plot!(pgsol,"bus_gfm",:udc)

plot(sol,idxs=[v_q])
plot!(pgsol,"bus_gfm",:u_i)

plot(sol,idxs=[Δp_cf])
plot!(pgsol,"bus_gfm",:Pdelta)

df = ModelingToolkit.defaults(odesys);
substitute((1.0 - q_icmax)*i_cdref - q_icmax*i_cdref*icmax/i_absref,df)
substitute((1.0 - q_icmax)*i_cqref - q_icmax*i_cqref*icmax/i_absref,df)
substitute((1.0 - q_idcrefτ)*i_dcrefτ - q_idcrefτ*(sign(i_dcrefτ)*idcmax),df)

#######Eigenvalue calculatoin####################
ew_noimax = CalcEigenValues(odesys,output = false,plot=true)

ic_imax = states(odesys) .=> sol[end]
@named new_odesys = ODESystem(eqs_gfm,t,ordered_states,tmp_params; defaults = ic_imax)
ew_imax = CalcEigenValues(new_odesys,plot=true)

f1 = my_rhs.(eqs_gfm);
states(odesys) .=> substitute(Num.(f1),new_odesys.defaults)

#Taken from states(simplified_sys)
char_states = ["x_id(t)",
"x_iq(t)",
"x_vd(t)",
"x_vq(t)",
"x_vdroop(t)",
"p_f(t)",
"Δθ(t)",
"Δv_dc(t)",
"i_dcrefτ(t)",
"Δp_cf(t)",
"q_icmax(t)",
"q_idcrefτ(t)",
]

using MATLAB
write_matfile("eigvals_GFM_DAIS.mat"; noimax = ew_noimax,imax=ew_imax,states=char_states)

# -> nach Khalil(S. 42-43): Eigenwerte gleich null bedeuten, dass es eine Vielzahl von Equilibrium points gibt (verschiedene Spannungen)



using PowerDynamics
using OrderedCollections: OrderedDict
using Distributed
import DiffEqBase: initialize_dae!
@everywhere using IfElse
include("C:/Users/liemann/github/PowerDynamicsDevelopment/src/include_custom_nodes_lines_utilities.jl")

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


p0 = -0.5;
q0 = -0.1;

buses=OrderedDict(
    "bus_gfm" => droop(Sbase = 100e6,Srated = 100e6, p0set = 0.5, u0set = 1.00,Kp_droop = pi,Kp_uset = 0.001, Ki_uset = 0.5,
                                Kdc = 100.0, gdc = Gdc, cdc = Cdc, xlf = Xlf, rf = R_f, xcf =  Xcf, Tdc = 0.05, Kp_u = 0.52,
                                Ki_u = 1.161022, Kp_i = 0.738891, Ki_i = 1.19, imax_csa = 1.0, imax_dc = 1.2, p_red = 1, LVRT_on = 0.0,
                                p_ind = collect(1:20)),
    "bus_load" => GeneralVoltageDependentLoadParam(P=p0, Q=q0, U=1.0, Ap=1.0, Bp=0.0, Aq=1.0, Bq=0.0,p_ind=[21,22]))

# HIER EXTREM KURZE LEITUNG UM DIE LAST DIREKT AM CONVERTER ZU SIMULIEREN
branches=OrderedDict(
    "line"=> PiModelLine(from= "bus_gfm", to = "bus_load",y=1.0/(0.000001+1im*0.000001), y_shunt_km=0.0, y_shunt_mk=0.0))
pg=  PowerGrid(buses, branches)

#nodes_postfault = deepcopy(pg.nodes)
#branches_postfault = deepcopy(pg.lines)
#nodes_postfault["bus_load"] = VoltageDependentLoad(P=-1.5, Q=-0.1, U=1.0, A=1.0, B=0.0,Y_n = complex(0.0))
#pg2 =  PowerGrid(nodes_postfault,branches_postfault)

U1,δ1,ic0,cu = PowerFlowClassic(pg,iwamoto = true,max_tol = 1e-7,ind_sl=1)
pg, ic = InitializeInternalDynamics(pg,ic0);
#pg.nodes["bus_gfm"].p0set
#display.(rhs(pg).syms .=> ic)

function loadstep(integrator)
        integrator.p[21] -= 0.2
        initialize_dae!(integrator,BrownFullBasicInit())
        auto_dt_reset!(integrator)
end

cb = PresetTimeCallback([1.0,2.0,3.0], loadstep)
params = getallParameters(pg.nodes["bus_gfm"])[3:22]
params = vcat(params,[p0,q0])
problem= ODEProblem{true}(rhs(pg),ic,(0.0,5.0),params)

sol1 = solve(problem, Rodas4(), callback = cb, tstops=[1.0], dtmax = 1e-3,force_dtmin=false,maxiters=1e6, initializealg = BrownFullBasicInit(),alg_hints=:stiff,abstol=1e-8,reltol=1e-8);
pgsol = PowerGridSolution(sol1,pg);
myplot(pgsol,"bus_gfm",[:Pf])
plotallvoltages(pgsol)
myplot(pgsol,"bus_gfm",:θ)
myplot(pgsol,"bus_gfm",[:P0,:Pf])
myplot(pgsol,"bus_gfm",[:i_abs])


plot(pgsol,"bus_gfm",:Pf)


vr1 =  ExtractResult(pgsol,:u_r_1)
vi1 =  ExtractResult(pgsol,:u_i_1)
v = hypot.(vr1,vi1)
iabs =  ExtractResult(pgsol,:i_abs_1)
p0 =  ExtractResult(pgsol,:P0_1)
q0 =  ExtractResult(pgsol,:Q0_1)
tsol = pgsol.dqsol.t
theta =  ExtractResult(pgsol,:θ_1)
udc =  ExtractResult(pgsol,:udc_1)
plot(udc)
write_matfile("sol_GFM_DAIS.mat",t=tsol,vgfm = v,iabs=iabs,p=p0,q=q0,theta=theta,udc=udc)