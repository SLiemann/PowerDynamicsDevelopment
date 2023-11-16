using ModelingToolkit
using DifferentialEquations
using Plots, LinearAlgebra
import DiffEqBase: initialize_dae!
# independent variable and its differential 1.0008309834581939  0.015337314442955285
@variables t
D = Differential(t)
# Differential states
diff_states_control = @variables x_id(t)=0.00016473250985093573 x_iq(t)= 4.1944765973302266e-5 x_vd(t)=0.0 x_vq(t)=0.0 x_vdroop(t)=1.0 p_f(t)=0.32946506815322535 θ(t)=0.003229761089173144 Δv_dc(t)=0.0 Δp_cf(t)=5.77922159553466e-5 i_dcrefτ(t)=0.3795228572230594
diff_states_grid = @variables id(t)=0.32919233551126365 iq(t)=0.08495313861172944 ild(t)=0.32980112926654404 ilq(t)=-0.10354143579833089 vfd(t)=0.9999947821426575 vfq(t)=0.0032297525914168615 iLloadd(t)=-0.0012663989066732345 iLloadq(t)=-0.09932032509395052
# algebraic states
alg_states = @variables vcd(t)=0.9974904947238704 vcq(t)=0.01361411340574629 i_cdlim(t)=0.33 i_cqlim(t)=0.2884955592153876 p_ref(t)=1/3 i_dcref(t)=0.3795228572230594 i_absref0(t) = 0.3399773971043972

tmp_params = @parameters cdc=0.0956 rdc=20.0 Tdc=0.05 kii=1.19 kip=0.738891 kiv=1.161022 kvp=0.52 pref0=0.32946506488000704 vdcref=1.0 icmax=1.0 kvidroop=0.5 kvpdroop=0.001 ωf=10*pi kd=pi vref=1.0 vfqref = 0.0 idcmax=1.2 kdc=100.0 Rf = 0.0005 Lf = 0.03141592653589793/(100*pi)  Cf =  1/(5.305164769729845*100*pi) Ll = 0.049751859510499465/(100*pi) Rl = 0.004975185951049947 Rload = 3.0 Lload = 10.0/(100*pi) ω0 = 100*pi
# discrete states
discrete_states = @variables q_icmax(t)=0.0 q_idcrefτ(t)=0.0

# Auxillary variables and equations
vfdθ = vfd*cos(θ) + vfq*sin(θ)
vfqθ = vfq*cos(θ) - vfd*sin(θ)
vcdθ = vcd*cos(θ) + vcq*sin(θ)
vcqθ = vcq*cos(θ) - vcd*sin(θ)
ildθ = ild*cos(θ) + ilq*sin(θ)
ilqθ = ilq*cos(θ) - ild*sin(θ)
idθ = id*cos(θ) + iq*sin(θ)
iqθ = iq*cos(θ) - id*sin(θ)

vcdm = vfdθ + x_id + (i_cdlim - idθ)*kip - iqθ*Lf*ω0
vcqm = vfqθ + x_iq + (i_cqlim - iqθ)*kip + idθ*Lf*ω0
vcdm0 = vcdm*cos(θ) - vcqm*sin(θ)
vcqm0 = vcqm*cos(θ) + vcdm*sin(θ)
vfdref = x_vdroop + kvpdroop*(vref - sqrt(vfdθ^2 +vfqθ^2))
i_cdref = ildθ + x_vd + (vfdref - vfdθ)*kvp - vfqθ*Cf*ω0
i_cqref = ilqθ + x_vq + (vfqref - vfqθ)*kvp + vfdθ*Cf*ω0
i_absref = sqrt(i_cdref^2 + i_cqref^2)
p_measθ = vfdθ*ildθ + vfqθ*ilqθ
i_dcref0 = p_ref/vdcref +  -Δv_dc*kdc + (1.0 + Δv_dc)/rdc + Δp_cf/vdcref
Δp_c = vcdθ*idθ + vcqθ*iqθ - p_measθ
ω = ω0 + kd*(pref0 - p_f)

eqs_gfm = [
    D(id) ~ (vcd  - Rf*id + Lf*ω*iq - vfd) / Lf
    D(iq) ~ (vcq  - Rf*iq - Lf*ω*id - vfq) / Lf
    D(ild) ~ (vfd - (Rload + Rl)*ild + Ll*ω*ilq + Rload*iLloadd) / Ll
    D(ilq) ~ (vfq - (Rload + Rl)*ilq - Ll*ω*ild + Rload*iLloadq) / Ll
    D(vfd) ~ (id - ild + Cf*ω*vfq) / Cf
    D(vfq) ~ (iq - ilq - Cf*ω*vfd) / Cf
    D(iLloadd) ~ (Rload*ild - Rload*iLloadd + Lload*ω*iLloadq) / Lload
    D(iLloadq) ~ (Rload*ilq - Rload*iLloadq - Lload*ω*iLloadd) / Lload

    0.0 ~ vcd - vcdm0*(1.0 + Δv_dc)/vdcref
    0.0 ~ vcq - vcqm0*(1.0 + Δv_dc)/vdcref  
    D(x_id) ~ (i_cdlim - idθ)*kii
    D(x_iq) ~ (i_cqlim - iqθ)*kii
    0.0 ~ i_cdlim  - (1.0 - q_icmax)*i_cdref - q_icmax*i_cdref*icmax/i_absref
    0.0 ~ i_cqlim  - (1.0 - q_icmax)*i_cqref - q_icmax*i_cqref*icmax/i_absref
    D(x_vd) ~ (vfdref - vfdθ)*kiv*(1.0 - q_icmax)
    D(x_vq) ~ (vfqref - vfqθ)*kiv*(1.0 - q_icmax)
    D(x_vdroop) ~ (vref - sqrt(vfdθ^2 + vfqθ^2))*kvidroop*(1.0 - q_icmax)
    D(p_f) ~ (p_measθ - p_f)*ωf
    D(θ) ~ (p_ref - p_f)*kd
    0.0 ~ p_ref - (1.0 - q_icmax)*pref0 - q_icmax*(i_cdlim*vfdθ + i_cqlim*vfqθ)
    D(Δv_dc) ~ (i_dcref -(1.0 +Δv_dc)/rdc - (vcd*id + vcq*iq)/(1.0 +Δv_dc))/cdc
    0.0 ~ i_dcref - (1.0 - q_idcrefτ)*i_dcrefτ - q_idcrefτ*(sign(i_dcrefτ)*idcmax)
    D(i_dcrefτ) ~ (i_dcref0 - i_dcrefτ)/Tdc
    D(Δp_cf) ~ (Δp_c - Δp_cf)*ωf
    D(q_icmax) ~ 0.0
    D(q_idcrefτ) ~ 0.0
    0.0 ~ i_absref0 - i_absref
]

all_states = [diff_states_control;diff_states_grid; alg_states;discrete_states]
t_events_on = [];
t_events_off = [];

function step_load_init(integ,u,p,ctx)
    integ.p[p.Rload] = integ.p[p.Rload]*5.0/(integ.p[p.Rload] + 5.0)
    initialize_dae!(integ,BrownFullBasicInit())
end

step_load_incr = [1.0,2.0,3.0,4.0] => (step_load_init,[],[Rload],nothing)

function step_load_decrease(integ,u,p,ctx)
    integ.p[p.Rload] = 3.0
    initialize_dae!(integ,BrownFullBasicInit())
end

step_load_decr =  [6.0] => (step_load_decrease,[],[Rload],nothing)

function affect_imax_on(integ,u,p,ctx)
    append!(t_events_on,integ.t)
    integ.u[u.q_icmax] = 1.0
    initialize_dae!(integ,BrownFullBasicInit())
end  

imax_on = (i_absref-icmax-q_icmax>= 0.0) => (affect_imax_on,[q_icmax],[],nothing)

function affect_imax_off(integ,u,p,ctx)
    append!(t_events_off,integ.t)
    integ.u[u.q_icmax] = 0.0
    initialize_dae!(integ,BrownFullBasicInit())
end  
imax_off = ((i_absref-icmax)*q_icmax < 0.0) => (affect_imax_off,[q_icmax],[],[])

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
all_disc_events = [step_load_incr,step_load_decr,imax_on,imax_off,idcmax_on,idcmax_off]

ordered_states = [id,iq,ild,ilq,vfd,vfq,iLloadd,iLloadq,vcd,vcq,x_id,x_iq,i_cdlim,i_cqlim,x_vd,x_vq,x_vdroop,p_f,θ,p_ref,Δv_dc,i_dcref,i_dcrefτ,Δp_cf,q_icmax,q_idcrefτ,i_absref0]

@named odesys =ODESystem(eqs_gfm,t,ordered_states,tmp_params; discrete_events = all_disc_events) #; discrete_events = all_disc_events
tmp = ModelingToolkit.get_defaults(odesys)
ind_states = indexin(states(odesys),collect(keys(tmp)))
u0 = states(odesys) .=> collect(values(tmp))[ind_states] 
ind_params = setdiff(indexin(parameters(odesys),collect(keys(tmp))),[nothing])
params = parameters(odesys) .=> collect(values(tmp))[ind_params] 
prob  = ODEProblem(odesys, u0,(0.0,10.0),params)
sol = solve(prob,Rodas4(),dtmax=1e-3,; tstops = [1.0,2.0,3.0], initializealg = BrownFullBasicInit(),alg_hints=:stiff,abstol=1e-8,reltol=1e-8);
Plots.plot(sol,idxs=[vcd])
t_events_off
t_events_on
plot(sol,idxs=[pref])
plot(sol,idxs=[id])

using NonlinearSolve
u0f = Float64.(collect(values(tmp))[ind_states])
pf = Float64.(collect(values(tmp))[ind_params])
nprob = SteadyStateProblem(prob.f, u0f,pf)
nsol = solve(nprob,DynamicSS(Rodas4()))

pf[4] = 0.9
u0f[19] = 0.0
nprob = SteadyStateProblem(prob.f, u0f,pf)
nsol = solve(nprob,DynamicSS(Rodas4()))
states(odesys) .=> u0f .- nsol.u

pf[4] = 0.9
u0f[19] = 1.0
nprob = SteadyStateProblem(prob.f, u0f,pf)
nsol = solve(nprob,DynamicSS(Rodas4()),reltol=1e-10,abstol=1e-10)
states(odesys) .=> nsol.u

u0_steady = states(odesys) .=> nsol.u
params_steady = parameters(odesys) .=> pf
prob_steady = ODEProblem(odesys, u0_steady,(0.0,1.0),params_steady)
sol_steady = solve(prob_steady,Rodas4())
plot(sol_steady,idxs=[p_f])
plot(sol_steady,idxs=[θ])


df = ModelingToolkit.defaults(odesys);
substitute((1.0 - q_icmax)*i_cdref - q_icmax*i_cdref*icmax/i_absref,df)
substitute((1.0 - q_icmax)*i_cqref - q_icmax*i_cqref*icmax/i_absref,df)
substitute((1.0 - q_idcrefτ)*i_dcrefτ - q_idcrefτ*(sign(i_dcrefτ)*idcmax),df)

#######Eigenvalue calculatoin####################
ew_noimax = CalcEigenValues(odesys,output = false,plot=true)

ic_imax = states(odesys) .=> sol[end]
@named new_odesys = ODESystem(eqs_gfm,t,ordered_states,tmp_params; defaults = ic_imax)
ew_imax = CalcEigenValues(new_odesys,plot=true)

####################################################



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


p0 = -1/3;
q0 = -0.1;

#Line parameters
RX = 0.1; 
Zline = 0.05; 
Xl = Zline/sqrt(1+RX^2)
Rl = RX*Xl

buses=OrderedDict(
    "bus_gfm" => droop(Sbase = 100e6,Srated = 100e6, p0set = 0.5, u0set = 1.00,Kp_droop = pi/100,Kp_uset = 0.001, Ki_uset = 0.5,
                                Kdc = 100.0, gdc = Gdc, cdc = Cdc, xlf = Xlf, rf = R_f, xcf =  Xcf, Tdc = 0.05, Kp_u = 0.52,
                                Ki_u = 1.161022, Kp_i = 0.738891, Ki_i = 1.19, imax_csa = 1.0, imax_dc = 1.2, p_red = 1, LVRT_on = 0.0,
                                p_ind = collect(1:20)),
    "bus_load" => GeneralVoltageDependentLoadParam(P=p0, Q=q0, U=1.0, Ap=1.0, Bp=0.0, Aq=1.0, Bq=0.0,p_ind=[21,22]))

# HIER EXTREM KURZE LEITUNG UM DIE LAST DIREKT AM CONVERTER ZU SIMULIEREN
branches=OrderedDict(
    "line"=> PiModelLine(from= "bus_gfm", to = "bus_load",y=1.0/(Rl+1im*Xl), y_shunt_km=0.0, y_shunt_mk=0.0))
pg=  PowerGrid(buses, branches)

U1,δ1,ic0,cu = PowerFlowClassic(pg,iwamoto = true,max_tol = 1e-7,ind_sl=1)
U1
δ1
pg, ic = InitializeInternalDynamics(pg,ic0);

rhs(pg).syms .=> ic


function loadstep(integrator)
    integrator.p[21] -= 0.2
    initialize_dae!(integrator,BrownFullBasicInit())
    auto_dt_reset!(integrator)
end

cb = PresetTimeCallback([1.0,2.0,3.0], loadstep)

function loadstep_decrease(integrator)
integrator.p[21] += 0.4
initialize_dae!(integrator,BrownFullBasicInit())
auto_dt_reset!(integrator)
end

cb2 = PresetTimeCallback([5.0], loadstep_decrease)

params = getallParameters(pg.nodes["bus_gfm"])[3:22]
params = vcat(params,[p0,q0])
problem= ODEProblem{true}(rhs(pg),ic,(0.0,6.0),params)

sol1 = solve(problem, Rodas4(), callback = CallbackSet(cb,cb2), tstops=[1.0,5.0], dtmax = 1e-3,force_dtmin=false,maxiters=1e6, initializealg = BrownFullBasicInit(),alg_hints=:stiff)
pgsol = PowerGridSolution(sol1,pg);
myplot(pgsol,"bus_gfm",[:Pf])
myplot(pgsol,"bus_gfm",:θ)
myplot(pgsol,"bus_gfm",[:P0,:Pf])
myplot(pgsol,"bus_gfm",[:i_abs])

myplot(pgsol,"bus_load",[:i])

plot(pgsol,"bus_gfm",:Pf)
pgsol(1,2,:i)

symbolsof.(values(pg.nodes))