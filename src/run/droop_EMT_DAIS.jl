using ModelingToolkit
using DifferentialEquations
using Plots, LinearAlgebra
import DiffEqBase: initialize_dae!
# independent variable and its differential 1.0008309834581939  0.015337314442955285
@variables t
D = Differential(t)
# Differential states
diff_states_control = @variables x_id(t)=0.0002499998300626934 x_iq(t)= 4.424770956641677e-5 x_vd(t)=0.0 x_vq(t)=0.0 x_vdroop(t)=1.0 p_f(t)=0.5 θ(t)=0.0 Δv_dc(t)=0.0 Δp_cf(t)=0.00012891554966654173 i_dcrefτ(t)=0.5501289155496666
diff_states_grid = @variables id(t)=1.0 iq(t)=1.0 ild(t)=1.0 ilq(t)=1.0 vfd(t)=1.0 vfq(t)=1.0 iLloadd(t)=1.0 iLloadq(t)=1.0
# algebraic states
alg_states = @variables vcd(t)=1.00 vcq(t)=0.0 i_cdlim(t)=0.5 i_cqlim(t)=0.08849555921538757 p_ref(t)=0.5 i_dcref(t)=0.5501289155496666 i_absref0(t) = sqrt(0.08849555921538757^2 + 0.5^2)

tmp_params = @parameters cdc=0.0956 rdc=20.0 Tdc=0.05 kii=1.19 kip=0.738891 kiv=1.161022 kvp=0.52 pref0=0.5 vdcref=1.0 icmax=1.0 kvidroop=0.5 kvpdroop=0.001 ωf=10*pi kd=pi vref=1.0 vfqref = 0.0 idcmax=1.2 kdc=100.0 Rf = 1.0 Lf = 1.0 Cf = 1.0 Ll = 1.0 Rl = 1.0 Rload = 1.0 Lload = 1.0 ω0 = 100*pi
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

eqs_gfm = [
    D(id) ~ (vcd  - Rf*id + Lf*ω0*iq - vfd) / Lf
    D(iq) ~ (vcq  - Rf*iq - Lf*ω0*id - vfq) / Lf
    D(ild) ~ (vfd - (Rload + Rl)*ild + Ll*ω0*ilq + Rload*iLloadd) / Ll
    D(ilq) ~ (vfq - (Rload + Rl)*ilq - Ll*ω0*ild + Rload*iLloadq) / Ll
    D(vfd) ~ (id - ild + Cf*ω0*vfq) / Cf
    D(vfq) ~ (iq - ilq - Cf*ω0*vfd) / Cf
    D(iLloadd) ~ (Rload*ild - Rload*iLloadd + Lload*ω0*iLloadq) / Lload
    D(iLloadq) ~ (Rload*ilq - Rload*iLloadq - Lload*ω0*iLloadd) / Lload

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

step_load_incr = [1.0,2.0,3.0] => (step_load_init,[],[Rload],nothing)

function step_load_decrease(integ,u,p,ctx)
    integ.p[p.Rload] = 2.0
    initialize_dae!(integ,BrownFullBasicInit())
end

step_load_decr =  [5.0] => (step_load_decrease,[],[Rload],nothing)

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

ordered_states = [vfd,vfq,vcd,vcq,x_id,x_iq,i_cdlim,i_cqlim,x_vd,x_vq,x_vdroop,p_f,θ,p_ref,Δv_dc,i_dcref,i_dcrefτ,Δp_cf,q_icmax,q_idcrefτ,i_absref0]
loose_states = [diff_states; alg_states;discrete_states]

@named odesys =ODESystem(eqs_gfm,t,ordered_states,tmp_params; discrete_events = all_disc_events)
tmp = ModelingToolkit.get_defaults(odesys)
ind_states = indexin(states(odesys),collect(keys(tmp)))
u0 = states(odesys) .=> collect(values(tmp))[ind_states] 
ind_params = setdiff(indexin(parameters(odesys),collect(keys(tmp))),[nothing])
params = parameters(odesys) .=> collect(values(tmp))[ind_params] 
prob  = ODEProblem(odesys, u0,(0.0,6.0),params)
sol = solve(prob,Rodas4(),dtmax=1e-4,; tstops = [1.0,2.0,3.0], initializealg = BrownFullBasicInit(),alg_hints=:stiff,abstol=1e-8,reltol=1e-8);
plot(sol,idxs=[q_icmax])
t_events_off
t_events_on
plot(sol,idxs=[p_f])
plot(sol,idxs=[vcq])

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

