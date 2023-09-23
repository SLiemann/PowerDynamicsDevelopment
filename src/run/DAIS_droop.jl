using ModelingToolkit

# independent variable and its differential
@variables t
D = Differential(t)
# Differential states
@variables x_id(t)=0.0002472492531056434 x_iq(t)=3.205784091951831e-5 x_vd(t)=0.0 x_vq(t)=0.0 x_vdroop(t)=1.0 p_f(t)=0.4944985062113048 Δθ(t)=0.0 Δv_dc(t)=0.0 Δp_cf(t)=0.00012431979665145088 i_dcrefτ(t)=0.5501243197966514
# algebraic states
@variables v_d(t)=0.982617027200252 v_q(t)=-0.048206051847367 v_cd(t)=1.0 v_cq(t)=0.0 i_cdlim(t) i_cqlim(t) c_awud(t)=0.0 c_awuq(t)=0.0 c_awuvdroop(t)=0.0 p_ref(t)=0.5 i_dcref(t)=0.5501243197966514
# parameters
@parameters rf=0.0005 xlf=0.03141592653589793 xcf=5.305164769729845 rload=2.0 xload=10.0 cdc=0.956 rdc=20.0 Tdc=0.05 kii=1.19 kip=0.738891 kiv=1.161022 kvp=0.52 pref0=0.5 vdcref=1.0 icmax=1.0 kvidroop=0.5 kvpdroop=0.001 ωf=10*pi kd vref=1.0 v_qref = 0.0 idcmax=1.2 kdc=100
# discrete states
@variables q_icmax(t)=0.0 q_idcrefτ(t)=0.0

# Auxillary variables and equations
v_dθ = v_d*cos(Δθ) + v_q*sin(Δθ)
v_qθ = v_q*cos(Δθ) - v_d*sin(Δθ)
v_cdθ = v_cd*cos(Δθ) + v_cq*sin(Δθ)
v_cqθ = v_cq*cos(Δθ) - v_cd*sin(Δθ)
i_cdθ = v_dθ/rload + v_qθ/xload
i_cqθ = v_qθ/rload - v_dθ/xload
i_cfdθ = v_qθ/xcf
i_cfqθ = -v_dθ/xcf
i_dθ = i_cdθ -i_cfdθ
i_qθ = i_cqθ -i_cfqθ
vcdm = v_dθ + x_id + (i_cdlim - i_cdθ)*kip - i_cqθ*xlf
vcqm = v_qθ + x_iq + (i_cqlim - i_cqθ)*kip + i_cdθ*xlf
v_dref = x_vdroop + kvpdroop*(vref - sqrt(v_dθ^2 +v_qθ^2))
i_cdref = i_dθ + x_vd + (v_dref - v_dθ)*kvp - v_qθ/xcf
i_cqref = i_qθ + x_vq + (v_qref - v_qθ)*kvp + v_dθ/xcf
i_absref = sqrt(i_cdref^2 + i_cqref^2)
p_measθ = v_dθ*i_dθ + v_qθ*i_qθ
i_cd = v_d/rload + v_q/xload
i_cq = v_q/rload - v_d/xload
i_dcref0 = p_ref/vdcref + (vdcref - 1.0 -Δv_dc)*kdc + (1.0 + Δv_dc)/rdc + Δp_cf/vdcref
Δp_c = v_cdθ*i_cdθ + v_cqθ*i_cqθ - p_measθ

eqs = [ 
    0.0 ~ v_d +  rf*(v_d/rload + v_q/xload - v_q/xcf) - xlf*(v_q/rload - v_d/xload - v_d/xcf)
    0.0 ~ v_q + xlf*(v_d/rload + v_q/xload - v_q/xcf) +  rf*(v_q/rload - v_d/xload - v_d/xcf)
    0.0 ~ v_cd - vcdm*(1.0 + Δv_dc)/vdcref
    0.0 ~ v_cq - vcqm*(1.0 + Δv_dc)/vdcref  
    D(x_id) ~ (i_cdlim - i_cdθ)*kii
    D(x_iq) ~ (i_cqlim - i_cqθ)*kii
    0.0 ~ i_cdlim  - (1.0 - q_icmax)*i_cdref - q_icmax*i_cdref*icmax/i_absref
    0.0 ~ i_cqlim  - (1.0 - q_icmax)*i_cqref - q_icmax*i_cqref*icmax/i_absref
    D(x_vd) ~ (v_dref - v_dθ)*kiv + c_awud
    D(x_vq) ~ (v_qref - v_qθ)*kiv + c_awuq
    0.0 ~ c_awud - (1.0 - q_icmax)*0.0 - q_icmax*(-(v_dref - v_dθ)*kiv)
    0.0 ~ c_awuq - (1.0 - q_icmax)*0.0 - q_icmax*(-(v_qref - v_qθ)*kiv)
    D(x_vdroop) ~ (vref - sqrt(v_dθ^2 + v_qθ^2))*kvidroop + c_awuvdroop
    0.0 ~ c_awuvdroop - (1.0 - q_icmax)*0.0 - q_icmax*(-(vref - sqrt(v_dθ^2 + v_qθ^2))*kvidroop)
    D(p_f) ~ (p_measθ - p_f)*ωf
    D(Δθ) ~ (p_ref - p_f)*kd
    0.0 ~ p_ref - (1.0 - q_icmax)*pref0 - q_icmax*(i_cdlim*v_dθ + i_cqlim*v_qθ)
    D(Δv_dc) ~ (i_dcref -(1.0 +Δv_dc)/rdc - (v_cd*i_cd + v_cq*i_cq)/(1.0 +Δv_dc))
    0.0 ~ i_dcref - (1.0 - q_idcrefτ)*i_dcrefτ - q_idcrefτ*(sign(i_dcrefτ)*idcmax)
    D(i_dcrefτ) ~ (i_dcref0 - i_dcrefτ)/Tdc
    D(Δp_cf) ~ (Δp_c - Δp_cf)*ωf
]

@named ode = ODESystem(eqs)
df = ModelingToolkit.defaults(ode);
substitute(i_dθ + x_vd + (v_dref - v_dθ)*kvp - v_qθ/xcf,ModelingToolkit.defaults(ode))
substitute(i_qθ + x_vq + (v_qref - v_qθ)*kvp + v_dθ/xcf,ModelingToolkit.defaults(ode))

substitute(i_dθ,df)


i_qθ + x_vq + (v_qref - v_qθ)*kvp + v_dθ/xcf
ModelingToolkit.defaults(ode)




using PowerDynamics
using OrderedCollections: OrderedDict
using Distributed
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


buses=OrderedDict(
    "bus_gfm" => droop(Sbase = 100e6,Srated = 100e6, p0set = 0.5, u0set = 1.00,Kp_droop = pi,Kp_uset = 0.001, Ki_uset = 0.5,
                                Kdc = 100.0, gdc = Gdc, cdc = Cdc, xlf = Xlf, rf = R_f, xcf =  Xcf, Tdc = 0.05, Kp_u = 0.52,
                                Ki_u = 1.161022, Kp_i = 0.738891, Ki_i = 1.19, imax_csa = 1.0, imax_dc = 1.2, p_red = 1, LVRT_on = 0.0,
                                p_ind = collect(6:25)),
    "bus_load" => VoltageDependentLoad(P=-0.5, Q=-0.1, U=1.0, A=0.0, B=1.0,Y_n = complex(0.0)))
branches=OrderedDict(
    "line"=> PiModelLine(from= "bus_gfm", to = "bus_load",y=1.0/(0.01+1im*0.1), y_shunt_km=0.0, y_shunt_mk=0.0))
pg=  PowerGrid(buses, branches)

U1,δ1,ic0,cu = PowerFlowClassic(pg,iwamoto = false,max_tol = 1e-4,ind_sl=1)
pg, ic = InitializeInternalDynamics(pg,ic0)