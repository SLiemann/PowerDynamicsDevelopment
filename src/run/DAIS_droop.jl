using ModelingToolkit


# independent variable and its differential
@variables t
D = Differential(t)
# Differential states
@variables x_id(t) x_iq(t) x_vd(t) x_vq(t) x_vdroop(t) x_id(t) p_f(t) Δθ(t) Δv_dc(t) Δp_cf(t) i_dcrefτ(t)
# algebraic states
@variables v_d(t) v_q(t) v_cd(t) v_cq(t) i_cdlim(t) i_cqlim(t) c_awud(t) c_awuq(t) c_awuvdroop(t) p_ref(t) i_dcref(t)
# parameters
@parameters rf xlf xcf rload xload cdc rdc Tdc kii kip kiv kvp pref0 vdcref=1.0 icmax=1.0 kvidroop kvpdroop ωf=10*pi kd vref=1.0 v_qref = 0.0 idcmax=1.25 kdc
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

display.(equations(ode))