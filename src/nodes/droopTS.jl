@DynamicNode droopTS(Sbase,Srated,p0set,u0set,Kp_droop,Kp_uset,Ki_uset,Kdc,gdc,cdc,xlf,rf,xcf,Tdc,Kp_u,Ki_u,Kp_i,Ki_i,imax_csa,imax_dc,LVRT_on,p_ind) begin
    MassMatrix(m_int =[true,true,true,true,true,true,true,true,true,true,false,false,true,false,false,true,true])
end begin
    @assert Sbase > 0 "Base apparent power of the grid in VA, should be >0"
    @assert Srated > 0 "Rated apperent power of the machine in VA, should be >0"
    @assert p0set >= 0 "Set point for active power in p.u, should be >0"
    @assert u0set > 0 "Set point for voltage in p.u, should be >0"
    @assert Kp_droop >= 0 "Droop constant for active power in p.u, should be >=0"
    @assert Kp_uset >= 0 "P-Gain for voltage control, should be >=0"
    @assert Ki_uset >= 0 "I-Gain for voltage control.u, should be >=0"
    @assert Kdc >= 0 "Droop constant for current control in p.u, should be >=0"
    @assert gdc >=0 "Conductance of DC-circuit in global p.u., should be >=0"
    @assert cdc >0 "Capacitance of DC-circuit in global p.u., should be >=0"
    @assert xlf >= 0 "filter inductive reactance in p.u., should be >=0"
    @assert rf >= 0 "filter resistance in p.u., should be >=0"
    @assert xcf >= 0 "filter capacitive reactance in p.u., should be >=0"
    @assert Tdc >= 0 "DC control time constant in s, should be > 0"
    @assert Kp_u >= 0 "Proportional gain for voltage control loop, should be >0"
    @assert Ki_u >= 0 "Integral gain for voltage control loop, should be >0"
    @assert Kp_i >= 0 "Proportional gain for current control loop, should be >0"
    @assert Ki_i >= 0 "Integral gain for current control loop, should be >0"
    @assert imax_csa >= 0 "max. current for current saturation algorithm (CSA) in p.u., should be >=0"
    @assert imax_dc >= 0 "max. current of dc source in p.u., should be >=0"

end [[θ,dθ],[udc,dudc],[idc0,didc0],[x_uabs,dx_uabs],[e_ud,de_ud],[e_uq,de_uq],[e_id,de_id],[e_iq,de_iq],[Pf,dPf],[Pdelta,dPdelta],[i_abs,di_abs],[iset_abs,diset_abs],[LVRT,dLVRT],[Q0,dQ0],[P0,dP0],[q_imax,dq_imax],[q_idcmax,dq_idcmax]] begin 
    p0set = p[p_ind[1]]    
    u0set = p[p_ind[2]]
    Kp_droop = p[p_ind[3]]
    Kp_uset = p[p_ind[4]]
    Ki_uset = p[p_ind[5]]
    Kdc = p[p_ind[6]] #
    gdc = p[p_ind[7]] #
    cdc = p[p_ind[8]] #
    xlf = p[p_ind[9]]
    rf = p[p_ind[10]]
    xcf = p[p_ind[11]]
    Tdc = p[p_ind[12]] #
    Kp_u = p[p_ind[13]]
    Ki_u = p[p_ind[14]]
    Kp_i = p[p_ind[15]]
    Ki_i = p[p_ind[16]]
    imax_csa = p[p_ind[17]]
    imax_dc = p[p_ind[18]]
    LVRT_on = p[p_ind[19]]

    #after filter
    umeas = u*(cos(-θ)+1im*sin(-θ))
    udmeas = real(umeas)
    uqmeas = imag(umeas)
    imeas = i*(cos(-θ)+1im*sin(-θ))/(Srated/Sbase) #1im*
    idmeas = real(imeas)
    iqmeas = imag(imeas)
    pmeas = real(umeas * conj(imeas))
    qmeas = imag(umeas * conj(imeas))

    #before filter
    idq =  imeas + umeas / (-1im * xcf) #/ (Srated/Sbase)
    id  = real(idq)
    iq  = imag(idq)

    E = umeas + (rf + 1im*xlf) * idq
    p_before_filter = real(conj(idq) * E)
    q_before_filter = imag(conj(idq) * E)
    ix  = p_before_filter   #AC/DC coupling

    #filtered power
    dPf = 10.0*pi*(pmeas - Pf)

    #Voltage control
    Δuabs = u0set - abs(u)
    dx_uabs = Ki_uset * Δuabs * (1.0 - q_imax)
    Uset = x_uabs + Kp_uset * Δuabs

    #Building voltage reference
    udset = Uset
    uqset = 0.0

    idset = idmeas - uqmeas / xcf + Kp_u * (udset - udmeas) + e_ud
    iqset = iqmeas + udmeas / xcf + Kp_u * (uqset - uqmeas) + e_uq

    #Current saturation algorithm
    diset_abs = iset_abs - hypot(idset,iqset)
    ϕ1 = atan(iqset,idset)
    idset_csa = idset * (1.0 - q_imax) + imax_csa *cos(ϕ1) * q_imax 
    iqset_csa = iqset * (1.0 - q_imax) + imax_csa *sin(ϕ1) * q_imax 

    de_ud = (udset - udmeas) * Ki_u * (1.0 - q_imax)
    de_uq = (uqset - uqmeas) * Ki_u * (1.0 - q_imax)

    #Current control
    de_id = (idset_csa - id) * Ki_i
    de_iq = (iqset_csa - iq) * Ki_i

    umd = udmeas - iq * xlf + Kp_i * (idset_csa - id) + e_id
    umq = uqmeas + id * xlf + Kp_i * (iqset_csa - iq) + e_iq

    #Coupling with DC voltage
    um_abs = (udc + 1.0) * hypot(umd,umq)
    ϕm = atan(umq,umd)
    umds = um_abs * cos(ϕm)
    umqs = um_abs * sin(ϕm)

    #Back transformation to global reference systems
    um = umds + 1im * umqs
    um = um*(cos(θ)+1im*sin(θ))
    umd0 = real(um)
    umq0 = imag(um)

    idq = id + 1im*iq
    idq = idq*(cos(θ)+1im*sin(θ)) #-1im*
    id = real(idq) #* (Srated/Sbase)
    iq = imag(idq) #* (Srated/Sbase)
    I_abs = hypot(id,iq)
    #Voltage equations for filter
    u0d =  umd0 - rf * id + xlf * iq
    u0q =  umq0 - rf * iq - xlf * id
    u0 = u0d + 1im * u0q

    du = u - u0 #algebraic constraint
    di_abs = i_abs - I_abs
   
    #Pref reduction & droop control 
    plim = idset_csa * udmeas + iqset_csa * uqmeas 
    pmax = IfElse.ifelse(plim > 1.0, 1.0 , IfElse.ifelse(plim < -1.0 , -1.0 , plim))
    dP = (p0set -pmax) * q_imax
    w = (p0set - dP - Pf) * Kp_droop
    dθ = w

    #DC current control & DC circuit
    dPdelta = 10.0*pi*(p_before_filter - pmeas - Pdelta)
    idc = -Kdc * udc + p0set - dP + (1.0+udc)*gdc + Pdelta
    didc0 = (idc - idc0) / Tdc
    idc0_lim = idc0 * (1.0 - q_idcmax) + sign(imax_dc) * q_idcmax
    dudc = (idc0_lim - gdc * (1.0+udc) - ix) / cdc

    dP0 = P0 - pmeas
    dQ0 = Q0 - qmeas
    dLVRT = LVRT_on
    dq_imax = 0.0
    dq_idcmax = 0.0
end

export droopTS
