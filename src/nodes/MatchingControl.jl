#= Sebastian Liemann, ie3 TU Dortmund, based on F. Milano, Power System Modelling and Scripting, Springer Verlag, 2010
@doc doc"""
```Julia
GridFormingConverter(Sbase,Srated,p0set,u0set,Kp_droop,Kq_droop,ωf,xlf,rf,xcf,Kp_u,Ki_u,Kp_i,Ki_i)
```

A node type that applies..

The model has the following internal dynamic variables:
* ``u`` is here an algebraic constraint
* ``θ`` representing the angle of the rotor with respect to the voltage angle ``ϕ``.

# Keyword Arguments
- `Sbase`: "Base apparent power of the grid in VA, should be >0"
- `Srated`: "Rated apperent power of the machine in VA, should be >0"


"""
=#
@DynamicNode MatchingControl(Sbase,Srated,p0set,u0set,Kp_uset,Ki_uset,Kdc,gdc,cdc,xlf,rf,xcf,Kp_u,Ki_u,Kp_i,Ki_i,imax_csa,p_ind) begin
    MassMatrix(m_int =[true,true,true,true,true,true,true])#,false,false,false,false
end begin
    @assert Sbase > 0 "Base apparent power of the grid in VA, should be >0"
    @assert Srated > 0 "Rated apperent power of the machine in VA, should be >0"
    @assert p0set >= 0 "Set point for active power in p.u, should be >0"
    @assert q0set >= 0 "Set point for reactive power in p.u, should be >0"
    @assert u0set > 0 "Set point for voltage in p.u, should be >0"
    @assert Kp_uset >= 0 "Droop constant for active power in p.u, should be >=0"
    @assert Ki_uset >= 0 "Droop constant for reactive power in p.u, should be >=0"
    @assert Kdc >= 0 "Droop constant for current control in p.u, should be >=0"
    @assert gdc >=0 "Conductance of DC-circuit in global p.u., should be >=0"
    @assert cdc >0 "Capacitance of DC-circuit in global p.u., should be >=0"
    @assert xlf >= 0 "filter inductive reactance in p.u., should be >=0"
    @assert rf >= 0 "filter resistance in p.u., should be >=0"
    @assert xcf >= 0 "filter capacitive reactance in p.u., should be >=0"
    @assert Kp_u >= 0 "Proportional gain for voltage control loop, should be >0"
    @assert Ki_u >= 0 "Integral gain for voltage control loop, should be >0"
    @assert Kp_i >= 0 "Proportional gain for current control loop, should be >0"
    @assert Ki_i >= 0 "Integral gain for current control loop, should be >0"
    #@assert imax >= 0 "Current threshold for virtual impedance in p.u., should be >=0"
    #@assert Kvi >= 0 "Gain for virtual impedance in p.u., should be >=0"
    #@assert σXR >= 0 "X/R ratio for for virtual impedance in p.u., should be >=0"
    @assert imax_csa >= 0 "max. current for current saturation algorithm (CSA) in p.u., should be >=0"

end [[θ,dθ],[udc,dudc],[x_uabs,dx_uabs],[e_ud,de_ud],[e_uq,de_uq],[e_id,de_id],[e_iq,de_iq]] begin
    Kp_uset = p[p_ind[1]]
    Ki_uset = p[p_ind[2]]
    ωf_P = p[p_ind[3]]
    ωf_Q = p[p_ind[4]]
    xlf = p[p_ind[5]]
    rf = p[p_ind[6]]
    xcf = p[p_ind[7]]
    Kp_u = p[p_ind[8]]
    Ki_u = p[p_ind[9]]
    Kp_i = p[p_ind[10]]
    Ki_i = p[p_ind[11]]
    #imax = p[p_ind[12]]
    #Kvi = p[p_ind[13]]
    #σXR = p[p_ind[14]]
    imax_csa = p[p_ind[12]]

    #after filter
    umeas = u*(cos(-θ)+1im*sin(-θ))
    udmeas = real(umeas)
    uqmeas = imag(umeas)
    imeas = i*(cos(-θ)+1im*sin(-θ))/(Srated/Sbase) #1im*
    idmeas = real(imeas)
    iqmeas = imag(imeas)
    pmeas = real(u * conj(i))
    qmeas = imag(u * conj(i))

    #before filter
    #The current of the capacitor has to be related, since rf,xlf and xcf are related to Sbase!!!
    idq =  imeas + umeas / (-1im * xcf) / (Srated/Sbase)
    id  = real(idq)
    iq  = imag(idq)
    E = u + (rf + 1im*xlf) * iqd
    p_before_filter = conj(idq) * E
    ix  = id #AC/DC coupling

    #DC current control
    idc = Kdc * (udc - 1.0) + p0set/1.0 + udc*gdc + (p_before_filter - pmeas)

    #DC circuit
    dudc = (idc - gdc * udc - ix) / (cdc)

    #Matching control
    dθ = (udc - 1.0) * 2.0 * pi * 50.0

    Δuabs = u0set - abs(u)
    dx_uabs = Ki_uset * Δuabs
    Uset = x_uabs + Kp_uset * Δuabs

    #Virtual Impedance
    #I_abs = hypot(id,iq)
    #Δi = (I_abs - imax)
    #x_vi = Δi * σXR * Kvi
    #r_vi = x_vi / σXR

    #Δud_vi = IfElse.ifelse(I_abs > imax,r_vi * id - x_vi * iq,0.0)
    #Δuq_vi = IfElse.ifelse(I_abs > imax,r_vi * iq + x_vi * id,0.0)

    #Building voltage reference
    udset = Uset #- Δud_vi
    uqset = 0.0 #- Δuq_vi

    #Voltage control
    de_ud = (udset - udmeas) * Ki_u
    de_uq = (uqset - uqmeas) * Ki_u

    idset = idmeas - uqmeas / xcf + Kp_u * (udset - udmeas) + e_ud
    iqset = iqmeas + udmeas / xcf + Kp_u * (uqset - uqmeas) + e_uq

    #Current saturation algorithm
    iset_abs = hypot(idset,iqset)
    iset_lim = IfElse.ifelse(iset_abs > imax_csa,imax_csa,iset_abs)
    ϕ1 = atan(iqset,idset)
    idset_csa = iset_lim*cos(ϕ1)
    iqset_csa = iset_lim*sin(ϕ1)

    #experimentell
    #anti_windup = IfElse.ifelse(iset_abs > imax_csa,true,false)
    #de_ud = IfElse.ifelse(anti_windup,0.0, (udset - udmeas) * Ki_u)
    #de_uq = IfElse.ifelse(anti_windup,0.0, (uqset - uqmeas) * Ki_u)

    #Current control
    de_id = (idset_csa - id) * Ki_i
    de_iq = (iqset_csa - iq) * Ki_i

    umd = udmeas - iq * xlf + Kp_i * (idset_csa - id) + e_id
    umq = uqmeas + id * xlf + Kp_i * (iqset_csa - iq) + e_iq

    #Coupling with DC voltage
    um_abs = udc * hypot(umd,umq)
    ϕm = atan(umq,umq)
    umds = um_abs * cos(ϕm)
    umqs = um_abs * sin(ϕm)

    #Back transformation to global reference systems
    um = umds + 1im * umqs
    um = um*(cos(θ)+1im*sin(θ))
    umd = real(um)
    umq = imag(um)

    idq = id + 1im*iq
    idq = idq*(cos(θ)+1im*sin(θ)) #-1im*
    id = real(idq) #* (Srated/Sbase)
    iq = imag(idq) #* (Srated/Sbase)

    #Voltage equations for filter
    u0d =  umd - rf * id + xlf * iq
    u0q =  umq - rf * iq - xlf * id
    u0 = u0d + 1im * u0q

    du = u - u0 #algebraic constraint
end

export MatchingControl
