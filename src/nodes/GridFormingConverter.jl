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
@DynamicNode GridFormingConverter(Sbase,Srated,p0set,q0set,u0set,Kp_droop,Kq_droop,ωf,xlf,rf,xcf,Kp_u,Ki_u,Kp_i,Ki_i,imax,Kvi,σXR) begin
    MassMatrix(m_int =[true,true,true,true,true,true,true,false,false])
end begin
    @assert Sbase > 0 "Base apparent power of the grid in VA, should be >0"
    @assert Srated > 0 "Rated apperent power of the machine in VA, should be >0"
    @assert p0set >= 0 "Set point for active power in p.u, should be >0"
    @assert q0set >= 0 "Set point for reactive power in p.u, should be >0"
    @assert u0set > 0 "Set point for voltage in p.u, should be >0"
    @assert Kp_droop >= 0 "Droop constant for active power in p.u, should be >=0"
    @assert Kq_droop >= 0 "Droop constant for reactive power in p.u, should be >=0"
    @assert ωf > 0 "Cut-off angular filter of meauserement (both Q/P) in rad, should be >0"
    @assert xlf >= 0 "filter inductive reactance in p.u., should be >=0"
    @assert rf >= 0 "filter resistance in p.u., should be >=0"
    @assert xcf >= 0 "filter capacitive reactance in p.u., should be >=0"
    @assert Kp_u >= 0 "Proportional gain for voltage control loop, should be >0"
    @assert Ki_u >= 0 "Integral gain for voltage control loop, should be >0"
    @assert Kp_i >= 0 "Proportional gain for current control loop, should be >0"
    @assert Ki_i >= 0 "Integral gain for current control loop, should be >0"
    @assert imax >= 0 "Current threshold for virtual impedance in p.u., should be >=0"
    @assert Kvi >= 0 "Gain for virtual impedance in p.u., should be >=0"
    @assert σXR >= 0 "X/R ratio for for virtual impedance in p.u., should be >=0"

end [[θ,dθ],[ω,dω],[Qm,dQm],[e_ud,de_ud],[e_uq,de_uq],[e_id,de_id],[e_iq,de_iq],[umabs,dumabs],[umangle,dumangle]] begin
    umeas = u*(cos(-θ)+1im*sin(-θ))
    udmeas = real(umeas)
    uqmeas = imag(umeas)
    imeas = 1im*i*(cos(-θ)+1im*sin(-θ))/(Srated/Sbase)
    idmeas = real(imeas)
    iqmeas = imag(imeas)
    pmeas = real(u * conj(i))
    qmeas = imag(u * conj(i))

    idq = i + u / (-1im * xcf) / (Srated/Sbase)
    idq = 1im*idq*(cos(-θ)+1im*sin(-θ))
    id  = real(idq)
    iq  = imag(idq)

    #Droop control
    dω = (p0set - pmeas) * Kp_droop * ωf - ω * ωf
    dθ = ω

    dQm = qmeas * ωf - Qm * ωf
    Uset = (q0set - Qm) * Kq_droop + u0set

    #Virtual Impedance
    I_abs = hypot(id,iq)
    Δud_vi = 0.0
    Δuq_vi = 0.0
    if I_abs >= imax
        Δi = I_abs - imax
        #r_vi = Δi * Kvi
        #x_vi = r_vi * σXR
        x_vi = Δi * σXR * Kvi
        r_vi = x_vi / σXR
        Δud_vi = r_vi * id - x_vi * iq
        Δuq_vi = r_vi * iq + x_vi * id
    end

    udset = Uset - Δud_vi #real(Uset * (cos(θ) + 1im * sin(θ)))
    uqset = 0.0 - Δuq_vi #imag(Uset * (cos(θ) + 1im * sin(θ)))
    #Voltage control
    de_ud = udset - udmeas
    de_uq = uqset - uqmeas

    idset = idmeas - uqmeas / xcf + Kp_u * de_ud + Ki_u * e_ud
    iqset = iqmeas + udmeas / xcf + Kp_u * de_uq + Ki_u * e_uq

    #Current control
    de_id = idset - id
    de_iq = iqset - iq

    umd = udmeas - iq * xlf + Kp_i * de_id + Ki_i * e_id
    umq = uqmeas + id * xlf + Kp_i * de_iq + Ki_i * e_iq

    #Back transformation to global reference systems
    um = umd + 1im * umq
    um = um*(cos(θ)+1im*sin(θ))
    umd = real(um)
    umq = imag(um)

    dumabs = umabs - abs(um)  #for output
    dumangle = umangle - angle(um)/pi*180.0 #for output

    idq = id + 1im*iq
    idq = -1im*idq*(cos(θ)+1im*sin(θ))
    id = real(idq)
    iq = imag(idq)

    #Voltage equations for filter
    u0d =  umd - rf * id + xlf * iq
    u0q =  umq - rf * iq - xlf * id
    u0 = u0d + 1im * u0q

    du = u - u0 #algebraic constraint
end

export GridFormingConverter
