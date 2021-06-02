#= Sebastian Liemann, ie3 TU Dortmund, based on F. Milano, Power System Modelling and Scripting, Springer Verlag, 2010
@doc doc"""
```Julia
GridFormingConverter(Sbase,Srated,p0set,u0set,Kp_droop,Kq_droop,ωf,lf,rf,cf,Kp_u,Ki_u,Kp_i,Ki_i)
```

A node type that applies the 6th-order (sometimes also called 4th-order if ω and δ are not counted)
synchronous machine model which is implemented according to
F. Milano, "Power System Modelling and Scripting", Springer Verlag, 2010
The main equations are on page 331, cf. Table 15.2 Model 6.b

The model has the following internal dynamic variables:
* ``u`` is here an algebraic constraint
* ``θ`` representing the angle of the rotor with respect to the voltage angle ``ϕ``.

# Keyword Arguments
- `Sbase`: "Base apparent power of the grid in VA, should be >0"
- `Srated`: "Rated apperent power of the machine in VA, should be >0"


"""
=#
@DynamicNode GridFormingConverter(Sbase,Srated,p0set,u0set,Kp_droop,Kq_droop,ωf,lf,rf,cf,Kp_u,Ki_u,Kp_i,Ki_i) begin
    MassMatrix(m_int =[true,true,true,true,true,true,true])
end begin
    @assert Sbase > 0 "Base apparent power of the grid in VA, should be >0"
    @assert Srated > 0 "Rated apperent power of the machine in VA, should be >0"
    @assert u0set > 0 "Set point for voltage in p.u, should be >0"
    @assert Kp_droop >= 0 "Droop constant for active power in p.u, should be >=0"
    @assert Kq_droop >= 0 "Droop constant for reactive power in p.u, should be >=0"
    @assert ωf > 0 "Cut-off angular filter of meauserement (both Q/P) in rad, should be >0"
    @assert lf >= 0 "filter inductance in p.u., should be >=0"
    @assert rf >= 0 "filter reactance in p.u., should be >=0"
    @assert cf >= 0 "filter capacitance in p.u., should be >=0"
    @assert Kp_u > 0 "Proportional gain for voltage control loop, should be >0"
    @assert Ki_u > 0 "Integral gain for voltage control loop, should be >0"
    @assert Kp_i > 0 "Proportional gain for current control loop, should be >0"
    @assert Ki_i > 0 "Integral gain for current control loop, should be >0"

    ω0 = 2.0*pi*50.0

end [[θ,dθ],[dω,ω],[dQ,Q],[de_ud,e_ud],[de_uq,e_uq][de_id,e_id],[de_iq,e_iq]] begin
    umeas = 1im*u*(cos(-θ)+1im*sin(-θ))
    udmeas = real(umeas)
    uqmeas = imag(umeas)
    imeas = 1im*i*(cos(-θ)+1im*sin(-θ))/(Srated/Sbase)
    idmeas = real(imeas)
    iqmeas = imag(imeas)
    pmeas = real(u * conj(i))
    qmeas = imag(u * conj(i))

    id  = imeasd + uqmeas * ω0 * cf
    iq  = imeasq - udmeas * ω0 * cf

    #Droop control
    dω = (p0set - pmeas) * Kp_droop * ωf - ω * ωf
    dθ = ω

    dQ = qmeas * ωf - Q * ωf
    Uset = (qmeas - Q) * Kq_droop + u0set

    #Virtual Impedance
    #TO-DO

    udset = real(Uset * (cos(θ) + 1im * sin(θ)))
    uqset = imag(Uset * (cos(θ) + 1im * sin(θ)))
    #Voltage control
    de_ud = udset - udmeas
    de_uq = uqset - uqmeas

    idset = idmeas + uqmeas * ω0 * cf + Kp_u * de_ud + Ki_u * e_ud
    iqset = iqmeas - udmeas * ω0 * cf + Kp_u * de_uq + Ki_u * e_uq

    #Current control
    de_id = idset - idmeas
    de_iq = iqset - iqmeas

    umd = udmeas + iqmeas * ω0 * lf + Kp_i * de_id + Ki_i * e_id
    umq = uqmeas - idmeas * ω0 * lf + Kp_i * de_iq + Ki_i * e_iq

    u0d =  umd + rf * id - ω0 * lf * iq
    u0q =  umq + rf * iq + ω0 * lf * id
    u0 = u0d + 1im * u0q

    du = u - -1im*u0*(cos(θ)+1im*sin(θ)) #algebraic constraint
end

export GridFormingConverter
