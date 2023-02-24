using IfElse
#= Sebastian Liemann, ie3 TU Dortmund, based on F. Milano, Power System Modelling and Scripting, Springer Verlag, 2010
@doc doc"""
```Julia
SixOrderMarconatoMachineAVROEL(Sbase,Srated,H, P, D, Ω, R_a,T_ds,T_qs,T_dss,T_qss,X_d,X_q,X_ds,X_qs,X_dss,X_qss,T_AA,V0,Ifdlim,L1,G1,Ta,Tb,G2,L2)
```

A node type that applies the 6th-order (sometimes also called 4th-order if ω and δ are not counted)
synchronous machine model which is implemented according to
F. Milano, "Power System Modelling and Scripting", Springer Verlag, 2010
The main equations are on page 331, cf. Table 15.2 Model 6.b

Also the model includes an Automatic Voltage Regulator and Overexcitation Limiter (OEL) from (see Fig. 2.6).
PES Technical Report 19 : "Test Systems for Voltage Stability Analysis and Security Assessment", 2015.

The model has the following internal dynamic variables:
* ``u`` is here an algebraic constraint
* ``e_ds`` transient magnetic state in d-axis
* ``e_qs`` transient magnetic state in q-axis
* ``e_dss`` subtransient magnetic state in d-axis
* ``e_qss`` subtransient magnetic state in q-axis
* ``ω`` representing the frequency of the rotator (not relative)
* ``θ`` representing the angle of the rotor with respect to the voltage angle ``ϕ``
* ``ifd`` field current - algebraic constraint.
* ``Timer`` Timer inside OEL.
* ``x1`` state inside transient gain reduction (PDT1).
* ``E_f`` Field voltage, output of Exciter.

# Keyword Arguments
- `Sbase`: "Base apparent power of the grid in VA, should be >0"
- `Srated`: "Rated apperent power of the machine in VA, should be >0"
- `H`: shaft inertia constant, given in [s],
- `P`: active (real) power output, also called the mechanical torque applied to the shaft, given in [pu]
- `D`: damping coefficient, given in [s] (here D(ω-1.) is used)
- `Ω`: rated frequency of the power grid, often ``2π⋅50Hz``
- `R_a` : armature resistance, given in [pu]
- `T_ds` : short-circuit transient time constant of d-axis, given in [s]
- `T_qs` : short-circuit transient time constant of q-axis, given in [s]
- `T_dss`: short-circuit subtransient transient time constant of d-axis, given in [s]
- `T_qss`: short-circuit subtransient transient time constant of s-axis, given in [s]
- `X_d`: synchronous reactance of d-axis, given in [pu]
- `X_q`: synchronous reactance of q-axis, given in [pu]
- `X_ds`: transient reactance of d-axis, given in [pu]
- `X_qs`: transient reactance of d-axis, given in [pu]
- `X_dss`: subtransient reactance of d-axis, given in [pu]
- `X_qss`: subtransient reactance of d-axis, given in [pu]
- `T_AA` : additional leakage time constant in d-axis, given in [s]
- `V0` : Reference Voltage of AVR
- `Ifdlim` : Maximum field current
- `L1` : Lower Limit of Timer-Integrator
- `G1` : Gain before Transient Gain Reduction (PDT1)
- `Ta` : Nominator Time Constant of PDT1
- `Tb` : Nominator Time Constant of PDT1
- `G2` : Gain of PT1
- `L2` : Upper Limit of Anti-Windup Integrator inside of the PT1

"""
=#
@DynamicNode SixOrderMarconatoMachineAVROEL(Sbase,Srated,H, P, D, Ω, R_a,T_ds,T_qs,T_dss,T_qss,X_d,X_q,X_ds,X_qs,X_dss,X_qss,T_AA,V0,Ifdlim,L1,G1,Ta,Tb,G2,L2) begin
    MassMatrix(m_int =[true,true,true,true,true,true,false,true,true,true])
end begin
    @assert Sbase > 0 "Base apparent power of the grid in VA, should be >0"
    @assert Srated > 0 "Rated apperent power of the machine in VA, should be >0"
    @assert H > 0 "inertia (H) should be >0"
    @assert P >= 0 "Active power (P) should be >=0"
    @assert D >= 0 "damping (D) should be >=0"
    #@assert E_f >= 0 "Field Voltage (E_f) should be >=0"
    @assert R_a >= 0 "armature resistance (R_a) should be >=0"
    @assert T_ds > 0 "time constant of d-axis (T_ds) should be >0"
    @assert T_qs > 0 "time constant of q-axis (T_qs) should be >0"
    @assert T_dss > 0 "time constant of d-axis (T_dss) should be >0"
    @assert T_qss > 0 "time constant of q-axis (T_qss) should be >0"
    @assert X_d >= 0 "reactance of d-axis (X_d) should be >=0"
    @assert X_q >= 0 "reactance of q-axis (X_q) should be >=0"
    @assert X_ds > 0 "transient reactance of d-axis (X_ds) should be >0"
    @assert X_qs > 0 "transient reactance of q-axis (X_qs) should be >0"
    @assert X_dss > 0 "subtransient reactance of d-axis (X_dss) should be >0"
    @assert X_qss > 0 "subtransient reactance of q-axis (X_qss) should be >0"
    @assert T_AA >= 0 "additional leakage time constant of d-axis (T_AA) should be >=0"

    #AVR & OEL
    @assert V0 >= 0 "Reference Voltage of AVR"
    @assert Ifdlim >= 0 "Maximum field current"
    @assert L1 < 0 "Lower Limit of Timer-Integrator"
    @assert G1 >= 0 "Gain before Transient Gain Reduction (PDT1)"
    @assert Ta >= 0 "Nominator Time Constant of PDT1"
    @assert Tb >= 0 "Nominator Time Constant of PDT1"
    @assert G2 >= 0 "Gain of PT1 (Exciter)"
    @assert L2 >= 0 "Upper Limit of Anti-Windup Integrator inside of the PT1 (Exciter)"

    #Converstion of short-circuit time constants to open-loop time constants
    T_d0s = T_ds*(X_d/X_ds)
    T_q0s = T_qs*(X_q/X_qs)
    T_d0ss = T_dss*(X_ds/X_dss)
    T_q0ss = T_qss*(X_qs/X_qss)

    #Auxillary variables
    γ_d = T_d0ss * X_dss * (X_d - X_ds) / (T_d0s * X_ds)
    γ_q = T_q0ss * X_qss * (X_q - X_qs) / (T_q0s * X_qs)

end [[θ,dθ],[ω, dω],[e_ds, de_ds],[e_qs, de_qs],[e_dss, de_dss],[e_qss, de_qss],[ifd,difd],[timer,dtimer],[x1,dx1],[E_f,dE_f]] begin
    #i_c = 1im*i*exp(-1im*θ)
    i_c = 1im*i*(cos(-θ)+1im*sin(-θ))/(Srated/Sbase)
    i_d = real(i_c)
    i_q = imag(i_c)
    pe = real(u * conj(i))

    de_ds = (1 / T_q0s) * (- e_ds + (X_q - X_qs - γ_q) * i_q)
    de_qs = (1 / T_d0s) * (- e_qs - (X_d - X_ds - γ_d) * i_d + (1 - T_AA/T_d0s) * E_f)

    de_dss = (1 / T_q0ss) * (- e_dss + e_ds + (X_qs - X_qss + γ_q) * i_q)
    de_qss = (1 / T_d0ss) * (- e_qss + e_qs - (X_ds - X_dss + γ_d) * i_d + (T_AA/T_d0s) * E_f)

    v_d = -R_a * i_d + (ω + 1.) * (e_dss  + X_qss * i_q)
    v_q = -R_a * i_q + (ω + 1.) * (e_qss - X_dss * i_d)

    v  = v_d + 1im*v_q
    #du = u - -1im*v*exp(1im*θ) #algebraic constraint
    du = u - -1im*v*(cos(θ)+1im*sin(θ)) #algebraic constraint

    Te = ((v_q + R_a * i_q) * i_q + (v_d + R_a * i_d) * i_d) / (ω + 1.0) 

    dθ = Ω * 2*pi * ω
    dω = ((P - D * ω)/(ω + 1.0) - Te) / (2*H)

    #Field current in a non-reciprocal system, otherwise would be: ifd = (E_f - T_d0s * de_qs) / (X_d - X_l)
    difd = ifd - (E_f - T_d0s * de_qs) / (X_d - 0.15)  #algebraic constraint for output

    #AVR error
    V_error = V0 - abs(u)
    #OEL
    ifd_error = ifd - Ifdlim
    array_out  = IfElse.ifelse(ifd_error > 0.0,ifd_error,IfElse.ifelse(ifd_error >= -0.1,0.0,-1.0))
    dtimer =IfElse.ifelse(timer<=L1,IfElse.ifelse(array_out<0.0,0.0,array_out),array_out)
    switch_output = IfElse.ifelse(timer<0.0,V_error,-ifd_error)

    #AVR - Transient Gain Reducution & Exciter
    min_out = min(V_error,switch_output)
    dx1 = (min_out*G1 - x1) / Tb
    PDT1_out = x1 + dx1*Ta
    e = G2*(PDT1_out - E_f)

    lowlimit  = IfElse.ifelse(E_f<=0.0,IfElse.ifelse(e<0.0,true,false),false)
    highlimit = IfElse.ifelse(E_f>= L2,IfElse.ifelse(e>0.0,true,false),false)
    dE_f = IfElse.ifelse(lowlimit==true,0.0,IfElse.ifelse(highlimit==true,0.0,e))
end

export SixOrderMarconatoMachineAVROEL
