# Sebastian Liemann, ie3 TU Dortmund, based on F. Milano, Power System Modelling and Scripting, Springer Verlag, 2010
@doc doc"""
```Julia
SixOrderMarcanatoMachine(H, P, D, Ω, E_f, R_a,T_ds,T_qs,T_dss,T_qss,X_d,X_q,X_ds,X_qs,X_dss,X_qss,T_AA)
```

A node type that applies the 4th-order synchronous machine model with frequency/angle and voltage dynamics,
which is implemented according to P. Sauer, "Power System Dynamics and Stability".
For an illustration of a synchronous machine schematic see P. Sauer, Fig. 3.1 on p. 25.

Additionally to ``u``, it has the internal dynamic variables
* ``ω`` representing the frequency of the rotator relative to the grid frequency ``Ω``, i.e. the real frequency ``ω_r`` of the rotator is given as ``\omega_r = \Omega + \omega`` and
* ``θ`` representing the relative angle of the rotor with respect to the voltage angle ``ϕ``.

# Keyword Arguments
- `H`: shaft inertia constant (given in [s]), defined according to P. Sauer, p. 33, eq. (3.60)
- `P`: active (real) power output, also called the mechanical torque applied to the shaft, given in [pu]
- `D`: damping coefficient (given in [s], see P. Sauer, eq. (5.156) where the damping torque is equal `Dω`)
- `Ω`: rated frequency of the power grid, often ``2π⋅50Hz``
- `E_f`: field voltage in [pu.]
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

"""
@DynamicNode SixOrderMarcanatoMachine(H, P, D, Ω, E_f, R_a,T_ds,T_qs,T_dss,T_qss,X_d,X_q,X_ds,X_qs,X_dss,X_qss,T_AA) begin
    @assert H > 0 "inertia (H) should be >0"
    @assert D >= 0 "damping (D) should be >=0"
    @assert E_f >= 0 "Field Voltage (E_f) should be >=0"
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

    #Converstion of short-circuit time constants to open-loop time constants
    T_d0s = T_ds*(X_d/X_ds)
    T_q0s = T_qs*(X_q/X_qs)
    T_d0ss = T_dss*(X_ds/X_dss)
    T_q0ss = T_qss*(X_qs/X_qss)

    #Auxillary variables
    γ_d = T_d0ss * X_dss * (X_d - X_ds) / (T_d0s * X_ds)
    γ_q = T_q0ss * X_qss * (X_q - X_qs) / (T_q0s * X_qs)

    Ω_H = Ω / (2*H)

end [[θ,dθ],[ω, dω],[e_ds, de_ds],[e_qs, de_qs],[e_dss, de_dss],[e_qss, de_qss]] begin
    i_c = 1im*i*exp(-1im*θ)
    #e_c = 1im*u*exp(-1im*θ)
    #p = real(u * conj(i))
    #e_dss = real(e_c)
    #e_qss = imag(e_c)
    i_d = real(i_c)
    i_q = imag(i_c)

    dθ = ω
    de_qs = (1 / T_d0s) * (- e_qs - (X_d - X_ds - γ_d) * i_d + (1 - T_AA/T_d0s) * E_f)
    de_ds = (1 / T_q0s) * (- e_ds + (X_q - X_qs - γ_q) * i_q)

    de_qss = (1 / T_d0ss) * (- e_qss + e_qs - (X_ds - X_dss + γ_d) * i_d + (T_AA/T_d0s) * E_f)
    de_dss = (1 / T_q0ss) * (- e_dss + e_ds + (X_qs - X_qss + γ_q) * i_q)

    v_d = -R_a * i_d + e_dss + X_qss * i_q
    v_q = -R_a * i_q + e_qss - X_dss * i_d

    v  = v_d + 1im*v_q
    du = -1im*v*exp(1im*θ) + u*1im*ω
    p  = (v_q + R_a * i_q) * i_q + (v_d + R_a * i_d) * i_d
    dω = (P - D*ω - p)*Ω_H
end

export SixOrderMarcanatoMachine
