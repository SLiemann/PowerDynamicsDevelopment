@DynamicNode SynchronousMachineGENSAL(Sbase, Srated, D, H, P, E_fd, R_a, x_d, x_q, x_ds, x_dss, x_qss, x_l, T_d0s, T_d0ss, T_q0ss, S10, S12) begin
MassMatrix(m_u = false, m_int =[true,true,true,true,true])
end begin
    # based on: PowerWorld Corporation
    a = 0.2; b = 2*S12-1.2*2*S10; c = 1.2*S10^2-S12^2
    A = (-b + sqrt(b^2 -4*a*c))/(2*a)
    B = 1/(A^2 - 2*A*S10 + S10^2)
    sf(x) = B*(x-A)^2
    ω0 = 2*pi*50

end [[ω, dω], [δ, dδ], [E_qs, dE_qs], [ψ_ds, dψ_ds], [ψ_qss, dψ_qss]] begin

    I_c = i*(cos(δ-pi/2) - 1im*sin(δ-pi/2)) # from Sauer/Pai
    I_d = real(I_c)
    I_q = imag(I_c)

    print("\n ω:  ", ω)
    print("\n δ:  ", δ)


#### Internal Model
    
    dψ_qss = 1/T_q0ss*(-I_q*(x_q-x_qss) - ψ_qss)
    dψ_ds = 1/T_d0ss*(E_qs - I_d*(x_ds-x_l) - ψ_ds)
    dE_qs = 1/T_d0s*(E_fd - E_qs*(1+sf(E_qs)) + (x_d-x_ds)*((x_ds-x_dss)/(x_ds-x_l)^2*(E_qs-ψ_ds-I_d*(x_ds-x_l)) + I_d))

    ψ_dss = ψ_ds*(x_ds-x_dss)/(x_ds-x_l) + E_qs*(x_dss-x_l)/(x_ds-x_l)

#### Electrical Torque

    #ψ_q = ψ_qss - I_q*x_dss
    #ψ_d = ψ_dss - I_d*x_dss
    #T_elec = ψ_d*I_q - ψ_q*I_d
    T_elec = ψ_dss*I_q - ψ_qss*I_d

    dδ = ω*ω0
    dω = 1/(2*H)*((P - D*ω)/(1+ω) - T_elec)

#### Network Interface Equations
# Taken from Powerworld


    v_d = -ψ_qss*(1+ω)
    v_q = ψ_dss*(1+ω)

    v_dq = v_d + 1im*v_q

    i_d = real(v_dq/(R_a+1im*x_dss))
    i_q = imag(v_dq/(R_a+1im*x_dss))

    i_dq = i_d+1im*i_q

    i_r = real(i_dq*(cos(δ-pi/2)+1im*sin(δ-pi/2)))
    i_i = imag(i_dq*(cos(δ-pi/2)+1im*sin(δ-pi/2)))

    u_ll = (i_r+1im*i_i)*(R_a+1im*x_qss)
    
    u_network = u_ll - I_c*(u_ll/(i_r+1im*i_i))

    du = u - u_network

    #v_d = -R_a * i_d + (ω + 1.) * (ψ_dss  + x_qss * i_q)
    #v_q = -R_a * i_q + (ω + 1.) * (ψ_qss - x_dss * i_d)

    #v  = v_d + 1im*v_q
    #du = u - -1im*v*(cos(δ)+1im*sin(δ)) #algebraic constraint

#### Network Interface Equations
#=
    v_d = -ψ_qss*(1+ω)
    v_q = ψ_dss*(1+ω)

    v  = v_d + 1im*v_q
    du = u - -1im*v*(cos(δ)+1im*sin(δ)) #algebraic constraint
=#
end
export SynchronousMachineGENSAL