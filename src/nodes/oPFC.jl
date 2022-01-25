#=
@doc doc"""
```Julia
SimpleRecoveryLoad(P0, Q0, Pd, Pt, Qt,Tp, Tq)
```
A node type that represents a simple exponential recovery load model.
Based on:
I. Hisken, M. Pai: "Trajectory Sensitivity Analysis of Hybrid Systems",
IEEE TRANSACTIONS ON CIRCUITS AND SYSTEMS, vol. 47, February 2000

# Keyword Arguments
- `P0`: Active power load demand [pu]
- `Q0`: Reactive power load demand [pu]
- `Pt`: Inital transient active power change [pu]
- `Qt`: Inital transient reactive power change  [pu]
- `Tp`: Load recovery constant p-axis [s]
- `Tq`: Load recovery constant q-axis [s]
"""
=#
@DynamicNode oPFC(Cd, Pdc, Ulow, Qn, t0,ϵ,p_ind)  begin
    MassMatrix(m_int = [true, true, true, true,true, true])
end  begin
    @assert Cd > 0.0 "Cd should be >0"
    #@assert Pdc > 0.0 "Pdc should be >0"
    @assert Ulow > 0.0 "Ulow should be >0"
    @assert t0 >= 0.0 "t0 should be >0"
    @assert ϵ >= 0.0 "ϵ should be >0"

end [[Uc, dUc],[qstate, dqstate],[P, dP],[Q, dQ],[tsim,dtsim],[Iload, dIload]] begin
    Cd = p[p_ind[1]]
    Pdc = p[p_ind[2]]
    Ulow = p[p_ind[3]]
    Qn = p[p_ind[4]]
    t0 = p[p_ind[5]]
    ϵ = p[p_ind[6]]

    s = u*conj(i)
    u_abs = abs(u)

    if 0.9 <= qstate <= 1.1 #equal
        dUc = 0.0
        Psoll = Pdc
        if u_abs <= 0.0
             u_abs = Inf
         end
        Qsoll = Qn/(u_abs^0.9)
        dtsim = 0.0
        ic = 0.0
    elseif 1.9 <= qstate <= 2.1 # low grid voltage
       dUc = -1.0 /(Cd*Uc)
       Psoll = 0.0
       Qsoll = 0.0
       dtsim = 0.0
       ic = 0.0
    elseif 2.9 <= qstate <= 3.1 #high grid voltage
       ic = 2*pi*50*u_abs*cos(2*pi*50*tsim+t0)*Cd
       dUc = ic/Cd
       Psoll = Uc*ic/1.2  + 1.0
       Qsoll = -Uc*ic/2 + Qn #-u[2]*dUc*Cd*0.5 (Psoll-1)*sin(PH)
       dtsim = 1.0
     elseif 3.9 <= qstate <= 4.1 #disconnected
       dUc = 0.0
       Psoll = 0.0
       Qsoll = 0.0
       dtsim = 0.0
       ic = 0.0
    else #should not be reached
        dUc = 0.0
        Psoll = 0.0
        Qsoll = 0.0
        ic = 0.0
    end

    dqstate = 0.0
    dP = (Psoll - P) / 0.01
    dQ = (Qsoll - Q) / 0.01
    dIload= ic

    du = P + im*Q - s
end

export oPFC
