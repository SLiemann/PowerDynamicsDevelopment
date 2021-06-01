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
@DynamicNode SimpleRecoveryLoad(P0, Q0, Pt, Qt, Tp, Tq)  begin
    MassMatrix(m_int = [true, true])
end  begin
    @assert Tp > 0 "Load recovery constant should be >0"
    @assert Tq > 0 "Load recovery constant should be >0"

end [[x_p, dx_p],[x_q, dx_q]] begin
    s = u*conj(i)
    Pd = real(s)
    Qd = imag(s)

    dx_p = (1/Tp)*(P0 -Pd)
    dx_q = (1/Tq)*(Q0 -Qd)

    du = x_p + Pt*abs(u)^2.0 + im*(x_q + Qt*abs(u)^2.0) - s
end

export SimpleRecoveryLoad
