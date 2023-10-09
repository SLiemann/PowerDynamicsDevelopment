@DynamicNode WeccPeLoad(P, Q, Vd1, Vd2, p_ind) begin
    MassMatrix(m_int =[false,false])
end  begin
    @assert Vd1 > 0
    @assert Vd2 > 0
end [[P0,dp],[Q0,dq]] begin
    s = u * conj(i)

    frac = p[p_ind[1]]

    Pv = P * frac
    Qv = Q * frac

    du = complex(Pv, Qv) - s

    dp = P0 - Pv
    dq = Q0 - Qv
end

export WeccPeLoad
