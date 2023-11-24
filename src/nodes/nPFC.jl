@DynamicNode nPFC(Cd, Pdc,p_offset,q_offset,p_ind)  begin
    MassMatrix(m_int = [true, true, true, true,true, true,true])
end  begin
    @assert Cd > 0.0 "Cd should be >0"
    @assert Pdc > 0.0 "Pdc should be >0"
    #@assert Voff > 0.0 "Ulow should be >0"

end [[vofft2, dvofft2],[tsum,dtsum],[ton,dton],[toff,dtoff],[p1,dp1],[q1,dq1],[Vabstoff,dVabstoff]] begin
    Cd = p[p_ind[1]]
    Pdc = p[p_ind[2]]
    #Vlow = p[p_ind[3]]

    s = u*conj(i)
    q_off = q_offset*(abs(u)^2)
    p_off = p_offset*(abs(u))
    
    dvofft2 = 0.0
    dtsum = 0.0
    dton = 0.0
    dtoff = 0.0
    dp1 = 0.0
    dq1 = 0.0
    dVabstoff = 0.0

    du = (p1+p_off) + im*(q1+q_off) - s
end

export nPFC
