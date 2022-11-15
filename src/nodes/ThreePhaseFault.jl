
@DynamicNode ThreePhaseFault(p_ind) begin
    MassMatrix()#m_int =[false,false]
end  begin
end [] begin#[[P,dP],[Q,dQ]]
    s = u * conj(i)
    rfault = p[p_ind[1]]
    xfault = p[p_ind[2]]
    yfault = 1.0 / (rfault + 1im * xfault)
    
    du = conj(-yfault * u) * u - s 
    #dP = P - real(conj(i)*u)
    #dQ = Q - imag(conj(i)*u)
end

export VoltageDependentLoad