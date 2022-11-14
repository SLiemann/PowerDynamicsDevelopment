
@DynamicNode ThreePhaseFault(Yfault,p_ind) begin
    MassMatrix()#m_int =[false,false]
end  begin
end [] begin#[[P,dP],[Q,dQ]]
    s = u * conj(i)
    active = p[p_ind]
    du = active .* conj(-Yfault * u) * u - s 
    #dP = P - real(conj(i)*u)
    #dQ = Q - imag(conj(i)*u)
end

export VoltageDependentLoad