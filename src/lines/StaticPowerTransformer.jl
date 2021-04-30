"""
```Julia
    StaticPowerTransformer(from,to,S_r,U_r,uk,XR_ratio,i0,Pv0,Sbase,Ubase,tap_side,tap_pos,tap_inc)
```

Transformer based on typical equipment parameters

# Arguments

- `from` : start node
- `to` : end node
- `S_r`: rated apparent Power in [W]
- `U_r`: rated RMS line-line voltage in [V]
- `uk`: short-circuit voltage in [%], e.g 15% => 0.15
- `XR_ratio`:X/R ratio of uk
- `i0`: no load current (core losses) in [%], e.g 6 % => 0.06
- `Pv0`: iron core losses in [W]
- `Sbase`: Base apprarent power in [W]
- `Ubase`: Base voltage in [V]
- `tap_side`: Side at which side the Transformer should be stepped, either LV or HV
- `tap_pos`: current tap position (integer), zero is neutral position
- `tap_inc`: voltag increase per tap in [%], e.g. 1 %

"""
@Line StaticPowerTransformer(from,to,S_r,U_r,uk,XR_ratio,i0,Pv0,Sbase,Ubase,tap_side,tap_pos,tap_inc) begin
    üLV    = 1.0
    üHV    = 1.0
    if tap_side == "LV"
        üLV = 1.0 / (1.0 + tap_pos * tap_inc / 100.0)
    elseif tap_side == "HV"
        üHV = 1.0 / (1.0 + tap_pos * tap_inc / 100.0)
    else
        error("Can not interprete tap_side (HV/LV): $tap_side")
    end

    #I_r   = S_r/(U_r*sqrt(3)) #rated current
    I_r   = S_r/(Ubase*sqrt(3)) #rated current

    #Calculatiing leakage reactance Xa and winding resistance Ra
    Ra = 0.
    Xa = 0.
    if XR_ratio == 0
        Ra = uk*U_r/(sqrt(3)*I_r);
        Xa = 0.
    elseif XR_ratio == Inf
        Ra = 0.
        Xa = uk*U_r/(sqrt(3)*I_r);
    else
        Ra = sqrt(((uk*U_r/(I_r*sqrt(3)))^2)/(1+XR_ratio^2));
        Xa = XR_ratio*Ra;
    end
    Ya = 1/(0.5*(Ra+1im*Xa)) #it is assumed that losses are equally (0.5) distributed over both sides
    Ybs = Ya; #Ybs should be changed, if losses are not equally distributed

    #Calculating magnetising reactance Xm and core resistance Rfe from iron losess
    Zm  = U_r/(sqrt(3)*i0/100.0*I_r) #no-load currents depends on complete magnetising impedance
    Rfe = U_r*U_r/(Pv0)
    Xm  = Inf
    if i0 != 0
        Xm = 1/(sqrt(1/Zm^2 - 1/Rfe^2))
    end
    Ym  = 1/Rfe + 1/(1im*Xm)

    Ybase = 1/(Ubase^2/Sbase) #this could be changed, if global base values are available
    Y     = 1.0/(Ya+Ybs+Ym)./Ybase
    Y = Y.*PiModel(Ya*Ybs, Ya*Ym, Ybs*Ym,üHV,üLV)
    # If current is flowing away from the source, it is negative at the source.
    voltage_vector = [source_voltage,destination_voltage]
    current_vector = Y * voltage_vector
end

export StaticPowerTransformer
