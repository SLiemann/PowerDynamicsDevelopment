"""
```Julia
    DyanmicpowerTransformer(from,to,S_r,U_HV,U_LV,uk,XR_ratio,i0,Pv0,Sbase,Ubase,tap_side,tap_pos,tap_inc)
```

Transformer based on typical equipment parameters

# Arguments

- `from` : start node
- `to` : end node
- `S_r`: rated apparent Power in [W]
- `U_HV`: rated RMS line-line voltage at high-voltage side in [V]
- `U_LV`: rated RMS line-line voltage at low-voltage side in [V]
- `uk`: short-circuit voltage in [%], e.g 15% => 0.15
- `XR_ratio`:X/R ratio of uk
- `i0`: no load current (core losses) in [%], e.g 6 % => 0.06
- `Pv0`: iron core losses in [W]
- `Sbase`: Base apprarent power in [W]
- `Ubase`: Base voltage in [V]
- `tap_side`: Side at which side the Transformer should be stepped, either LV or HV
- `tap_pos`: current tap position (integer)
- `tap_inc`: voltag increase per tap in [%], e.g. 1 %
- `tap_max`: maximum tap position (both for pos. and neg. tapping)
- `v_ref`: reference voltage for voltage regulation in p.u.
- `v_dead`: deadband around v_ref in p.u.
- `tap_time`: time after the transformer taps, after voltage deadband is left
"""
@Line DynamicPowerTransformer(from,to,S_r,U_HV,U_LV,uk,XR_ratio,i0,Pv0,Sbase,Ubase,tap_side,tap_pos,tap_inc,tap_max,v_ref,v_dead,tap_time) begin
end [[timer, dtimer][tap,dtap]]begin
    abs_uHV = abs(v_s[1]+1im*v_s[2])
    abs_uLV = abs(v_d[1]+1im*v_d[2])

    if abs_uLV < v_ref - v_dead || abs_uHV < v_ref + v_dead
        dtimer = 1. # start timer
    else
        dtimer = 0. # stop timer
        timer  = 0. # reset timer
    end
    timer_temp = timer
    tap_temp   = tap

    if tap_side == "LV"
        if timer_temp >= tap_time && tap < tap_max
            tap_temp += 1
            timer = 0. #reset timer
        end
        U_LV = (1+(tap_pos+tap_temp)*tap_inc/100)*U_LV
    elseif tap_side == "HV"
        if timer_temp >= tap_time && tap > -tap_max
            tap_temp -= 1
            timer = 0. #reset timer
        end
        U_HV = (1+(tap_pos+tap_temp)*tap_inc/100)*U_HV
    else
        error("Can not interprete tap_side (HV/LV): $tap_side")
    end
    ü     = U_HV/U_LV #transformer ratio
    I_r   = S_r/(U_HV*sqrt(3)) #rated current

    #Calculatiing leakage reactance Xa and winding resistance Ra
    Ra = 0.
    Xa = 0.
    if XR_ratio == 0
        Ra = uk*U_HV/(sqrt(3)*I_r);
        Xa = 0.
    elseif XR_ratio == Inf
        Ra = 0.
        Xa = uk*U_HV/(sqrt(3)*I_r);
    else
        Ra = sqrt(((uk*U_HV/(I_r*sqrt(3)))^2)/(1+XR_ratio^2));
        Xa = XR_ratio*Ra;
    end
    Ya = 1/(0.5*(Ra+1im*Xa)) #it is assumed that losses are equally (0.5) distributed over both sides
    Ybs = Ya; #Ybs should be changed, if losses are not equally distributed

    #Calculating magnetising reactance Xm and core resistance Rfe from iron losess
    Zm  = U_HV/(sqrt(3)*i0/100.0*I_r) #no-load currents depends on complete magnetising impedance
    Rfe = U_HV*U_HV/(Pv0)
    Xm  = Inf
    if i0 != 0
        Xm = 1/(sqrt(1/Zm^2 - 1/Rfe^2))
    end
    Ym  = 1/Rfe + 1/(1im*Xm)

    Ybase = 1/(Ubase^2/Sbase) #this could be changed, if global base values are available
    Y     = 1/(Ya+Ybs+Ym)./Ybase
    Y = Y.*PiModel(Ya*Ybs, Ya*Ym, Ybs*Ym,ü,1)

    dtap    = (tap_temp-tap)*dtap
    # If current is flowing away from the source, it is negative at the source.
    voltage_vector = [source_voltage,destination_voltage]
    current_vector = Y * voltage_vector
end

export DynamicPowerTransformer
