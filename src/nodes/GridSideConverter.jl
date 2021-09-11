@doc """
```Julia
GridSideConverterOfDERs(mode, p_ref, q_ref, v_ref, idmax, iqmax, imax, Kp, Tp, Kq, Tq, Kv, Tv, Kgsc, Tgsc, δqv, v1_max, v1_min)
```

A node that represents the power control structure of the grid-side converter of distributed energy resources.
The modulation is based on the paper "Impact of Varying Shares of Distributed Energy Resources
on Voltage Stability in Electric Power Systems" by S. Liemann, L. Robitzky and C. Rehtanz.

The model has the following internal dynamic variables:
* ``x_st`` state of the Integrator for active power control
* ``y_st`` state of the Integrator for reactive power control (Mode 1 & 3)
* ``z_st`` state of the Integrator for reactive power control (Mode 2)
* ``id`` state of the first order time delay (P-T1) of id0
* ``iq`` state of the first order time delay (P-T1) of iq0

# Keyword Arguments
- `mode`:   mode 1 -> fix value q_ref is used for reactive power control
            mode 2 -> local voltage control determines the setpoint for reactive power
            mode 3 -> Q(v)-curve to control reactive power
- `p_ref`: setpoint for active power control
- `q_ref`: setpoint for reactive power control
- `v_ref`: setpoint for reactive power control, reference variable is voltage
- `idmax`: max id current
- `iqmax`: max iq current
- `imax`: max idq current, idq = id + 1im*iq
- `Kp`: parameter for active power PI-controller
- `Tp`: parameter for active power PI-controller
- `Kq`: parameter for reactive power PI-controller (Mode 1 & 3)
- `Tq`: parameter for reactive power PI-controller (Mode 1 & 3)
- `Kv`: parameter for reactive power PI-controller (Mode 2)
- `Tv`: parameter for reactive power PI-controller (Mode 2)
- `Kgsc`: parameter for first order time dealy (P-T1) of idq
- `Tgsc`: parameter for first order time dealy (P-T1) of idq
- `δqv`: width of the dead band of the Q(V)-curve
- `v1_max`: voltage at which max inductive power is induced by the Q(V)-curve
- `v1_min`: voltage at which max capacitive power is induced by the Q(V)-curve

"""

@DynamicNode GridSideConverter(mode, p_ref, q_ref, v_ref, idmax, iqmax, imax, Kp, Tp, Kq, Tq, Kv, Tv, Kgsc, Tgsc, δqv, v1_max, v1_min, q_max) begin
    MassMatrix(m_u = false, m_int =[true,true,true,true,true])
end begin
    @assert Tp > 0 "Time has to be >0"
    @assert Tq > 0 "Time has to be >0"
    @assert Tv > 0 "Time has to be >0"
    @assert Tgsc > 0 "Time has to be >0"
    @assert idmax > 0 "Upper current limit has to be >0"
    @assert iqmax > 0 "Upper current limit has to be >0"
    @assert idmax > 0 "Upper current limit has to be >0"
    @assert δqv >0 "Width of the dead band of the Q(V)-curve"
end [[x_st, dx_st], [y_st, dy_st], [z_st, dz_st], [id, did], [iq, diq]] begin

    p = real(u * conj(i))
    Δp = p_ref - p
    dx_st = Δp
    id0 = x_st/Tp + Kp*Δp

    if id0 > idmax
        id0 = idmax
        dx_st = 0                  #Anti-Reset-Windup
    elseif id0 < -idmax
        id0 = -idmax
        dx_st = 0                  #Anti-Reset-Windup
    end
    if mode == 1                    #a fix q-value is used to control the reactive power
        q = imag(u * conj(i))
        Δq = q - q_ref
        dy_st = Δq
        dz_st = 0
        iq0 = Kq*Δq + y_st/Tq
    elseif mode == 2                #local voltage control determines the setpoint for reactive power
        Δu = u - v_ref 
        dz_st = Δu
        dy_st = 0
        iq0 = Kv*Δu + z_st/Tv
    elseif mode == 3                #Q(v)-curve to control reactive power
        q = imag(u * conj(i))
        v = abs(u)

        #### Following implementation of a smooth Q(V)-curve is from "Integration of Voltage Dependent Power Injections of Distributed Generators into
        #### the Power Flow by using a Damped Newton Method"

        # Hilfsvariablen
        # v1_min = 0.95; v1_max = 1.05; δqv = 0.01; v = 1.02
        k = 600 # dieser Wert bestimmt wie nah die smooth function der sich der nicht-stetigen anpasst. Wert 600 ist aus dem Paper übernommen
        V1 = v1_min; V2 = 1.0 - δqv/2; V3 = 1.0 + δqv/2; V4 = v1_max;
        y1 = 1.0; y2 = 0.0; y3 = -1.0;
        γ1 = (y2-y1)/(V2 - V1) # Steigung erster Abschnitt
        γ2 = (y3-y2)/(V4 - V3) # Steigung zweiter Abschnitt

        # Q(V)-Funktion
        Q(v) = y1 + γ1/k*(log(1+exp(k*(v-V1))) - log(1+exp(k*(v-V2)))) + γ2/k*(log(1+exp(k*(v-V3))) - log(1+exp(k*(v-V4))))

        #Δq = q+Q(v)*q_max # Hier muss wahrscheinlich ein minus hin! Ausprobieren!
        Δq = q-Q(v)*q_max 
        dy_st = Δq
        dz_st = 0
        iq0 = Kq*Δq + y_st/Tq
    end

    if iq0 > iqmax
        iq0 = iqmax
        dy_st = 0                      #Anti-Reset-Windup
        dz_st = 0                      #Anti-Reset-Windup
    elseif iq0 < -iqmax
        iq0 = -iqmax
        dy_st = 0                      #Anti-Reset-Windup
        dz_st = 0                      #Anti-Reset-Windup
    end


    idq0 = id0+1im*iq0
    if (abs(id0) > imax || abs(iq0) > imax)
        id0 = imax; iq0 = 0;                   # Dieser Fall sollte eig nie auftreten. Zur Sicherheit hinzugefügt
    elseif abs(idq0) > imax                    #dynamic limiter
        if (abs(u) < 1.1) && (abs(u) > 0.9)    #normal case
            if iq0 <= 0
                iq0 = -sqrt(imax^2 - id0^2)
            else
                iq0 = sqrt(imax^2 - id0^2)
            end
        else                                    #FRT
            id0 = sqrt(imax^2-iq0^2)
        end
    end

    did = id0*Kgsc/Tgsc - id/Tgsc
    diq = iq0*Kgsc/Tgsc - iq/Tgsc
    idq = id + 1im*iq
    du = i - idq*exp(1im*angle(u))
end
export GridSideConverter
