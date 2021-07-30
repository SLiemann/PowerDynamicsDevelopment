@doc """
```Julia
GridSideConverterOfDERs(mode, p_ref, q_ref, v_ref, idmax, iqmax, imax, Kp, Tp, Kq, Tq, Kv, Tv, Kgsc, Tgsc, delta_qv, v1_max, v1_min)
```

A node that represents the power control structure of the grid-side converter of distributed energy resources.
The modulation is based on the paper "Impact of Varying Shares of Distributed Energy Resources
on Voltage Stability in Electric Power Systems" by S. Liemann, L. Robitzky and C. Rehtanz.

The model has the following internal dynamic variables:
* ``x`` state of the Integrator for active power control
* ``y`` state of the Integrator for reactive power control (Mode 1 & 3)
* ``z`` state of the Integrator for reactive power control (Mode 2)
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
- `delta_qv`: width of the dead band of the Q(V)-curve
- `v1_max`: voltage at which max inductive power is induced by the Q(V)-curve
- `v1_min`: voltage at which max capacitive power is induced by the Q(V)-curve

"""

@DynamicNode GridSideConverter(mode, p_ref, q_ref, v_ref, idmax, iqmax, imax, Kp, Tp, Kq, Tq, Kv, Tv, Kgsc, Tgsc, delta_qv, v1_max, v1_min, q_max) begin
    MassMatrix(m_u = false, m_int =[true,true,true,true,true])
end begin
    @assert Tp > 0 "Time has to be >0"
    @assert Tq > 0 "Time has to be >0"
    @assert Tv > 0 "Time has to be >0"
    @assert Tgsc > 0 "Time has to be >0"
    @assert idmax > 0 "Upper current limit has to be >0"
    @assert iqmax > 0 "Upper current limit has to be >0"
    @assert idmax > 0 "Upper current limit has to be >0"
    @assert delta_qv >0 "Width of the dead band of the Q(V)-curve"
end [[dx, x], [dy, y], [dz, z], [did, id], [diq, iq]] begin

    p = real(u * conj(i))
    Δp = p_ref - p
    dx = ΔP
    id0 = x/Tp + Kp*ΔP

    if id0 > idmax
        id0 = idmax
        dx = 0                  #Anti-Reset-Windup
    elseif id0 < -idmax
        id0 = -idmax
        dx = 0                  #Anti-Reset-Windup
    end
    if mode == 1                    #a fix q-value is used to control the reactive power
        q = imag(u * conj(i))
        Δq = q - q_ref
        dy = Δq
        dz = 0
        iq0 = Kq*Δq + y/Tq
    elseif mode == 2                #local voltage control determines the setpoint for reactive power
        Δu = u - v_ref 
        dz = Δu
        dy = 0
        iq0 = Kv*Δu + z/Tv
    elseif mode == 3                #Q(v)-curve to control reactive power
        q = imag(u * conj(i))
        v = abs(u)
        if v < v1_min                 #Q(v)-characteristic
            Q = 1
        elseif v < 1-delta_qv/2
            m = -1/(1-delta_qv/2-v1_min)
            b = 1-m*v1_min
            Q = m*v+b
        elseif v < 1+delta_qv/2
            Q = 0
        elseif v < v1_max
            m = -1/(v1_max-1-delta_qv/2)
            b = -1-m*v1_max
            Q = m*v+b
        else
            Q = -1
        end
        Δq = q+Q*q_max
        dy = Δq
        dz = 0
        iq0 = Kq*Δq + y/Tq
    end

    if iq0 > iqmax
        iq0 = iqmax
        dy = 0                      #Anti-Reset-Windup
        dz = 0                      #Anti-Reset-Windup
    elseif iq0 < -iqmax
        iq0 = -iqmax
        dy = 0                      #Anti-Reset-Windup
        dz = 0                      #Anti-Reset-Windup
    end


    idq0 = id0+1im*iq0
    if abs(idq0) > imax                        #dynamic limiter
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
    du = i - i_dq*exp(1im*angle(u))
end
export GridSideConverter
