import Base: @__doc__
import PowerDynamics: AbstractLine, PiModel, PiModelLine, construct_edge, dimension
using NetworkDynamics: ODEEdge
using LinearAlgebra: Diagonal
using IfElse

begin
    @__doc__ struct OLTC <: AbstractLine
            from
            to
            Sbase
            Srated
            uk
            XR_ratio
            i0
            Pv0
            tap_side
            tap_pos
            tap_inc
            tap_max
            tap_min
            p_ind
            OLTC(; from, to, Sbase, Srated, uk, XR_ratio, i0, Pv0, tap_side, tap_pos, tap_inc, tap_max, tap_min, p_ind) = new(from, to, Sbase, Srated, uk, XR_ratio, i0, Pv0, tap_side, tap_pos, tap_inc, tap_max, tap_min, p_ind)
        end
    function construct_edge(par::OLTC)
        from = par.from
        to = par.to
        Sbase = par.Sbase
        Srated = par.Srated
        uk = par.uk
        XR_ratio = par.XR_ratio
        i0 = par.i0
        Pv0 = par.Pv0
        tap_side = par.tap_side
        tap_pos = par.tap_pos
        tap_inc = par.tap_inc
        tap_max = par.tap_max
        tap_min = par.tap_min
        p_ind = par.p_ind
        function rhs!(de,e, v_s, v_d, p, t)
            source_voltage = v_s[1] + v_s[2] * im
            destination_voltage = v_d[1] + v_d[2] * im
            U_LV = abs(destination_voltage)
            ON_OFF = IfElse.ifelse(0.99<U_LV <1.01)

            Δtap = p[p_ind[1]]
            Δtap = IfElse.ifelse(tap_pos + Δtap >= tap_max, tap_max - tap_pos, IfElse.ifelse(tap_pos + Δtap <= tap_min, tap_min - tap_pos, p[p_ind]))
            üLV = IfElse.ifelse(tap_side == "LV", 1.0 / (1.0 + ((tap_pos + Δtap) * tap_inc) / 100.0), 1.0)
            üHV = IfElse.ifelse(tap_side == "HV", 1.0 / (1.0 + ((tap_pos + Δtap) * tap_inc) / 100.0), 1.0)
            Ra = IfElse.ifelse(XR_ratio == 0.0, uk, IfElse.ifelse(XR_ratio == Inf, 0.0, uk / sqrt(1 + XR_ratio ^ 2)))
            Xa = IfElse.ifelse(XR_ratio == 0.0, 0.0, IfElse.ifelse(XR_ratio == Inf, uk, XR_ratio * Ra))
            Ya = 1.0 / (0.5 * (Ra + (1im) * Xa))
            Ybs = Ya
            Zm = 1.0 / (i0 / 100.0)
            Rfe = Srated / Pv0
            Xm = IfElse.ifelse(i0 != 0.0, 1.0 / sqrt(1.0 / Zm ^ 2 - 1.0 / Rfe ^ 2), Inf)
            Ym = 1.0 / Rfe + 1.0 / ((1im) * Xm)
            Ybase = Sbase / Srated
            Y = (1.0 / (Ya + Ybs + Ym)) ./ Ybase
            Y = Y .* PiModel(Ya * Ybs, Ya * Ym, Ybs * Ym, üHV, üLV)
            voltage_vector = [source_voltage, destination_voltage]
            current_vector = Y * voltage_vector
            de[1] = real(current_vector[1])
            de[2] = -(real((1im) * current_vector[1]))
            de[3] = real(current_vector[2])
            de[4] = -(real((1im) * current_vector[2]))
            de[5] = ON_OFF
        end
        return ODEEdge(f! = rhs!, dim = 5, mass_matrix= Diagonal([0,0,0,0,1]),sym=Symbol[:id, :iq, :id_r, :iq_r,:timer])
    end
end
