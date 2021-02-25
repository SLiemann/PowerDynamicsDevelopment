using LinearAlgebra: I, Diagonal
import PowerDynamics: PiModel

begin
    @doc """
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
    struct DynamicPowerTransformer <: AbstractLine
            from
            to
            S_r
            U_HV
            U_LV
            uk
            XR_ratio
            i0
            Pv0
            Sbase
            Ubase
            tap_side
            tap_pos
            tap_inc
            tap_max
            v_ref
            v_dead
            tap_time
            DynamicPowerTransformer(; from,to,S_r,U_HV,U_LV,uk,XR_ratio,i0,Pv0,Sbase,Ubase,tap_side,tap_pos,tap_inc,tap_max,v_ref,v_dead,tap_time) =
                new(from,to,S_r,U_HV,U_LV,uk,XR_ratio,i0,Pv0,Sbase,Ubase,tap_side,tap_pos,tap_inc,tap_max,v_ref,v_dead,tap_time)
        end
    function construct_edge(par::DynamicPowerTransformer)
            from  = par.from
            to    = par.to
            S_r   = par.S_r
            U_HV  = par.U_HV
            U_LV  = par.U_LV
            uk    = par.uk
            XR_ratio = par.XR_ratio
            i0    = par.i0
            Pv0   = par.Pv0
            Sbase = par.Sbase
            Ubase = par.Ubase
            tap_side = par.tap_side
            tap_pos  = par.tap_pos
            tap_inc  = par.tap_inc
            tap_max  = par.tap_max
            v_ref    = par.v_ref
            v_dead   = par.v_dead
            tap_time = par.tap_time
        function rhs!(de, e, v_s, v_d, p, t)
            #dtimer,timer,tap,Y = DPT_Y(de[5],e[5],e[7],S_r,U_HV,U_LV,uk,XR_ratio,i0,Pv0,Sbase,Ubase,tap_side,tap_pos,tap_inc,tap_max,v_ref,v_dead,tap_time)
            abs_uHV = abs(v_s[1]+1im*v_s[2])
            abs_uLV = abs(v_d[1]+1im*v_d[2])

            dtimer_temp = de[5]
            timer_temp = e[5]
            tap_temp   = e[7]

            if abs_uLV < v_ref - v_dead || abs_uHV < v_ref + v_dead
                dtimer_temp = 1. # start timer
            else
                dtimer_temp = 0. # stop timer
                timer_temp  = 0. # reset timer
            end

            U_LV_temp = U_LV
            U_HV_temp = U_HV
            if tap_side == "LV"
                if timer_temp >= tap_time && tap_temp < tap_max
                    tap_temp += 1
                    timer_temp = 0. #reset timer
                end
                U_LV_temp = (1+(tap_pos+tap_temp)*tap_inc/100)*U_LV
            elseif tap_side == "HV"
                if timer_temp >= tap_time && tap_temp > -tap_max
                    tap_temp -= 1
                    timer = 0. #reset timer
                end
                U_HV_temp = (1+(tap_pos+tap_temp)*tap_inc/100)*U_HV
            else
                error("Can not interprete tap_side (HV/LV): $tap_side")
            end
            ü     = U_LV_temp/U_HV_temp #transformer ratio
            I_r   = S_r/(U_HV_temp*sqrt(3)) #rated current

            #Calculatiing leakage reactance Xa and winding resistance Ra
            Ra = 0.
            Xa = 0.
            if XR_ratio == 0
                Ra = uk*U_HV_temp/(sqrt(3)*I_r);
                Xa = 0.
            elseif XR_ratio == Inf
                Ra = 0.
                Xa = uk*U_HV_temp/(sqrt(3)*I_r);
            else
                Ra = sqrt(((uk*U_HV_temp/(I_r*sqrt(3)))^2)/(1+XR_ratio^2));
                Xa = XR_ratio*Ra;
            end
            Ya = 1/(0.5*(Ra+1im*Xa)) #it is assumed that losses are equally (0.5) distributed over both sides
            Ybs = Ya; #Ybs should be changed, if losses are not equally distributed

            #Calculating magnetising reactance Xm and core resistance Rfe from iron losess
            Zm  = U_HV_temp/(sqrt(3)*i0/100.0*I_r) #no-load currents depends on complete magnetising impedance
            Rfe = U_HV_temp*U_HV_temp/(Pv0)
            Xm  = Inf
            if i0 != 0
                Xm = 1/(sqrt(1/Zm^2 - 1/Rfe^2))
            end
            Ym  = 1/Rfe + 1/(1im*Xm)

            Ybase = 1/(Ubase^2/Sbase) #this could be changed, if global base values are available
            Y     = 1/(Ya+Ybs+Ym)./Ybase
            Y = Y.*PiModel(Ya*Ybs, Ya*Ym, Ybs*Ym,ü,1)

            # the different minus signs are due to the PowerDynamics sign convention for currents
            #e[5] = timer
            source_voltage = v_s[1] + v_s[2] * im
            destination_voltage = v_d[1] + v_d[2] * im
            voltage_vector = [source_voltage, destination_voltage]
            current_vector = Y * voltage_vector
            try
                de[1] = real(current_vector[1])
                de[2] = imag(current_vector[1])
                de[3] = real(current_vector[2])
                de[4] = imag(current_vector[2])
                de[5] = dtimer_temp
                de[6] = timer_temp
                de[7] = tap_temp
                return nothing
            catch e
                throw(e)
            end
        end
        return ODEEdge(f! = rhs!, dim=7, mass_matrix=Diagonal([0,0,0,0,1,0,0]), sym=Symbol[:id, :iq, :id_r, :iq_r,:dtime,:time,:tap])
    end
    symbolsof(::DynamicPowerTransformer) = begin
            [:id, :iq, :id_r, :iq_r,:dtime,:time,:tap]
    end
    dimension(::DynamicPowerTransformer) = begin
            7
    end
end

function DPT_Y(dtimer,timer,tap,S_r,U_HV,U_LV,uk,XR_ratio,i0,Pv0,Sbase,Ubase,tap_side,tap_pos,tap_inc,tap_max,v_ref,v_dead,tap_time)
        abs_uHV = abs(v_s[1]+1im*v_s[2])
        abs_uLV = abs(v_d[1]+1im*v_d[2])

        dtimer_temp = dtimer
        timer_temp = timer
        tap_temp   = tap

        if abs_uLV < v_ref - v_dead || abs_uHV < v_ref + v_dead
            dtimer_temp = 1. # start timer
        else
            dtimer_temp = 0. # stop timer
            timer_temp  = 0. # reset timer
        end

        if tap_side == "LV"
            if timer_temp >= tap_time && tap_temp < tap_max
                tap_temp += 1
                timer_temp = 0. #reset timer
            end
            U_LV = (1+(tap_pos+tap_temp)*tap_inc/100)*U_LV
        elseif tap_side == "HV"
            if timer_temp >= tap_time && tap_temp > -tap_max
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

        return dtimer,timer,tap,Y
end

export DynamicPowerTransformer
