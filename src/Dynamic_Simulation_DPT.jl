using PowerDynamics: SlackAlgebraic, FourthOrderEq, VoltageDependentLoad, PiModelLine, StaticLine, Transformer, PowerGrid#, write_powergrid, Json, Inc, find_operationpoint, ChangeInitialConditions, LineFailure, PowerPerturbation, simulate
using PowerDynamics: symbolsof, initial_guess, guess, RLLine
using PowerDynamics: StaticPowerTransformer, DynamicPowerTransformer, SixOrderMarcanatoMachine
using PowerDynamics
using PowerDynamics: rhs
using OrderedCollections: OrderedDict
using Plots
import PowerDynamics: PiModel
using DifferentialEquations

include("operationpoint/PowerFlow.jl")

begin

    buses=OrderedDict(
        "bus1"=> SlackAlgebraic(U=1),
        "bus2"=> VoltageDependentLoad(P=0.3, Q=-0.3, U=1.0, A=0., B=0.),
        #"bus3"=> FourthOrderEq(T_d_dash=6.1, D=2, X_d=1.05, X_q=0.98, Ω=50, X_d_dash=0.185, T_q_dash=0.4, X_q_dash=0.36, P=0.5, H=6.54, E_f= 1.625))
        "bus3"=> SixOrderMarcanatoMachine(H = 5, P=0.5, D=5, Ω=50, E_f=1.4, R_a = 0.0,T_ds=1.136,T_qs=0.8571,T_dss=0.04,T_qss=0.06666,X_d=1.1,X_q=0.7,X_ds=0.25,X_qs=0.3,X_dss=0.2,X_qss=0.2,T_AA=0.))
        #"bus4"=> VoltageDependentLoad(P=0.5, Q=0.5, U=1.0, A=0.5, B=0.2))

    #branches=OrderedDict(
    #    "branch1"=> PiModelLine(from= "bus1", to = "bus2",y=1/(0.05+1im*0.15), y_shunt_km=0., y_shunt_mk=0.),
    #    "branch2"=> PiModelLine(from= "bus2", to = "bus3",y=1/(0.05+1im*0.15), y_shunt_km=0., y_shunt_mk=0.),
    #    "branch3"=> StaticPowerTransformer(from="bus3",to="bus4",S_r=100e6,U_HV=380e3,U_LV=110e3,uk=0.1581138,XR_ratio=5,i0=6.35,Pv0=100e3,Sbase=100e6,Ubase=380e3,tap_side = "LV",tap_pos = 5,tap_inc = 1))

    branches=OrderedDict(
        "branch1"=> PiModelLine(from= "bus1", to = "bus2",y=1/(0.05+1im*0.15), y_shunt_km=0., y_shunt_mk=0.),
        "branch2"=> PiModelLine(from= "bus2", to = "bus3",y=1/(0.05+1im*0.15), y_shunt_km=0., y_shunt_mk=0.))
        #"branch3"=> DynamicPowerTransformer(from="bus3",to="bus4",S_r=100e6,U_HV=380e3,U_LV=110e3,uk=0.1581138,XR_ratio=5,i0=6.35,Pv0=100e3,Sbase=100e6,Ubase=380e3,tap_side = "LV",tap_pos = 5,tap_inc = 1,tap_max=10,v_ref=1.,v_dead=0.05,tap_time=5.))

    #branches=OrderedDict(
    #    "branch1"=> RLLine(from= "bus1", to = "bus2",R= 0.05, L = 0.15, ω0 = 1.),
    #    "branch2"=> RLLine(from= "bus2", to = "bus3",R= 0.05, L = 0.15, ω0 = 1.),
    #    "branch3"=> DynamicPowerTransformer(from="bus3",to="bus4",S_r=100e6,U_HV=380e3,U_LV=110e3,uk=0.1581138,XR_ratio=5,i0=6.35,Pv0=100e3,Sbase=100e6,Ubase=380e3,tap_side = "LV",tap_pos = 5,tap_inc = 1,tap_max=10,v_ref=1.,v_dead=0.05,tap_time=5.))

    pg = PowerGrid(buses, branches)
    Unodes = [380e3,380e3,380e3]
    Qmax   = [Inf, Inf, Inf]
    Qmin   = -Qmax
    U,δ1,ic = PowerFlowClassic(pg,Unodes,iwamoto = true, Qmax = Qmax, Qmin = Qmin)
    #pg = powergrid
end
#ic[7] =0.5
#ic2 = find_valid_initial_condition(powergrid,ic)
#x0_rand = zeros(25,1)
#x0_rand[1:10] .= ic

# Pm richtig initialisieren? 
begin
    include("operationpoint/InitializeInternals.jl")
    Uc = U.*exp.(1im*δ1/180*pi)
    Ykk = NodalAdmittanceMatrice(pg,Unodes,380e3)
    I_c = Ykk*Uc
    PG, ic0 = InitializeInternalDynamics(pg,I_c,ic)
    ODEProb = ODEProblem{true}(rhs(PG),ic0,[0,100.])
    test = solve(ODEProb,Rodas4())
    plot(test)
end

#=
begin
    ic2 = find_valid_initial_condition(pg,ic)
    ODEProb = ODEProblem{true}(rhs(pg),ic2,[0,100.])
    test = solve(ODEProb,Rodas4())
    plot(test)
end
=#
#Δic = ic .- pg_st.vec
#S  = Uc.*(conj.(Ykk)*conj.(Uc))
