using PowerDynamics: SlackAlgebraic, FourthOrderEq, VoltageDependentLoad, PiModelLine, StaticLine, Transformer, PowerGrid, write_powergrid, Json, Inc, find_operationpoint, ChangeInitialConditions, LineFailure, PowerPerturbation, simulate
using PowerDynamics: symbolsof
using PowerDynamics: VoltageDependentLoadstep, StaticPowerTransformer
using PowerDynamics
using PowerDynamics: initial_guess
using PowerDynamics: guess
import PowerDynamics: PiModel
using OrderedCollections: OrderedDict
using OrdinaryDiffEq: ODEProblem, Rodas4
using Plots

begin
    include("plotting.jl")
    include("PowerFlow.jl")

    buses=OrderedDict(
        "bus1"=> SlackAlgebraic(U=1),
        "bus2"=> VoltageDependentLoad(P=0.3, Q=-0.3, U=1.0, A=0., B=0.),
        "bus3"=> FourthOrderEq(T_d_dash=6.1, D=2, X_d=1.05, X_q=0.98, Ω=50, X_d_dash=0.185, T_q_dash=0.4, X_q_dash=0.36, P=0.5, H=6.54, E_f= 1),
        "bus4"=> VoltageDependentLoad(P=0.5, Q=0.5, U=1.0, A=0., B=0.))

    branches=OrderedDict(
        "branch1"=> PiModelLine(from= "bus1", to = "bus2",y=1/(0.05+1im*0.15), y_shunt_km=0., y_shunt_mk=0.),
        "branch2"=> PiModelLine(from= "bus2", to = "bus3",y=1/(0.05+1im*0.15), y_shunt_km=0., y_shunt_mk=0.),
        "branch3"=> StaticPowerTransformer(from="bus3",to="bus4",S_r=100e6,U_HV=380e3,U_LV=110e3,uk=0.1581138,XR_ratio=5,i0=6.35,Pv0=100e3,Sbase=100e6,Ubase=380e3,tap_side = "LV",tap_pos = 5,tap_inc = 1))
        #"branch3"=> Transformer(from= "bus3", to = "bus4",y=1/(0.05+1im*0.15),t_ratio = 1.05))

    powergrid = PowerGrid(buses, branches)
    Unodes = [380e3,380e3,380e3,110e3]
    U,δ1,ic = PowerFlowClassic(powergrid,Unodes,iwamoto = false)
end
#TO-Do: Woher kommt der Unterschied in der LF-Rechnung zwischen Julia und PF?

Uc = U.*exp.(1im*δ1/180*pi)
Ykk = NodalAdmittanceMatrice(powergrid,Unodes,380e3)
S  = Uc.*(conj.(Ykk)*conj.(Uc))
