using PowerDynamics: SlackAlgebraic, FourthOrderEq, VoltageDependentLoad, PiModelLine, Transformer
using PowerDynamics
using OrderedCollections: OrderedDict
using Plots

begin
    include("operationpoint/PowerFlow_PD.jl");

    buses=OrderedDict(
        "bus1"=> SlackAlgebraic(U=1),
        "bus2"=> VoltageDependentLoad(P=-0.3, Q=0.3, U=1.0, A=0.0, B=0.0),
        "bus3"=> FourthOrderEq(T_d_dash=6.1, D=2, X_d=1.05, X_q=0.98, Ω=50, X_d_dash=0.185, T_q_dash=0.4, X_q_dash=0.36, P=0.5, H=6.54, E_f= 1),
        "bus4"=> VoltageDependentLoad(P=-0.5, Q=-0.5, U=1.0, A=0.0, B=0.0))

    branches=OrderedDict(
        "branch1"=> PiModelLine(from= "bus1", to = "bus2",y=1.0/(0.05+1im*0.15), y_shunt_km=0., y_shunt_mk=0.),
        "branch2"=> PiModelLine(from= "bus2", to = "bus3",y=1.0/(0.05+1im*0.15), y_shunt_km=0., y_shunt_mk=0.),
        "branch3"=> Transformer(from= "bus3", to = "bus4",y=1.0/(0.05+1im*0.15), t_ratio = 380e3/110e3))

    powergrid = PowerGrid(buses, branches)
    Unodes = [380e3,380e3,380e3,110e3]
    Ubase  = 380e3
    Qmax   = [Inf, Inf, 0.5, Inf]
    Qmin   = -Qmax
    U,δ1,ic = PowerFlowClassic(powergrid,Unodes,Ubase,iwamoto = true, Qmax = Qmax, Qmin = Qmin)
end

Uc = U.*exp.(1im*δ1/180*pi)
Ykk = NodalAdmittanceMatrice(powergrid,Unodes,380e3)
S  = Uc.*(conj.(Ykk)*conj.(Uc))
