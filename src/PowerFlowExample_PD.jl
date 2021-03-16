using PowerDynamics: SlackAlgebraic, FourthOrderEq, VoltageDependentLoad, PiModelLine, PVAlgebraic, Transformer
using PowerDynamics
using OrderedCollections: OrderedDict
using Plots

begin
    include("operationpoint/PowerFlow_PD.jl");
    Ubase = 380e3
    Sbase = 100e6

    buses=OrderedDict(
        "bus1"=> SlackAlgebraic(U=0.98),
        "bus2"=> VoltageDependentLoad(P=-0.3, Q=0.3, U=1.0, A=0.0, B=0.0),
        "bus3"=> FourthOrderEq(T_d_dash=6.1, D=2, X_d=1.05, X_q=0.98, Ω=50, X_d_dash=0.185, T_q_dash=0.4, X_q_dash=0.36, P=0.5, H=6.54, E_f= 1),
        "bus4"=> VoltageDependentLoad(P=-0.5, Q=-0.5, U=1.0, A=0.5, B=0.3),
        "bus5"=> PVAlgebraic(P=-0.0, V = 1.01))

    branches=OrderedDict(
        "branch1"=> PiModelLine(from= "bus1", to = "bus2",y=1.0/(0.05+1im*0.15), y_shunt_km=0., y_shunt_mk=0.),
        "branch2"=> PiModelLine(from= "bus2", to = "bus3",y=1.0/(0.05+1im*0.15), y_shunt_km=0., y_shunt_mk=0.),
        "branch3"=> Transformer(from= "bus3", to = "bus4",y=1.0/(0.05+1im*0.15), t_ratio = 1.0/0.93),
        "branch4"=> PiModelLine(from= "bus4", to = "bus5",y=1.0/((0.05+1im*0.15)*(Ubase/110e3)^2), y_shunt_km=0., y_shunt_mk=0.))
        #note that branch 4 has to be transformed to higher voltage level!
    powergrid = PowerGrid(buses, branches)

    Qmax   = [Inf, Inf, 0.5, Inf, Inf]
    Qmin   = -Qmax
    U,δ1,ic = PowerFlowClassic(powergrid,iwamoto = true, Qmax = Qmax, Qmin = Qmin, Qlimit_iter_check = 2)
end

Uc = U.*exp.(1im*δ1/180*pi)
Ykk = NodalAdmittanceMatrice(powergrid)
S  = Uc.*(conj.(Ykk)*conj.(Uc))
