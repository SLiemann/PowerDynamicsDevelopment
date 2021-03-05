using PowerDynamics: SlackAlgebraic, FourthOrderEq, VoltageDependentLoad, PiModelLine, StaticLine, Transformer, PowerGrid#, write_powergrid, Json, Inc, find_operationpoint, ChangeInitialConditions, LineFailure, PowerPerturbation, simulate
using PowerDynamics: symbolsof, initial_guess, guess, RLLine
using PowerDynamics: StaticPowerTransformer, DynamicPowerTransformer, SixOrderMarconatoMachine
using PowerDynamics
using PowerDynamics: rhs, State
using OrderedCollections: OrderedDict
using Plots
import PowerDynamics: PiModel
using DifferentialEquations
using CSV #read PF DataFrames
using DataFrames #for CSV

include("operationpoint/PowerFlow.jl")

begin

    buses=OrderedDict(
        "bus1"=> SlackAlgebraic(U=1),
        "bus2"=> VoltageDependentLoad(P=-0.3, Q=0.3, U=1.0, A=0., B=0.,Y_n = complex(0.0)),
        #"bus3"=> SixOrderMarconatoMachine(H = 5, P=0.5, D=0., Ω=50*2*pi, E_f=1.4, R_a = 0.1,T_ds=1.136,T_qs=0.8571,T_dss=0.04,T_qss=0.06666,X_d=1.1,X_q=0.7,X_ds=0.25,X_qs=0.25,X_dss=0.2,X_qss=0.2,T_AA=0.))
        "bus3"=> FourthOrderEq(H=5, P=0.5, D=0., Ω=50, E_f=1.0, T_d_dash=1.136 ,T_q_dash=0.8571 ,X_q_dash =0.25 ,X_d_dash=0.25,X_d=1.1, X_q=0.7))
    branches=OrderedDict(
        "branch1"=> PiModelLine(from= "bus1", to = "bus2",y=1.0/(0.05+1im*0.15), y_shunt_km=0., y_shunt_mk=0.),
        "branch2"=> PiModelLine(from= "bus2", to = "bus3",y=1.0/(0.05+1im*0.15), y_shunt_km=0., y_shunt_mk=0.))
        #"branch3"=> DynamicPowerTransformer(from="bus3",to="bus4",S_r=100e6,U_HV=380e3,U_LV=110e3,uk=0.1581138,XR_ratio=5,i0=6.35,Pv0=100e3,Sbase=100e6,Ubase=380e3,tap_side = "LV",tap_pos = 5,tap_inc = 1,tap_max=10,v_ref=1.,v_dead=0.05,tap_time=5.))

    pg = PowerGrid(buses, branches)
    Unodes = [380e3,380e3,380e3]
    Qmax   = [Inf, Inf, Inf]
    Qmin   = -Qmax
    U,δ1,ic = PowerFlowClassic(pg,Unodes,iwamoto = true, Qmax = Qmax, Qmin = Qmin)
end
ic = initial_guess(pg)
ic0 = find_valid_initial_condition(pg,ic)
begin
    include("operationpoint/InitializeInternals.jl")
    Uc = U.*exp.(1im*δ1/180*pi)
    Ykk = NodalAdmittanceMatrice(pg,Unodes,380e3)
    I_c = Ykk*Uc
    PG, ic0 = InitializeInternalDynamics(pg,I_c,ic)
    #ODEProb = ODEProblem{true}(rhs(PG),ic0,[0,15.])
    Zbase = (380e3^2)/(100e6)
    SS = NodeShortCircuit(node="bus2",Y = 1/(1im*1000/Zbase),tspan_fault=(1.0, 1.15))
    PT = PowerPerturbation(node="bus2", fault_power = -0.3 ,tspan_fault=(1.0, 5.0))
    PGsol = simulate(PT,PG,ic0,(0.,5.))
    #PG_state = State(PG,ic0)
    #PGsol = solve(PG,PG_state,[0.,1.])
end
plot(PGsol,collect(keys(PG.nodes)), :v,size = (1000, 500),legend = (0.3, 0.3))
begin
    plot(PGsol,collect(keys(PG.nodes)), :v,size = (1000, 500),legend = (0.5, 0.5))
    test = DataFrame(CSV.File("C:\\Users\\liemann\\Desktop\\PF_test.csv"; header=false, delim=';', type=Float64))
    plot!(test.Column1,test.Column2,label = "PF-bus1")
    plot!(test.Column1,test.Column3,label = "PF-bus2")
    plot!(test.Column1,test.Column4,label = "PF-bus3")
end
begin
    plot(PGsol,["bus3"], :θ)
    test2 = DataFrame(CSV.File("C:\\Users\\liemann\\Desktop\\PF_pol.csv"; header=false, delim=';', type=Float64))
    plot!(test2.Column1,(test2.Column2))
    xlims!((0.9,10.))
end
begin
    plot(PGsol,["bus3"], :ω)
    test3 = DataFrame(CSV.File("C:\\Users\\liemann\\Desktop\\PF_dreh.csv"; header=false, delim=';', type=Float64))
    plot!(test3.Column1,(test3.Column2.-1))
end
#ylims!((1.01,1.015))
#ylims!((0.75,1.015))
xlims!((0.9,1.5))
plot(PGsol,["bus3"], :θ)
plot(PGsol,["bus3"], :ω)
plot(PGsol,["bus3"], :p)
plot(PGsol,["bus3"], :q)
plot(PGsol,["bus3"], :e_dss)
plot(PGsol,["bus3"], :e_qss)
