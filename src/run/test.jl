using PowerDynamics
using OrderedCollections: OrderedDict
#using Plots
import PowerDynamics: PiModel
using DifferentialEquations
using CSV #read PF DataFrames
using DataFrames #for CSV
using Traceur
using Distributed
using Makie
@everywhere using IfElse


include("C:/Users/liemann/github/PowerDynamicsDevelopment/src/utility/utility_functions.jl")
include("C:/Users/liemann/github/PowerDynamicsDevelopment/src/operationpoint/InitializeInternals.jl")
include("C:/Users/liemann/github/PowerDynamicsDevelopment/src/operationpoint/simulate_fault.jl")

begin
    include("C:/Users/liemann/github/PowerDynamicsDevelopment/src/operationpoint/PowerFlow.jl")
    Ubase = 380e3
    Sbase = 100e6
    Zbase = (Ubase^2)/(Sbase)
    buses=OrderedDict(
        "bus1"=> SlackAlgebraic(U=0.98),
        "bus2"=> VoltageDependentLoad(P=-0.3, Q=0.3, U=1.0, A=0., B=0.,Y_n = complex(0.0)),
        "bus3"=> VoltageDependentLoad(P= 0.0, Q=0.0, U=1.0, A=0.0, B=0.0,Y_n = complex(0.0)),
        "bus4"=> VoltageDependentLoad(P=-0.5, Q=-0.5, U=1.0, A=0.5, B=0.3,Y_n = complex(0.0)),
        "bus5"=> SixOrderMarconatoMachineAVROEL(Sbase=100e6,Srated=100e6,H = 5, P=0.5, D=0., Ω=50, R_a = 0.1,
                                             T_ds=1.136,T_qs=0.8571,T_dss=0.04,T_qss=0.06666,
                                             X_d=1.1,X_q=0.7,X_ds=0.25,X_qs=0.25,X_dss=0.2,
                                             X_qss=0.2,T_AA=0.,V0 = 1.0, Ifdlim = 1.8991,
                                             L1 = -11.0, G1 = 70.0, Ta = 10.0, Tb = 20.0,
                                             G2 = 10.0, L2 = 4.0))

    branches=OrderedDict(
        "branch1"=> PiModelLine(from= "bus1", to = "bus2",y=1.0/(0.05+1im*0.15), y_shunt_km=0., y_shunt_mk=0.),
        "branch2"=> PiModelLine(from= "bus2", to = "bus3",y=1.0/(0.05+1im*0.15), y_shunt_km=0., y_shunt_mk=0.),
        "branch3"=> StaticPowerTransformer(from="bus3",to="bus4",Sbase=Sbase,Srated=100e6,uk=0.1581138,XR_ratio=3,
                                           i0=6.35,Pv0=300e3,tap_side = "LV",tap_pos = 3,tap_inc = 1.0),
        "branch4"=> StaticPowerTransformer(from="bus4",to="bus5",Sbase=Sbase,Srated=200e6,uk=0.10,XR_ratio=7,
                                           i0=0.0,Pv0=100e3,tap_side = "LV",tap_pos = 0,tap_inc = 1.0))
    pg = PowerGrid(buses, branches)
    Qmax   = [Inf, Inf, Inf,Inf, Inf]
    Qmin   = -Qmax
    U,δ1,ic = PowerFlowClassic(pg,iwamoto = true, Qmax = Qmax, Qmin = Qmin, Qlimit_iter_check = 2)
end

begin
    Ykk = NodalAdmittanceMatrice(pg)
    Uc = U.*exp.(1im*δ1/180*pi)
    I_c = Ykk*Uc
    S = conj(Ykk*Uc).*Uc
    dt = 1e-3
    tspan = (0.0,10.0)
    pg, ic0 = InitializeInternalDynamics(pg,I_c,ic)
    SS = NodeShortCircuit(node="bus2",Y = 1/(1im*250/Zbase),tspan_fault=(1.0, 1.15))
    PGsol = my_simulate(SS,pg,ic0,(0.,30.0))
end

begin
    plot(PGsol,collect(keys(pg.nodes)), :v,size = (2000, 1000),legend = (0.5, 0.5))
    test = DataFrame(CSV.File("C:\\Users\\liemann\\Desktop\\PF_test.csv"; header=false, delim=';', type=Float64))
    plot!(test.Column1,test.Column2,label = "PF-bus1")
    plot!(test.Column1,test.Column3,label = "PF-bus2")
    plot!(test.Column1,test.Column4,label = "PF-bus3")
    plot!(test.Column1,test.Column5,label = "PF-bus4")
    plot!(test.Column1,test.Column6,label = "PF-bus5")
end

xlims!((0.9,1.2))
