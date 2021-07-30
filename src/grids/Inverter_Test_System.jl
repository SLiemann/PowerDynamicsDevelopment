using PowerDynamics
using OrderedCollections: OrderedDict
using Distributed
@everywhere using IfElse
using ModelingToolkit

#include("C:/Users/Micha/Documents/Master-Studium/HiWi/PowerDynamicsDevelopment/PowerDynamicsDevelopment/src/nodes/GridSideConverter.jl")
#include("C:/Users/Micha/Documents/Master-Studium/HiWi/PowerDynamicsDevelopment/PowerDynamicsDevelopment/src/lines/StaticPowerTransformer.jl")

function Inverter_Test_System()
    Sbase = 100e6
    Ubase = 380e3
    Ibase = Sbase/Ubase/sqrt(3)
    Zbase = Ubase^2/Sbase

    buses=OrderedDict(
        "bus1" => SlackAlgebraic(U=1.0),
        "bus2" => VoltageDependentLoad(P=0.0, Q=0.0, U=1.0, A=0., B=0.,Y_n = complex(0.0)),
        "bus3" => VoltageDependentLoad(P=-10.0, Q = -3.286841, U=1.0, A=0.0, B=1.0,Y_n = complex(0.0)),
        "bus4" => GridSideConverter(mode=3, p_ref=12.0, q_ref=3.0, v_ref=1.0,
                                            idmax=1.5, iqmax=1.5, imax=1.5,
                                            Kp=1.0, Tp=0.1, Kq=1.0, Tq=0.1,
                                            Kv=1.0, Tv=0.1, Kgsc=1.0, Tgsc=1.0,
                                            delta_qv=0.01, v1_max=1.05, v1_min=0.95, q_max = 5),
        "busv" => VoltageDependentLoad(P=0.0, Q=0.0, U=1.0, A=0., B=0.,Y_n = complex(0.0)))

    Z_EHV_Line = (9.6 + 1im*64)/Zbase
    B_half     = 1im*1498.54*1e-6 / 2.0 *Zbase #already an admittance: 1498.54 = 2*pi*50*4.77001*10e-6
    branches=OrderedDict(
        "Line_1-2"=> PiModelLine(from= "bus1", to = "bus2",y=1.0/Z_EHV_Line, y_shunt_km=B_half, y_shunt_mk=B_half),
        "Line_1-v"=> PiModelLine(from= "bus1", to = "busv",y=1.0/(Z_EHV_Line/2.0), y_shunt_km=B_half, y_shunt_mk=0.0),
        "Line_v-2"=> PiModelLine(from= "bus2", to = "busv",y=1.0/(Z_EHV_Line/2.0), y_shunt_km=0.0, y_shunt_mk=B_half),
        "branch3"=> PiModelLine(from= "bus2", to = "busv",y=1.0/(Z_EHV_Line/2.0), y_shunt_km=0.0, y_shunt_mk=B_half),
        "branch4"=> PiModelLine(from= "bus2", to = "busv",y=1.0/(Z_EHV_Line/2.0), y_shunt_km=0.0, y_shunt_mk=B_half))
        return PowerGrid(buses, branches)
end

function Inverter_Test_System2()
    buses=OrderedDict(
    #"bus1"=> FourthOrderEq(T_d_dash=7.4, D=2, X_d=0.8979, X_q=0.646, Ω=50, X_d_dash=0.2995, T_q_dash=0.1, X_q_dash=0.646, P=2.32, H=5.148, E_f=1),
    "bus1"=> GridSideConverter(mode=3, p_ref=2.32, q_ref=3.0, v_ref=1.0,
        idmax=1.5, iqmax=1.5, imax=1.5,
        Kp=1.0, Tp=0.1, Kq=1.0, Tq=0.1,
        Kv=1.0, Tv=0.1, Kgsc=1.0, Tgsc=1.0,
        delta_qv=0.01, v1_max=1.05, v1_min=0.95, q_max = 0.1),
    "bus2"=> SlackAlgebraic(U=1),
    "bus3"=> FourthOrderEq(T_d_dash=6.1, D=2, X_d=1.05, X_q=0.98, Ω=50, X_d_dash=0.185, T_q_dash=0.4, X_q_dash=0.36, P=-0.942, H=6.54, E_f= 1),
    "bus4"=> PQAlgebraic(P=-0.478, Q=-0.0),
    "bus5"=> PQAlgebraic(P=-0.076, Q=-0.016),
    "bus6"=> FourthOrderEq(T_d_dash=4.75, D=2, X_d=1.25, X_q=1.22, Ω=50, X_d_dash=0.232, T_q_dash=1.6, X_q_dash=0.715, P=-0.122, H=5.06, E_f= 1),
    "bus7"=> PQAlgebraic(P=-0.0, Q=-0.0),
    "bus8"=> FourthOrderEq(T_d_dash=4.75, D=2, X_d=1.25, X_q=1.22, Ω=50, X_d_dash=0.232, T_q_dash=1.6, X_q_dash=0.715, P=0.0, H=5.06, E_f= 1),
    "bus9"=> PQAlgebraic(P=-0.295, Q=-0.166),
    "bus10"=> PQAlgebraic(P=-0.09, Q=-0.058),
    "bus11"=> PQAlgebraic(P=-0.035, Q=-0.018),
    "bus12"=> PQAlgebraic(P=-0.061, Q=-0.016),
    "bus13"=> PQAlgebraic(P=-0.135, Q=-0.058),
    "bus14"=> PQAlgebraic(P=-0.149, Q=-0.05));

branches=OrderedDict(
    "branch1"=> PiModelLine(from= "bus1", to = "bus2",y=4.999131600798035-1im*15.263086523179553, y_shunt_km=0.0528/2, y_shunt_mk=0.0528/2),
    "branch2"=> PiModelLine(from= "bus1", to = "bus5",y=1.025897454970189-1im*4.234983682334831, y_shunt_km=0.0492/2, y_shunt_mk=0.0492/2),
    "branch3"=> PiModelLine(from= "bus2", to = "bus3",y=1.1350191923073958-1im*4.781863151757718, y_shunt_km=0.0438/2, y_shunt_mk=0.0438/2),
    "branch4"=> PiModelLine(from= "bus2", to = "bus4",y=1.686033150614943-1im*5.115838325872083, y_shunt_km=0.0340/2, y_shunt_mk=0.0340/2),
    "branch5"=> PiModelLine(from= "bus2", to = "bus5",y=1.7011396670944048-1im*5.193927397969713, y_shunt_km=0.0346/2, y_shunt_mk=0.0346/2),
    "branch6"=> PiModelLine(from= "bus3", to = "bus4",y=1.9859757099255606-1im*5.0688169775939205, y_shunt_km=0.0128/2, y_shunt_mk=0.0128/2),
    "branch7"=> StaticLine(from= "bus4", to = "bus5",Y=6.840980661495672-1im*21.578553981691588),
    "branch8"=> Transformer(from= "bus4", to = "bus7", y=0.0-1im*4.781943381790359, t_ratio=0.978),
    "branch9"=> Transformer(from= "bus4", to = "bus9", y=0.0-1im*1.7979790715236075, t_ratio=0.969),
    "branch10"=> Transformer(from= "bus5", to = "bus6", y=0.0-1im*3.967939052456154, t_ratio=0.932),
    "branch11"=> StaticLine(from= "bus6", to = "bus11",Y=1.9550285631772604-1im*4.0940743442404415),
    "branch12"=> StaticLine(from= "bus6", to = "bus12",Y=1.525967440450974-1im*3.1759639650294003),
    "branch13"=> StaticLine(from= "bus6", to = "bus13",Y=3.0989274038379877-1im*6.102755448193116),
    "branch14"=> StaticLine(from= "bus7", to = "bus8",Y=0.0-1im*5.676979846721544),
    "branch15"=> StaticLine(from= "bus7", to = "bus9",Y=0.0-1im*9.09008271975275),
    "branch16"=> StaticLine(from= "bus9", to = "bus10",Y=3.902049552447428-1im*10.365394127060915),
    "branch17"=> StaticLine(from= "bus9", to = "bus14",Y=1.4240054870199312-1im*3.0290504569306034),
    "branch18"=> StaticLine(from= "bus10", to = "bus11",Y=1.8808847537003996-1im*4.402943749460521),
    "branch19"=> StaticLine(from= "bus12", to = "bus13",Y=2.4890245868219187-1im*2.251974626172212),
    "branch20"=> StaticLine(from= "bus13", to = "bus14",Y=1.1369941578063267-1im*2.314963475105352));


return PowerGrid(buses, branches)
end