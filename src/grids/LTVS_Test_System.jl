using PowerDynamics
using OrderedCollections: OrderedDict


function LTVS_Test_System()
    Ubase = 380e3
    Sbase = 100e6
    Zbase = Ubase^2/Sbase
    buses=OrderedDict(
        "bus1" => SlackAlgebraic(U=1.0),
        "busv"=> VoltageDependentLoad(P=0.0 Q=0.0, U=1.0, A=0., B=0.,Y_n = complex(0.0)),
        "bus2" => VoltageDependentLoad(P=0.0 Q=0.0, U=1.0, A=0., B=0.,Y_n = complex(0.0)),
        "bus3" => VoltageDependentLoad(P=-11.0, Q = -3.0, U=1.0, A=0.0, B=0.0,Y_n = complex(0.0)),
        "bus4" => SixOrderMarconatoMachineAVROEL(H = 3, P=6.0, D=0., Ω=50, R_a = 0.0,
                                             T_ds=0.9545455,T_qs=0.3T_dss=0.0333333,T_qss=0.025,
                                             X_d=2.2,X_q=2.0,X_ds=0.3,X_qs=0.4,X_dss=0.2,
                                             X_qss=0.2,T_AA=0.,V0 = 1.0, Ifdlim = 3.0618,
                                             L1 = -18.0, G1 = 120.0, Ta = 5.0, Tb = 12.5,
                                             G2 = 10.0, L2 = 5.0))
    Z_EHV_Line = (9.6 + 1im*64)/Zbase
    B_half     = 1im*1498.54*10e-6 / 2.0 *Zbase #already an admittance
    branches=OrderedDict(
        "Line_1-2"=> PiModelLine(from= "bus1", to = "bus2",y=1.0/Z_EHV_Line, y_shunt_km=B_half, y_shunt_mk=B_half),
        "Line_1-v"=> PiModelLine(from= "bus1", to = "busv",y=2.0/Z_EHV_Line, y_shunt_km=B_half/2.0, y_shunt_mk=B_half/2.0,),
        "Line_v-2"=> PiModelLine(from= "busv", to = "bus2",y=2.0/Z_EHV_Line, y_shunt_km=B_half/2.0,, y_shunt_mk=B_half/2.0,),
        "branch3"=> StaticPowerTransformer(from="bus3",to="bus4",S_r=100e6,U_r=380e3,uk=0.1581138,XR_ratio=3,
                                           i0=6.35,Pv0=100e3,Sbase=Sbase,Ubase=Ubase,tap_side = "LV",tap_pos = 3,tap_inc = 1.0),
        "branch4"=> PiModelLine(from= "bus4", to = "bus5",y=1.0/((0.05+1im*0.15)*(Ubase/110e3)^2), y_shunt_km=0., y_shunt_mk=0.))
        return PowerGrid(buses, branches)
end
