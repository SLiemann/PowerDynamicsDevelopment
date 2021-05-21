using PowerDynamics
using OrderedCollections: OrderedDict

function OLTC_Hybrid_Sensi()
    Ubase = 380e3
    Sbase = 100e6
    Zbase = (Ubase^2) / (Sbase)
    buses = OrderedDict(
        "bus1" => SlackAlgebraic(U = 1.05),
        "bus2" => SimpleRecoveryLoad(
            P0 = -0.4,
            Q0 = 0.0,
            Pt = -0.4,
            Qt = 0.0,
            Tp = 5.0,
            Tq = 5.0,
        ),
    )

    branches = OrderedDict(
        #"branch1"=> DynamicPowerTransformer(from="bus1",to="bus2",Sbase=Sbase,Srated=100e6,uk=0.1581138,XR_ratio=3, i0=6.35,Pv0=300e3,
        #                                   tap_side = "LV",tap_pos = 0,tap_inc = 1.0,tap_delay = 5.0,tap_max = 10,tap_min = -10,
        #                                   deadband_low = 0.99,deadband_high = 1.01, timer_start = -1,Δtap = 0.0,low_high = 0.0))
        "branch1" => StaticPowerTransformerTapParam(
            from = "bus1",
            to = "bus2",
            Sbase = Sbase,
            Srated = 100e6,
            uk = 0.1581138,
            XR_ratio = 3,
            i0 = 6.35,
            Pv0 = 300e3,
            tap_side = "LV",
            tap_pos = 0,
            tap_inc = 1.0,
            tap_max = 5,
            tap_min = -10,
        ),
    )
    pg = PowerGrid(buses, branches)
end
