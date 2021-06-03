using PowerDynamics
using OrderedCollections: OrderedDict

function GFC_Test_Grid()
    Ubase = 320e3
    Sbase = 1000e6
    Zbase = (Ubase^2) / (Sbase)
    buses = OrderedDict(
        "bus1" => SlackAlgebraic(U = 1.00),
        "bus2" => VoltageDependentLoad(
            P = 0.0,
            Q = 0.0,
            U = 1.0,
            A = 0.0,
            B = 0.0,
            Y_n = 0.0,
        ),
        "bus3" => GridFormingConverter(
            Sbase = Sbase,
            Srated = Sbase,
            p0set = 0.5,
            u0set = 1.01,
            Kp_droop = 0.04,
            Kq_droop = 0.04,
            ωf = 10.0*2*pi,
            lf = 0.15,
            rf = 0.005,
            cf = 15.51,
            Kp_u = 1.0,
            Ki_u = 1.16,
            Kp_i = 0.73,
            Ki_i = 1.19,
        ),
    )

    branches = OrderedDict(
        "branch1" => PiModelLine(from="bus1",to="bus2",y=1.0/(1im*0.1),y_shunt_km=0.0,y_shunt_mk=0.0),
        "branch2" => PiModelLine(from="bus2",to="bus3",y=1.0/(1im*0.1),y_shunt_km=0.0,y_shunt_mk=0.0),
    )
    pg = PowerGrid(buses, branches)
end
