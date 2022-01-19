using PowerDynamics
using OrderedCollections: OrderedDict

function Testgrid_oPFC(;u_new = 1.0)
    buses = OrderedDict(
        "bus1" => SlackAlgebraic(U = u_new),
        "bus2" => oPFC(Cd = 1.0 /(100*pi*0.036)*2, Pdc = 1.0, Ulow = 0.3, Qn= -0.21, t0 = 0.0),
    )

    branches = OrderedDict(
        "branch1" => PiModelLine(
            from = "bus1",
            to = "bus2",
            y = 1.0/(0.00+1im*0.20),#1.97712-1im*19.78677,
            y_shunt_km = 0.0,
            y_shunt_mk = 0.0,),
        )
    pg = PowerGrid(buses, branches)
end
