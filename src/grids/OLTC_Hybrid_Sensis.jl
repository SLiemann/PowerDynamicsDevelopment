using PowerDynamics
using OrderedCollections: OrderedDict
using ModelingToolkit

include("C:/Users/liemann/github/PowerDynamicsDevelopment/src/operationpoint/PowerFlow.jl")
include("C:/Users/liemann/github/PowerDynamicsDevelopment/src/utility/utility_functions.jl")
include("C:/Users/liemann/github/PowerDynamicsDevelopment/src/operationpoint/InitializeInternals.jl")

function OLTC_Hybrid_Sensi(;x_grid = 0.25)
    Ubase = 380e3
    Sbase = 100e6
    Zbase = (Ubase^2) / (Sbase)
    buses = OrderedDict(
        "bus1" => SlackAlgebraic(U = 1.05),
        "bus2" => VoltageDependentLoad(P=0.0,Q=0.0,U=1.0,A=0.0,B=0.0,Y_n = 0.0),
        "bus3" => VoltageDependentLoad(P=0.0,Q=0.0,U=1.0,A=0.0,B=0.0,Y_n = complex(0.0)),
        "bus4" => SimpleRecoveryLoad(
            P0 = -0.4,
            Q0 = 0.0,
            Pt = -0.4,
            Qt = 0.0,
            Tp = 5.0,
            Tq = 5.0,
        ),
    )

    branches = OrderedDict(
        "branch1" => PiModelLineParam(from="bus1",to="bus2",y=1.0/(1im*x_grid),y_shunt_km=0.0,y_shunt_mk=0.0,p_ind=1),
        "branch2" => StaticPowerTransformerTapParam(
            from = "bus2",
            to = "bus3",
            Sbase = Sbase,
            Srated = 100e6,
            uk = 0.01,
            XR_ratio = 999,
            i0 = 0.0,
            Pv0 = 0.0,
            tap_side = "LV",
            tap_pos = 3,
            tap_inc = 1.25,
            tap_max = 8,
            tap_min = -8,
            p_ind=2
        ),
        "branch3" => PiModelLine(from="bus3",to="bus4",y=1.0/(1im*0.8),y_shunt_km=0.0,y_shunt_mk=0.0),
    )
    pg = PowerGrid(buses, branches)
end

function GetInitializedOLTCHisken()
    pg = OLTC_Hybrid_Sensi()
    U,δ,ic0 = PowerFlowClassic(pg,iwamoto = true)
    Ykk = NodalAdmittanceMatrice(pg)
    Uc = U.*exp.(1im*δ/180*pi)
    I_c = Ykk*Uc
    S = conj(Ykk*Uc).*Uc
    return InitializeInternalDynamics(pg,I_c,ic0)
end

function GetMTKOLTCSystem(tspan,p)
    pg, ic = GetInitializedOLTCHisken()
    prob   = ODEProblem(rhs(pg),ic,tspan,p)
    new_f = ODEFunction(prob.f.f, syms = prob.f.syms, mass_matrix = Int.(prob.f.mass_matrix))
    ODEProb = ODEProblem(new_f,ic,tspan,p)
    return modelingtoolkitize(ODEProb)
end
