using PowerDynamics
using DifferentialEquations
using Plots
using OrderedCollections

begin
    include("C:/Users/liemann/github/PowerDynamicsDevelopment/src/include_custom_nodes_lines_utilities.jl")
    include("C:/Users/liemann/github/PowerDynamicsDevelopment/src/operationpoint/InitializeInternals.jl")
end
begin
    Ubase = 1e3;
    Sbase = 0.5e6;
    Zbase = Ubase^2/Sbase;

    Rdc = 1.2
    Rdc_pu = 1.2/Zbase
    Gdc = 1.0/Rdc_pu
    Xcdc = 1.0/(100*pi*0.008) / Zbase
    Cdc = 1.0/(100*pi*Xcdc)

    R_f = 0.001 /Zbase;
    L_f = 200*10^-6;
    Xlf = L_f * 100*pi /Zbase
    C_f = 300*10^-6;
    Xcf = 1.0/(100*pi*C_f) /Zbase
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
        "bus3" => MatchingControl(
            Sbase = 100e6,
            Srated = 100e6,
            p0set = 0.8, # based on Sbase!
            u0set = 1.00,
            Kp_uset = 0.001,
            Ki_uset = 0.5,
            Kdc = 1600.0,
            gdc = Gdc,
            cdc = Cdc,
            xlf = Xlf,    #0.15
            rf =  R_f, #0.0005
            xcf = Xcf, #15.1515151515
            Kp_u = 0.52, #1.0
            Ki_u = 232.2, #1.161022
            Kp_i = 0.73, # 0.738891
            Ki_i = 0.0059,
            imax_csa = 10.0,
            p_ind = [1:16]
        ),
    )

    branches = OrderedDict(
        "branch1" => PiModelLine(
            from = "bus1",
            to = "bus2",
            y = 1.0/(0.005+1im*0.05),#1.97712-1im*19.78677,
            y_shunt_km = 0.0,
            y_shunt_mk = 0.0,
        ),
        "branch2"=> PiModelLine(
            from = "bus2",
            to = "bus3",
            y = 1.0/(0.005+1im*0.05),#1.97712-1im*19.78677,
            y_shunt_km = 0.0,
            y_shunt_mk = 0.0,
        ))
    pg = PowerGrid(buses, branches)
end
begin
    U,δ,ic0 = PowerFlowClassic(pg, iwamoto = true, max_tol = 1e-7)
    pg1 ,ic = InitializeInternalDynamics(pg,ic0)
    params = GFC_params()
    prob = ODEProblem(rhs(pg1),ic,(0.0,2.0),params, initializealg = BrownFullBasicInit())
end
