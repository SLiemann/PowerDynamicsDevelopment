using PowerDynamics
using DifferentialEquations
using Plots
using OrderedCollections

begin
    include("C:/Users/liemann/github/PowerDynamicsDevelopment/src/include_custom_nodes_lines_utilities.jl")
    include("C:/Users/liemann/github/PowerDynamicsDevelopment/src/operationpoint/InitializeInternals.jl")
end
function getMachtingGrid(;y_new = 0.0)
    Ubase = 1e3;
    Sbase = 100e6;
    Zbase = Ubase^2/Sbase;

    Rdc = 1.2
    Rdc_pu = 1.2/Zbase
    Gdc = 1.0/Rdc_pu
    Xcdc = 1.0/(100*pi*0.008*200) / Zbase
    Cdc = 1.0/Xcdc#(100*pi*Xcdc)

    R_f = 0.001/200/Zbase;
    L_f = 1*10^-6;
    Xlf = L_f * 100*pi /Zbase
    C_f = 200*300*10^-6;
    Xcf = 1.0/(100*pi*C_f) /Zbase
    buses = OrderedDict(
        "bus1" => SlackAlgebraic(U = 1.00),
        "bus2" => VoltageDependentLoad(
            P = 0.0,
            Q = 0.0,
            U = 1.0,
            A = 0.0,
            B = 0.0,
            Y_n = y_new,
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
            xlf = Xlf,    #  0.15
            rf = R_f, # 0.0005
            xcf =  Xcf ,# 15.1515151515
            Tdc = 0.05,
            Kp_u = 0.52*200, #1.0
            Ki_u = 232.2, #1.161022#
            Kp_i = 0.73/200, # 0.738891
            Ki_i = 0.0059,
            imax_csa = 1.90,
            p_ind = collect(1:14)
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
    pg = getMachtingGrid()
    U,δ,ic0 = PowerFlowClassic(pg, iwamoto = true, max_tol = 1e-7)
    pg1 ,ic = InitializeInternalDynamics(pg,ic0)
    params = getallParameters(pg.nodes["bus3"])[5:18]
    prob = ODEProblem(rhs(pg1),ic,(0.0,5.0),params)#initializealg = BrownFullBasicInit()
    #sol = solve(prob, Rodas4(),dtmax = 1e-3,progress=true)
    #pgsol = PowerGridSolution(sol,pg1)
    pgsol = simMatching(prob)
end

plot(pgsol.dqsol,vars=(7:13))
plot(pgsol.dqsol,vars=(9))
plot(pgsol,"bus1",:iabs)
plot(pgsol,"bus3",:udc)
plot(pgsol,"bus3",:idc0)
plot(pgsol,"bus3",:θ)
plot(pgsol,"bus3",:x_uabs)

rhs(pg).syms

θ = ExtractResult(pgsol,:θ_3).*180/pi
plot(t,θ)

t = pgsol.dqsol.t
uc = pgsol(pgsol.dqsol.t,"bus3",:u)
i_c = pgsol(pgsol.dqsol.t,"bus3",:i)'
s = uc.*conj(i_c)
p = real(s)
q = imag(s)
plot(t,p)
plot(t,q)
plot(t,abs.(i_c))

function simMatching(prob)
    pg_new = getMachtingGrid(y_new = 1/0.5)
    tstep = [1.0]
    function fault_state(integrator)
        new_f = rhs(pg_new)
        op_prob = ODEProblem(new_f, integrator.sol[end], (0.0, 1e-6),integrator.p, initializealg = BrownFullBasicInit())
        ic_new = solve(op_prob,Rodas5())
        integrator.f = new_f
        integrator.cache.tf.f = integrator.f
        integrator.u = ic_new.u[end]
    end

    #=function postfault_state(integrator)
        sol = integrator.sol
        ic_tmp = deepcopy(integrator.sol.u[indexin(tstep[1],integrator.sol.t)[1]])
        ic_tmp = getPreFaultVoltages(pg_new,ic_tmp,deepcopy(sol[end]))
        #ic_tmp = getPreFaultAlgebraicStates(pg_new,ic_tmp,deepcopy(sol[end]))
        op_prob = ODEProblem(prob.f, ic_tmp, (0.0, 1e-6),integrator.p, initializealg = BrownFullBasicInit())
        ic_new = solve(op_prob,Rodas5())
        integrator.f = prob.f
        integrator.cache.tf.f = integrator.f
        integrator.u = ic_new.u[end]
        event_recorder = vcat(event_recorder,[integrator.t 1 integrator.p' 2 1])
    end=#

    cb = DiscreteCallback(((u,t,integrator) -> t in tstep[1]), fault_state)
    #cb1 = DiscreteCallback(((u,t,integrator) -> t in tstep[2]), postfault_state)

    sol = solve(prob, Rodas4(), tstops= tstep,callback = cb, dtmax = 1e-3,progress=true)
    #return sol, event_recorder
    return PowerGridSolution(sol,pg_new)
end
