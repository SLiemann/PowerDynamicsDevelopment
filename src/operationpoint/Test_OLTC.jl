using PowerDynamics: SlackAlgebraic, FourthOrderEq, VoltageDependentLoad, PiModelLine, StaticLine, Transformer, PowerGrid#, write_powergrid, Json, Inc, find_operationpoint, ChangeInitialConditions, LineFailure, PowerPerturbation, simulate
using PowerDynamics: symbolsof, initial_guess, guess, RLLine
using PowerDynamics: StaticPowerTransformer, DynamicPowerTransformer, SixOrderMarconatoMachine,SixOrderMarconatoMachineSin
using PowerDynamics
using PowerDynamics: rhs, State
using OrderedCollections: OrderedDict
using Plots
import PowerDynamics: PiModel
using DifferentialEquations
using CSV #read PF DataFrames
using DataFrames #for CSV
using ModelingToolkit

begin
    include("operationpoint/PowerFlow.jl")
    Ubase = 380e3
    Sbase = 100e6
    buses=OrderedDict(
        "bus1"=> SlackAlgebraic(U=0.98),
        "bus2"=> VoltageDependentLoad(P=-0.3, Q=0.3, U=1.0, A=0., B=0.,Y_n = complex(0.0)),
        "bus3"=> SixOrderMarconatoMachineSin(H = 5, P=0.5, D=0., Ω=50, E_f=1.4, R_a = 0.1,T_ds=1.136,T_qs=0.8571,T_dss=0.04,T_qss=0.06666,X_d=1.1,X_q=0.7,X_ds=0.25,X_qs=0.25,X_dss=0.2,X_qss=0.2,T_AA=0.),
        #"bus3"=> VoltageDependentLoad(P=+0.5, Q=0.0, U=1.0, A=0., B=0.,Y_n = complex(0.0)),
        "bus4"=> VoltageDependentLoad(P=-0.5, Q=-0.5, U=1.0, A=0.5, B=0.3,Y_n = complex(0.0)),
        "bus5"=> VoltageDependentLoad(P= 0.0, Q=0.1, U=1.0, A=0.5, B=0.3,Y_n = complex(0.0)))

    branches=OrderedDict(
        "branch1"=> PiModelLine(from= "bus1", to = "bus2",y=1.0/(0.05+1im*0.15), y_shunt_km=0., y_shunt_mk=0.),
        "branch2"=> PiModelLine(from= "bus2", to = "bus3",y=1.0/(0.05+1im*0.15), y_shunt_km=0., y_shunt_mk=0.),
        "branch3"=> StaticPowerTransformer(from="bus3",to="bus4",S_r=100e6,U_r=380e3,uk=0.1581138,XR_ratio=3,i0=6.35,Pv0=100e3,Sbase=Sbase,Ubase=Ubase,tap_side = "LV",tap_pos = 3,tap_inc = 1.0),
        "branch4"=> PiModelLine(from= "bus4", to = "bus5",y=1.0/((0.05+1im*0.15)*(Ubase/110e3)^2), y_shunt_km=0., y_shunt_mk=0.))
    pg = PowerGrid(buses, branches)
    Qmax   = [Inf, Inf, 0.5,Inf, Inf]
    Qmin   = -Qmax
    U,δ1,ic = PowerFlowClassic(pg,iwamoto = true, Qmax = Qmax, Qmin = Qmin, Qlimit_iter_check = 2)
end
#ic = initial_guess(pg)
#ic0 = find_valid_initial_condition(pg,ic)
begin
    include("operationpoint/InitializeInternals.jl")
    Uc = U.*exp.(1im*δ1/180*pi)
    Ykk = NodalAdmittanceMatrice(pg)
    I_c = Ykk*Uc
    PG, ic0 = InitializeInternalDynamics(pg,I_c,ic)
    ODEProb = ODEProblem(rhs(PG),ic0,(0.0,10.0))
    new_f = ODEFunction(ODEProb.f.f, syms = ODEProb.f.syms, mass_matrix = Int.(ODEProb.f.mass_matrix))
    ODEProb = ODEProblem(new_f,ic,(0,1.0), 1) # 1 is a dummy parameter
    mtsys = modelingtoolkitize(ODEProb)
    Jakob = ModelingToolkit.calculate_jacobian(mtsys)
    Jakob2 = calc_my_jac(mtsys)
    Jakob .== Jakob2
    #sol = solve(ODEProb,Rodas4(),dt=1e-4)
    #Zbase = (380e3^2)/(100e6)
    #SS = NodeShortCircuit(node="bus2",Y = 1/(1im*250/Zbase),tspan_fault=(1.0, 1.15))
    #PT = PowerPerturbation(node="bus2", fault_power = -0.4 ,tspan_fault=(1.0, 5.0))
    #PGsol = my_simulate(SS,PG,ic0,(0.,5.0))
    #PG_state = State(PG,ic0)
    #PGsol = solve(PG,PG_state,(0.,1.))
    #plot(PGsol,collect(keys(PG.nodes)), :v,size = (1000, 500),legend = (0.5, 0.5))
end

begin
plot(sol)
plot(PGsol,collect(keys(PG.nodes)), :v,size = (1000, 500),legend = (0.5, 0.5))
begin
    plot(PGsol,collect(keys(PG.nodes)), :v,size = (1000, 500),legend = (0.5, 0.5))
    test = DataFrame(CSV.File("C:\\Users\\liemann\\Desktop\\PF_test.csv"; header=false, delim=';', type=Float64))
    plot!(test.Column1,test.Column2,label = "PF-bus1")
    plot!(test.Column1,test.Column3,label = "PF-bus2")
    plot!(test.Column1,test.Column4,label = "PF-bus3")
    plot!(test.Column1,test.Column5,label = "PF-bus4")
    plot!(test.Column1,test.Column6,label = "PF-bus5")
end
begin
    plot(PGsol,["bus3"], :θ)
    test2 = DataFrame(CSV.File("C:\\Users\\liemann\\Desktop\\PF_pol.csv"; header=false, delim=';', type=Float64))
    plot!(test2.Column1,(test2.Column2))
    #xlims!((0.9,10.))
end
begin
    plot(PGsol,["bus3"], :ω)
    test3 = DataFrame(CSV.File("C:\\Users\\liemann\\Desktop\\PF_dreh.csv"; header=false, delim=';', type=Float64))
    plot!(test3.Column1,(test3.Column2.-1))
    #ylims!((-0.05,0.05))
end
#ylims!((1.01,1.015))
#ylims!((0.75,1.015))
xlims!((0.99,1.05))

plot(PGsol,["bus1"], :p)
plot!(PGsol,["bus1"], :q)


function OLTC_test(powergrid::PowerGrid, ic0::Array,branch_oltc)

    problem = ODEProblem{true}(rhs(powergrid), ic0, (0.0,50.0))
    np_pg = deepcopy(powergrid)
    timer_start = -1
    timer_now   = 0.0
    tap = PG.lines[branch_oltc].tap_pos
    OLTC = deepcopy(PG.lines[branch_oltc])
    function TapState(integrator)
        timer_start = integrator.t
        sol1 = integrator.sol
        tap += 1
        node = StaticPowerTransformer(from=OLTC.from,to=OLTC.to,S_r=OLTC.S_r,U_r=OLTC.U_r,uk=OLTC.uk,XR_ratio=OLTC.XR_ratio,i0=OLTC.i0,Pv0=OLTC.Pv0,Sbase=OLTC.Sbase,Ubase=OLTC.Ubase,tap_side = OLTC.tap_side,tap_pos = tap,tap_inc = OLTC.tap_inc)
        np_pg.lines[branch_oltc] = node

        ode = rhs(np_pg)
        op_prob = ODEProblem(ode, sol1[end], (0.0, 1e-6), nothing, initializealg = BrownFullBasicInit())
        x2 = solve(op_prob,Rodas5())
        x2 = x2.u[end]

        integrator.f = rhs(np_pg)
        integrator.u = x2#sol1[end]
    end

    function voltage_deadband(u,t,integrator)
         0.98 <= sqrt(u[13]*u[13] + u[14]*u[14]) <= 1.02
    end

    function timer_off(integrator)
        timer_start = -1
    end

    function voltage_outside(u,t,integrator)
         sqrt(u[13]*u[13] + u[14]*u[14]) < 0.98
    end

    function timer_on(integrator)
        if timer_start == -1
            timer_start = integrator.t
        end
    end

    function timer_hit(u,t,integrator)
        if timer_start == -1
            return false
        else
            return t-timer_start > 5
        end
    end

    cb1 = DiscreteCallback(voltage_deadband, timer_off)
    cb2 = DiscreteCallback(voltage_outside, timer_on)
    cb3 = DiscreteCallback(timer_hit, TapState)

    sol = solve(problem, Rodas4(), callback = CallbackSet(cb1,cb2,cb3), dt = 1e-3,adaptive=false)

    return PowerGridSolution(sol, powergrid)
end
rhs()

begin
    PGsol2 = OLTC_test(PG,ic0,"branch3")
    plot(PGsol2,collect(keys(PG.nodes)), :v,size = (1000, 500),legend = (0.6, 0.75))
    test = DataFrame(CSV.File("C:\\Users\\liemann\\Desktop\\PF_test.csv"; header=false, delim=';', type=Float64))
    plot!(test.Column1,test.Column2,label = "PF-bus1")
    plot!(test.Column1,test.Column3,label = "PF-bus2")
    plot!(test.Column1,test.Column4,label = "PF-bus3")
    plot!(test.Column1,test.Column5,label = "PF-bus4")
    plot!(test.Column1,test.Column6,label = "PF-bus5")
end
begin
    plot(PGsol2,["bus3"], :θ)
    test2 = DataFrame(CSV.File("C:\\Users\\liemann\\Desktop\\PF_pol.csv"; header=false, delim=';', type=Float64))
    plot!(test2.Column1,(test2.Column2))
end

function calc_my_jac(sys)
    rhs1 = [eq.rhs for eq ∈ equations(sys)]
    dvs = states(sys)
    Jakob = [Num(expand_derivatives(Differential(v)(lhs))) for lhs in rhs1, v in dvs]
end

function calc_my_hesse(sys)
    rhs1 = [eq.rhs for eq ∈ equations(sys)]
    dvs = [Differential(s)^2 for s in states(sys)]
    hesse = [Num(expand_derivatives(v(lhs))) for lhs in rhs1, v in dvs]
end

test = calc_my_jac(mtsys)
test2 = calc_my_hesse(mtsys)
