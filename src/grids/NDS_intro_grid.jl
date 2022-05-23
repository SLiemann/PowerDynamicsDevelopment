using PowerDynamics
using OrderedCollections: OrderedDict
using Plots

include("C:/Users/liemann/github/PowerDynamicsDevelopment/src/include_custom_nodes_lines_utilities.jl")

function GetTestGrid(;y_new = 0.0)
    buses = OrderedDict(
        "bus1" => SlackAlgebraic(U = 1.0),
        "bus2" => VoltageDependentLoad(P=-0.50,Q=0.1,U=1.0,A=1.0,B=0.0,Y_n = 0.0),
        "bus3" => SixOrderMarconatoMachine(Sbase=600e6,Srated=600e6,H = 0.50, P=0.5, D=0.0, Ω=50,E_f = 1.0, R_a = 0.0,
                                             T_ds=0.9545455,T_qs=0.3,T_dss=0.0333333,T_qss=0.025,
                                             X_d=2.2,X_q=2.0, X_ds=0.3,X_qs=0.4, X_dss=0.2,
                                             X_qss=0.2,T_AA=0.0),
        "bus4" => VoltageDependentLoad(P=-0.25,Q=0.0,U=1.0,A=1.0,B=0.0,Y_n = y_new),
    )
    branches = OrderedDict(
        "branch12" => PiModelLine(from="bus1",to="bus2",y=1.0/(1im*0.1),y_shunt_km=0.0,y_shunt_mk=0.0),
        "branch13" => PiModelLine(from="bus1",to="bus3",y=1.0/(1im*0.1),y_shunt_km=0.0,y_shunt_mk=0.0),
        "branch23" => PiModelLine(from="bus2",to="bus3",y=1.0/(1im*0.1),y_shunt_km=0.0,y_shunt_mk=0.0),
        "branch24" => PiModelLine(from="bus2",to="bus4",y=1.0/(1im*0.1),y_shunt_km=0.0,y_shunt_mk=0.0),
    )
    return pg = PowerGrid(buses, branches)
end

function simTestGrid()
    #Regular and fault powergrid
    pg_regular = GetTestGrid()
    pg_fault = GetTestGrid(y_new = 8.0)

    #parameters for simulation
    tstep = [0.1,0.20] #time points of ault
    timespan = (0.0,1.5) #duration of simulation

    #Find initial condictions
    U,δ,ic0 = PowerFlowClassic(pg_regular, iwamoto = true, max_tol = 1e-7)
    pg_regular,ic = InitializeInternalDynamics(pg_regular,ic0)

    #get equations
    regular = rhs(pg_regular)
    error = rhs(pg_fault)

    #Function-wrapper to distinguish between regular and fault powergrid
    _f = (dx, x, p, t) -> p ? regular(dx,x,nothing,t) : error(dx,x,nothing,t)
    #Create ODEFunction for ODEProblem
    f = ODEFunction(_f, mass_matrix = regular.mass_matrix, syms = regular.syms)
    #Build solvable ODEProblem with initial conditions, timespan and parameter (regular grid first)
    problem = ODEProblem{true}(f, ic, timespan, true)

    function errorState(integrator) #here, mainly the parameter is set to take the equations of the fault powergrid in the function wrapper
        integrator.p = false
        if integrator.opts.adaptive
            auto_dt_reset!(integrator)
            set_proposed_dt!(integrator, integrator.dt)
        end
        #Calculation of new initial conditions
        sol = integrator.sol
        op_prob = ODEProblem(error, sol.u[end], (0.0, 1e-6),[], initializealg = BrownFullBasicInit())
        ic_new = solve(op_prob,Rodas4())
        integrator.u = ic_new.u[end]
    end

    function regularState(integrator) #similar to errorstate but now the regular state again
        integrator.p = true
        if integrator.opts.adaptive
            auto_dt_reset!(integrator)
            set_proposed_dt!(integrator, integrator.dt)
        end
        #Calculation of new initial conditions with information about pre-fault values
        sol = integrator.sol
        ic_tmp = deepcopy(integrator.sol.u[indexin(tstep[1],integrator.sol.t)[1]])
        ic_tmp = getPreFaultAlgebraicStates(pg_regular,ic_tmp,deepcopy(sol[end]))
        op_prob = ODEProblem(regular, ic_tmp, (0.0, 1e-6),[], initializealg = BrownFullBasicInit())
        ic_new = solve(op_prob,Rodas4())
        integrator.u = ic_new.u[end]
    end

    #Callbacks
    #PresetTimeCallback creates callbacks at the specified points in time
    #If the point in time is hit it evaluates the given function
    cb1 = PresetTimeCallback([tstep[1]], errorState)
    cb2 = PresetTimeCallback([tstep[2]], regularState)

    #Solving the problem with given solver, callbacks and other parameter to control the calculation
    sol = solve(problem, Rodas4(),callback = CallbackSet(cb1,cb2), dtmax = 1e-3,progress=true, initializealg = BrownFullBasicInit())

    #returning the solution
    return PowerGridSolution(sol,pg_regular)
end

pgsol = simTestGrid()

plot(pgsol,["bus1","bus3","bus3","bus4"],:v, legend = (0.9,0.5), label = ["U1" "U2" "U3" "U4"])

plot(pgsol,"bus1",:v, legend = (0.9,0.5), label = "U1")
plot(pgsol,"bus2",:v, legend = (0.9,0.5), label = "U2")
plot(pgsol,"bus3",:v, legend = (0.9,0.5), label = "U3")
plot(pgsol,"bus4",:v, legend = (0.9,0.5), label = "U4")
plot(pgsol,"bus3",:θ, legend = (0.9,0.5), label = "θ")
plot(pgsol,"bus3",:ω, legend = (0.9,0.5), label = "omega")
plot(pgsol,["bus3"],:p, legend = (0.9,0.5), label = "active power")
plot!(pgsol,["bus3"],:q, legend = (0.9,0.5), label = "reactive power")








Uc1 = U.*exp.(1im*δ/180*pi)
