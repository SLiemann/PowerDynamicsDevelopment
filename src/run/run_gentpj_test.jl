using PowerDynamics
using DifferentialEquations
using OrderedCollections: OrderedDict
using Plots

include("C:/Users/liemann/github/PowerDynamicsDevelopment/src/include_custom_nodes_lines_utilities.jl")
function Testgrid()
    buses=OrderedDict(
            "bus1" => SlackAlgebraic(U=1.0),
            "bus2" => VoltageDependentLoad(P=-0.5,Q=0.0,U=1.0,A=0,B=0),
           "bus3" => gentpj(Sbase=100e6, Srated=100e6, H=6.0, P=1.0, D=0.1, Ω=50, E_fd=1.0, R_a=0, T_d0s=7.0, T_q0s=1.5, T_d0ss=0.05, T_q0ss=0.05, X_d=2.2, X_q=2.0, X_ds=0.3, X_qs=0.4, X_dss=0.2, X_qss=0.2, X_l=0.15, S_10=0.1, S_12=0.3,K_is=0.0))


    branches=OrderedDict(
        "Line_1-2"=> PiModelLine(from= "bus1", to = "bus2",y=1.0/(0.005+1im*0.05), y_shunt_km=0.0, y_shunt_mk=0.0),
        "Line_2-3"=> PiModelLine(from= "bus2", to = "bus3",y=1.0/(0.005+1im*0.05), y_shunt_km=0.0, y_shunt_mk=0.0),)
    return PowerGrid(buses, branches)
end


U,δ1,ic0,uc = PowerFlowClassic(pg,iwamoto = true,max_tol=1e-8)
pg_tmp,ic = InitializeInternalDynamics(pg,ic0)


prob = ODEProblem{true}(rhs(pg_tmp),ic,(0.0,100.0))
sol = solve(prob,Rosenbrock23());
plot(sol,vars=(12))
plot(sol,vars=(13))
for i=1:13
    display(plot(sol,vars=(i)))
end


function simTestGrid(;p_new = -0.5)
    #Regular and fault powergrid
    pg_regular = Testgrid()

    #parameters for simulation
    #tstep = [5,10]
    tstep = [1.00,5.0]
    #tstep = [0.1,0.2] #time points of fault
    timespan = (0.0,10.0) #duration of simulation #1.5

    #Find initial condictions
    U,δ,ic0 = PowerFlowClassic(pg_regular, iwamoto = true, max_tol = 1e-7)
    pg_regular,ic = InitializeInternalDynamics(pg_regular,ic0)
    pg_fault = deepcopy(pg_regular)
    pg_fault.nodes["bus2"] = VoltageDependentLoad(P=p_new,Q=0.0,U=1.0,A=0,B=0)

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
    #sol = solve(problem, Rodas4(),callback = CallbackSet(cb1,cb2), dtmax = 1e-3,progress=true, initializealg = BrownFullBasicInit())
    sol = solve(problem, Rodas4(),callback=CallbackSet(cb1,cb2), dtmax = 1e-3,progress=true, initializealg = BrownFullBasicInit())

    #returning the solution
    return PowerGridSolution(sol,pg_regular)
end

pgsol = simTestGrid(p_new = -0.6);

plot!(pgsol,["bus1","bus3","bus3"],:v, legend = (0.9,0.5), label = ["U1" "U2" "U3"])

plot(pgsol,"bus1",:v, legend = (0.9,0.5), label = "U1")
plot(pgsol,"bus2",:v, legend = (0.9,0.5), label = "U2")
plot(pgsol,"bus3",:v, legend = (0.9,0.5), label = "U3")
plot(pgsol,"bus3",:θ, legend = (0.9,0.5), label = "θ")
plot(pgsol,["bus3"],:p, legend = (0.9,0.5), label = "active power")
plot(pgsol,["bus3"],:q, legend = (0.9,0.5), label = "reactive power")
plot(pgsol,"bus3",:ω, legend = (0.9,0.5), label = "omega")

