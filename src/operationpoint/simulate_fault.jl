function my_simulate(np::AbstractPerturbation, powergrid::PowerGrid, x1::Array, timespan; solve_kwargs...)
    @assert first(timespan) <= np.tspan_fault[1] "fault cannot begin in the past"
    @assert np.tspan_fault[2] <= last(timespan) "fault cannot end in the future"

    np_powergrid = np(powergrid)

    problem = ODEProblem{true}(rhs(powergrid), x1, timespan)
    #include("operationpoint/PowerFlow.jl")

    function errorState(integrator)
        sol1 = integrator.sol
        #x2 = find_valid_initial_condition(np_powergrid, sol1[end]) # Jump the state to be valid for the new system.
        ode = rhs(np_powergrid)
        op_prob = ODEProblem(ode, sol1[end], (0.0, 1e-6), nothing, initializealg = BrownFullBasicInit())
        x2 = solve(op_prob,Rodas5())
        x2 = x2.u[end]

        integrator.f = rhs(np_powergrid)
        integrator.u = x2#sol1[end]
    end

    function regularState(integrator)
        sol2 = integrator.sol
        #x3 = find_valid_initial_condition(powergrid, sol2[end]) # Jump the state to be valid for the new system.
        ode = rhs(powergrid)
        ic_tmp = deepcopy(integrator.sol.u[indexin(np.tspan_fault[1],integrator.sol.t)[1]]) #get ic from pre-fault
        ic_tmp = getPreFaultVoltages(pg,ic_tmp,deepcopy(sol2[end])) #change only voltages
        op_prob = ODEProblem(ode, ic_tmp, (0.0, 1e-6), nothing, initializealg = BrownFullBasicInit())
        x3 = solve(op_prob,Rodas5())
        x3 = x3.u[end]

        integrator.f = rhs(powergrid)
        integrator.u = x3#sol2[end]
    end

    t1 = np.tspan_fault[1]
    t2 = np.tspan_fault[2]

    cb1 = DiscreteCallback(((u,t,integrator) -> t in np.tspan_fault[1]), errorState)
    cb2 = DiscreteCallback(((u,t,integrator) -> t in np.tspan_fault[2]), regularState)

    sol = solve(problem, Rodas4(), callback = CallbackSet(cb1, cb2), tstops=[t1, t2], solve_kwargs...)

    return PowerGridSolution(sol, powergrid)
end
