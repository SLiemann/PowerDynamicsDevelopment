using PowerDynamics
using DifferentialEquations
using Plots
using ModelingToolkit
using ForwardDiff
using DiffResults


begin
    include("C:/Users/liemann/github/PowerDynamicsDevelopment/src/include_costum_nodes_lines_utilities.jl")
    include("C:/Users/liemann/github/PowerDynamicsDevelopment/src/grids/OLTC_Hybrid_Sensis.jl")
    #include("C:/Users/liemann/github/PowerDynamicsDevelopment/src/sensitivity_analyses/Local_Sensitivity.jl")
end

pg = OLTC_Hybrid_Sensi()
U,δ,ic0 = PowerFlowClassic(pg,iwamoto = false)
pg,ic = InitializeInternalDynamics(pg,ic0)
begin
    Tp = 5.0
    timer_start = -1.0
    event_recorder = Array{Float64,2}(undef,0,3+length(GetParametersOLTCHisken(5)))
    function TapState(integrator)
        timer_start = integrator.t
        tap_max = pg.lines["branch2"].tap_max
        tap_min = pg.lines["branch2"].tap_min
        tap_pos = pg.lines["branch2"].tap_pos
        if integrator.p[2] + tap_pos < tap_max && integrator.p[2] + tap_pos > tap_min
            integrator.p[2] += 1.0
            p_tmp = deepcopy(integrator.p)
            ic_tmp = deepcopy(integrator.sol[end])
            op_prob = ODEProblem(integrator.f, ic_tmp, (0.0, 1e-6), p_tmp, initializealg = BrownFullBasicInit())
            ic_new = solve(op_prob,Rodas5())
            integrator.u = ic_new.u[end]
        end
    end

    function voltage_deadband(u, t, integrator)
        1.04 <= sqrt(u[7] * u[7] + u[8] * u[8])
    end

    function timer_off(integrator)
        if integrator.p[5] != -1
            integrator.p[5] = -1
        end
    end

    function voltage_outside(u, t, integrator)
        sqrt(u[7] * u[7] + u[8] * u[8]) < 1.04 && t>10.0
    end

    function timer_on(integrator)
        if integrator.p[5] == -1
            integrator.p[5] = integrator.t
        end
    end

    function timer_hit(u, t, integrator)
        if integrator.p[5] == -1
            return false
        else
            return t - integrator.p[5] > 20.0
        end
    end

    function newGridImpedanze(integrator)
        integrator.p[1] = 0.65
        p_tmp = deepcopy(integrator.p)
        ic_tmp = deepcopy(integrator.sol[end])
        op_prob = ODEProblem(integrator.f, ic_tmp, (0.0, 1e-6), p_tmp, initializealg = BrownFullBasicInit())
        ic_new = solve(op_prob,Rodas5())
        integrator.u = ic_new.u[end]
    end

    function check_voltage_low(u,t,integrator)
            sqrt(u[7]*u[7] + u[8]*u[8]) < 0.45
    end
    function check_voltage_high(u,t,integrator)
            sqrt(u[7]*u[7] + u[8]*u[8]) > 1.4
    end

    function stop_integration(integrator)
        println("Terminated at t = $(integrator.t)")
        terminate!(integrator)
        #necessary, otherwise PowerGridSolution throws error
        integrator.sol = DiffEqBase.solution_new_retcode(integrator.sol, :Success)
    end

    cb1 = DiscreteCallback(voltage_deadband, timer_off)
    cb2 = DiscreteCallback(voltage_outside, timer_on)
    cb3 = DiscreteCallback(timer_hit, TapState)
    cb4 = DiscreteCallback(((u,t,integrator) -> t in 10.0), newGridImpedanze)
    cb5 = DiscreteCallback(check_voltage_low, stop_integration)
    cb6 = DiscreteCallback(check_voltage_high, stop_integration)

    params = GetParametersOLTCHisken(Tp)
    prob = ODEProblem(rhs(pg), ic, tspanOLTCHisken(),[params;-1])
    integrator = init(
    #sol = solve(
        prob,
        Rodas4(),
        callback = CallbackSet(cb1, cb2, cb3,cb4,cb5),
        dt = 1e-2,
        adaptive = false,
        tstops=[10.0],
        maxiters = 1e5,
        progress = true,
        save_everystep=false,
    )
end


function call_integrator(x)
    #_prob = remake(prob,u0=x[1:end-4],p=x[end-3:end])
    set_u!(integrator,ForwardDiff.value.(x)[1:end-5])
    integrator.p = ForwardDiff.value.(x)[end-4:end]
    step!(integrator)
    return integrator.sol[end]
end


function make_Jac(integrator, ic,tmax)
    Jac = Matrix{Float64}[]
    usol   = Vector{Float64}[]
    p =  integrator.p
    u = integrator.u
    result = DiffResults.JacobianResult(u,[u;p])

    while integrator.t <= tmax
        ForwardDiff.jacobian!(result,call_integrator,[u;p])
        push!(Jac,DiffResults.jacobian(result))
        push!(usol,DiffResults.value(result))
        u = usol[end]
        p = integrator.p
        display(integrator.t)
    end
    return integrator.sol,usol,Jac
end

sol, u1, J1 = make_Jac(integrator,ic,30.0)
pgsol = PowerGridSolution()
plot(sol,vars=(5))
xlims!((0,30))

tmp1 = rand(1500)
for i=1:1500
    tmp1[i] = J1[i][5,5]
end
plot(tmp1)

u -> rhs(pg)(zero(u),u,params,0.0)
f = rhs(pg)
tmp =  f(zero(ic0),ic0,params,0.0)
dx = ForwardDiff.gradient(u -> Subsystem1(zero([1,1]),[1.0,1.0],[1.0,1.0],0.0),[1.0,1.0])
