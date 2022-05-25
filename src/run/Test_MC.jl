using PowerDynamics
using DifferentialEquations
using Plots
using OrderedCollections
using MAT

begin
    include("C:/Users/liemann/github/PowerDynamicsDevelopment/src/include_custom_nodes_lines_utilities.jl")
    include("C:/Users/liemann/github/PowerDynamicsDevelopment/src/operationpoint/InitializeInternals.jl")
end
function getMachtingGrid(;y_new = 0.0)
    Ubase = 1e3;
    Sbase = 100e6;
    Zbase = Ubase^2/Sbase;

    Zbase_dc = (3*1e3*sqrt(2/3))^2/Sbase

    Rdc = 1.2
    Rdc_pu = 1.2/Zbase_dc
    Gdc = 1.0/Rdc_pu
    Xcdc = 1.0/(100*pi*0.008*200) /Zbase_dc
    Cdc = 1.0/Xcdc/100/pi#(100*pi*Xcdc)

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
        "bus3" => dVOC(
            Sbase = 100e6,
            Srated = 100e6,
            p0set = 0.9, # based on Sbase!
            q0set = 0.0,
            u0set = 1.00,
            eta = pi*Zbase*2/3,
            alpha = 0.1*2/3*1000^2,
            Kdc = 100.0, #1600
            gdc = Gdc,
            cdc = Cdc,
            xlf = Xlf,    #  0.15
            rf = R_f, # 0.0005
            xcf =  Xcf ,# 15.1515151515
            Tdc = 0.05,
            Kp_u = 0.52,#*200, #1.0
            Ki_u = 1.161022,#*200.0, #1.161022#
            Kp_i = 0.738891,#, /200.0  301.6510
            Ki_i = 1.19,#/200.0 , 485.815
            imax_csa = 1.25,
            p_red = 1.0,
            ϵ = 1e-4,
            p_ind = collect(1:16)
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

        "branch2"=> StaticPowerTransformer(from="bus2",to="bus3",Sbase=Sbase,Srated=100e6,uk=0.1500833,XR_ratio=30.0,
                                           i0=0.0,Pv0=0.0,tap_side = "LV",tap_pos = 0,tap_inc = 1.0))
    pg = PowerGrid(buses, branches)
end
function simMatching(prob_sim)
    pg_new = getMachtingGrid(y_new = 1444)#1/0.035/2
    pg_pre = getMachtingGrid()
    tstep = [0.1,0.45]
    function fault_state(integrator)
        new_f = rhs(pg_new)
        op_prob = ODEProblem(new_f, integrator.sol[end], (0.0, 1e-6),integrator.p, initializealg = BrownFullBasicInit())
        ic_new = solve(op_prob,Rodas5())
        integrator.f = new_f
        integrator.cache.tf.f = integrator.f
        integrator.u = ic_new.u[end]
    end

    function postfault_state(integrator)
        sol = integrator.sol
        ic_tmp = deepcopy(integrator.sol.u[indexin(tstep[1],integrator.sol.t)[1]])
        #ic_tmp = getPreFaultVoltages(pg_new,ic_tmp,deepcopy(sol[end]))
        ic_tmp = getPreFaultAlgebraicStates(pg_pre,ic_tmp,deepcopy(sol[end]))
        #display(ic_tmp .- sol.u[end])
        op_prob = ODEProblem(rhs(pg_pre), ic_tmp, (0.0, 1e-6),integrator.p, initializealg = BrownFullBasicInit())
        ic_new = solve(op_prob,Rodas4())
        #display(ic_tmp .- ic_new)
        #display(rhs(pg_pre).syms)
        integrator.f = rhs(pg_pre)
        integrator.cache.tf.f = integrator.f
        integrator.u = ic_new.u[end]
    end

    cb = DiscreteCallback(((u,t,integrator) -> t in tstep[1]), fault_state)
    cb1 = DiscreteCallback(((u,t,integrator) -> t in tstep[2]), postfault_state)

    sol = solve(prob, Rodas4(), tstops= tstep,callback = CallbackSet(cb,cb1), dtmax = 1e-3,progress=true, initializealg = BrownFullBasicInit())
    #sol.retcode = :Success
    #return sol, event_recorder
    return PowerGridSolution(sol,pg_new)
end

begin
    include("C:/Users/liemann/github/PowerDynamicsDevelopment/src/include_custom_nodes_lines_utilities.jl")
    pg = getMachtingGrid()
    U,δ,ic0 = PowerFlowClassic(pg, iwamoto = true, max_tol = 1e-7)
    pg1 ,ic = InitializeInternalDynamics(pg,ic0)
    params = getallParameters(pg1.nodes["bus3"])[6:21]
    prob = ODEProblem(rhs(pg1),ic,(0.0,1.0),params)#initializealg = BrownFullBasicInit()
    #sol = solve(prob, Rodas4(),dtmax = 1e-3,progress=true)
    #pgsol = PowerGridSolution(sol,pg1)
    #pgsol = simMatchingNEW(pg1,(0.0,0.3),ic)
    pgsol = simMatching(prob)
end

plot(pgsol,vars=(8))

plot(pgsol.dqsol,vars=(7:13))
plot(pgsol.dqsol,vars=(9))
plot(pgsol,["bus3"],:p)
plot(pgsol,["bus3"],:q)
plot(pgsol,["bus3"],:iabs)
plot(pgsol,["bus3"],:i_abs)

  #,ylim=(0.91,0.915)
plot(pgsol,["bus3"],:v)
plot(pgsol,["bus3"],:Pf)
plot(pgsol,["bus3"],:vd)

plot(pgsol,"bus3",:udc) #,xlim=(0.9,1.6)
plot(pgsol,"bus3",:idc0)
plot(pgsol,"bus3",:θ)
plot!(pgsol,"bus3",:idc0)
plot(pgsol,"bus3",:Pdelta)
plot(pgsol,"bus3",:e_ud)
plot(pgsol,"bus3",:id0)
plot!(pgsol,"bus3",:iq0)
plot!(pgsol,"bus3",:ids)
plot!(pgsol,"bus3",:iqs)
plot(pgsol,"bus3",:x_uabs)

file = matopen("C:\\Users\\liemann\\Desktop\\simtime.mat")
t = read(file, "simtime")'
close(file)

file = matopen("C:\\Users\\liemann\\Desktop\\Vconv.mat")
v = read(file, "Vconv")
close(file)


file = matopen("C:\\Users\\liemann\\Desktop\\PQ.mat")
PQ = read(file, "PQ")
close(file)
plot(pgsol,["bus3"],:p)
plot!(t',PQ[:,1],label = "MATLAB")
xlims!((0,1))
ylims!((-0.1,1.3))

plot(pgsol,["bus3"],:q)
plot!(t',PQ[:,2],label = "MATLAB")
xlims!((0,1))
ylims!((-1.5,1))

plot(pgsol,["bus3"],:v)
plot!(t',v,label = "MATLAB")
xlims!((0.05,1))
ylims!((0.0,1.1))

ic_try = deepcopy(pgsol.dqsol.u[end])
pg_try = deepcopy(pg1)
#pg_try = getMachtingGrid(y_new = 1/0.035/4)

params = deepcopy(getallParameters(pg.nodes["bus3"])[5:18])
prob = ODEProblem(rhs(pg_try),ic_try,(0.0,0.1),params)

integrator = init(prob,Rodas4(), initializealg = BrownFullBasicInit())
step!(integrator)

sol = solve(prob, Rodas4(), initializealg = BrownFullBasicInit())
pgsol_try = PowerGridSolution(sol,pg_try)
plot(pgsol_try,"bus3",:v)
#
integrator2 = reinit!(integrator)


function simMatchingNEW(pg_pre,timespan,x1)
    pg_new = getMachtingGrid(y_new = 1/0.035/4)
    #pg_pre = getMachtingGrid()
    params = getallParameters(pg_pre.nodes["bus3"])[5:18]

    tstep = [0.1,0.25]

    regular = rhs(pg_pre)
    error = rhs(pg_new)

    if regular.mass_matrix != error.mass_matrix || length(regular.syms) != length(error.syms)
        error("Change of MassMatrix or system size in abstract pertubation not supported!")
    end

    # wrap f and introduce parameter: if p=true no error, if p=false errorstate
    _f = (dx, x, p, t) -> Int(p[1])==1 ? regular(dx,x,p,t) : error(dx,x,p,t) #[2:end]

    f = ODEFunction(_f, mass_matrix = regular.mass_matrix, syms = regular.syms)

    problem = ODEProblem{true}(f, x1, timespan, [1;params])

    function errorState(integrator)
        display("Error")
        display(integrator.p[1])
        integrator.p[1] = 2
        if integrator.opts.adaptive
            auto_dt_reset!(integrator)
            set_proposed_dt!(integrator, integrator.dt)
        end

        # reset the adaptive timestepping
        sol = integrator.sol
        #ic_tmp = deepcopy(integrator.sol.u[indexin(tstep[1],integrator.sol.t)[1]])
        #ic_tmp = getPreFaultAlgebraicStates(pg_pre,ic_tmp,deepcopy(sol[end]))
        op_prob = ODEProblem(error, sol.u[end], (0.0, 1e-6),integrator.p, initializealg = BrownFullBasicInit())
        ic_new = solve(op_prob,Rodas4())
        integrator.u = ic_new.u[end]

    end

    function regularState(integrator)
        display("Regular")
        display(integrator.p[1])
        integrator.p[2] = 1
        if integrator.opts.adaptive
            auto_dt_reset!(integrator)
            set_proposed_dt!(integrator, integrator.dt)
        end
        #reinit!(integrator, initializealg = BrownFullBasicInit())
        # reset the adaptive timestepping
        sol = integrator.sol
        ic_tmp = deepcopy(integrator.sol.u[indexin(tstep[1],integrator.sol.t)[1]])
        ic_tmp = getPreFaultVoltages(pg_pre,ic_tmp,deepcopy(sol.u[end]))
        #ic_tmp = getPreFaultAlgebraicStates(pg_pre,ic_tmp,deepcopy(sol[end]))
        op_prob = ODEProblem(regular, sol.u[end], (0.0, 1e-3),integrator.p, initializealg = BrownFullBasicInit())
        ic_new = solve(op_prob,Rodas4())
        integrator.u = ic_new.u[end]
        #display(integrator.u)
    end

    cb1 = PresetTimeCallback([tstep[1]], errorState)
    cb2 = PresetTimeCallback([tstep[2]], regularState) # tstops= tstep

    sol = solve(problem, Rodas4(),callback = CallbackSet(cb1,cb2), dtmax = 1e-3,progress=true, initializealg = BrownFullBasicInit())

    return sol
    #return PowerGridSolution(sol,pg_new)
end
