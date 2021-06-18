using PowerDynamics
#using OrderedCollections: OrderedDict
using Plots
#import PowerDynamics: PiModel
using DifferentialEquations
using CSV #read PF DataFrames
using DataFrames #for CSV
#using Distributed

Ubase = 380e3
Sbase = 100e6
Zbase = (Ubase^2)/Sbase

begin
    include("C:/Users/liemann/github/PowerDynamicsDevelopment/src/include_costum_nodes_lines_utilities.jl")
    include("C:/Users/liemann/github/PowerDynamicsDevelopment/src/grids/LTVS_Test_System.jl")
end
pgsol  = run_LTVS_simulation(pg,ic0,(0.0,70.0))
plot(pgsol,collect(keys(pg.nodes))[1:end-1],:v,legend = (0.3, 0.3))
plot(pgsol,"bus4",:i_abs)
plot!(pgsol,"bus4",:iabs)

begin
    plot(pgsol,collect(keys(pg.nodes))[1:end-1], :v,size = (1000, 500))
    test = DataFrame(CSV.File("C:\\Users\\liemann\\Desktop\\Basisszenario\\data.csv"; header=false, delim=';', type=Float64))
    plot!(test.Column1,test.Column2,label = "PF-bus1")
    plot!(test.Column1,test.Column3,label = "PF-bus2")
    plot!(test.Column1,test.Column4,label = "PF-bus3")
    plot!(test.Column1,test.Column5,label = "PF-bus4")
end
ylims!((0.85,1.01))
xlims!((0,8.0))
begin
    plot(pgsol,"bus4", :iabs,size = (1000, 500))
    plot!(pgsol,"bus4", :i_abs)
    test = DataFrame(CSV.File("C:\\Users\\liemann\\Desktop\\Basisszenario\\data_c.csv"; header=false, delim=';', type=Float64))
    plot!(test.Column1,test.Column2,label = "I")
    plot!(test.Column1,test.Column3,label = "I0")
end

ylims!((5,9))
xlims!((0,8.0))


pg = GFC_LTVS_Test_System()
#Load Flow
Qmax   = [Inf, Inf, Inf,Inf, sqrt(1-0.95^2), Inf]
Qmin   = -Qmax
U,δ1,ic = PowerFlowClassic(pg,iwamoto = true, Qmax = Qmax, Qmin = Qmin, Qlimit_iter_check = 2,max_tol = 1e-6)
Ykk = NodalAdmittanceMatrice(pg)
Uc = U.*exp.(1im*δ1/180*pi)
I_c = Ykk*Uc
S = conj(Ykk*Uc).*Uc
pg, ic0 = InitializeInternalDynamics(pg,I_c,ic)


plot(pgsol,"bus3",:p,legend = false)
plot!([0.0,120.0],[6.8,6.8])


function run_LTVS_simulation(pg::PowerGrid,ic1::Array{Float64,1},tspan::Tuple{Float64,Float64})
    tfault = [2.0, 2.15]
    Zfault = 20.0
    pg_fault = deepcopy(pg)
    pg_fault.nodes["busv"] = VoltageDependentLoad(P=0.0, Q=0.0, U=1.0, A=0., B=0.,Y_n = complex(1.0/(Zfault/Zbase)))
    nodes_postfault = deepcopy(pg.nodes)
    branches_postfault = deepcopy(pg.lines)
    delete!(nodes_postfault,"busv")
    delete!(branches_postfault,"Line_1-v")
    delete!(branches_postfault,"Line_v-2")
    pg_postfault = PowerGrid(nodes_postfault,branches_postfault)

    params = GFC_LTVS_params()
    problem = ODEProblem{true}(rhs(pg),ic1,tspan,params)
    timer_start = -1.0
    timer_now   = 0.0
    branch_oltc = "branch4"
    tap = pg.lines[branch_oltc].tap_pos
    OLTC = deepcopy(pg.lines[branch_oltc])
    tap_max = 20.0
    postfault_state = false
    fault_state = false
    #index_U_oltc = getNodeVoltageSymbolPosition(pg,pg.lines[branch_oltc].to)
    index_U_oltc = PowerDynamics.variable_index(pg.nodes,pg.lines[branch_oltc].to,:u_r)
    index_U_load = PowerDynamics.variable_index(pg.nodes,"bus3",:u_r)
    function TapState(integrator)
        timer_start = integrator.t
        sol1 = integrator.sol
        if tap < tap_max
            tap += 1
            node = StaticPowerTransformer(from=OLTC.from,to=OLTC.to,Srated=OLTC.Srated,
                                          uk=OLTC.uk,XR_ratio=OLTC.XR_ratio,i0=OLTC.i0,
                                          Pv0=OLTC.Pv0,Sbase=OLTC.Sbase,
                                          tap_side = OLTC.tap_side, tap_pos = tap,tap_inc = OLTC.tap_inc)
            if postfault_state
                np_pg = deepcopy(pg_postfault)
            elseif fault_state
                np_pg = deepcopy(pg_fault)
            else
                np_pg = deepcopy(pg)
            end
            np_pg.lines[branch_oltc] = node
            ode = rhs(np_pg)
            op_prob = ODEProblem(ode, sol1[end], (0.0, 1e-6), params, initializealg = BrownFullBasicInit())
            x2 = solve(op_prob,Rodas4())
            x2 = x2.u[end]

            integrator.f = ode
            integrator.cache.tf.f = integrator.f
            integrator.u = x2#sol1[end]
        end
    end

    function voltage_deadband(u,t,integrator)
         0.99 <= sqrt(u[index_U_oltc]*u[index_U_oltc] + u[index_U_oltc+1]*u[index_U_oltc+1]) <= 1.01
    end

    function timer_off(integrator)
        timer_start = -1
    end

    function voltage_outside(u,t,integrator)
         sqrt(u[index_U_oltc]*u[index_U_oltc] + u[index_U_oltc+1]*u[index_U_oltc+1]) < 0.99
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

    function errorState(integrator)
        sol1 = integrator.sol
        ode = rhs(pg_fault)
        op_prob = ODEProblem(ode, sol1[end], (0.0, 1e-6), params, initializealg = BrownFullBasicInit())
        x2 = solve(op_prob,Rodas5())
        x2 = x2.u[end]

        integrator.f = ode
        integrator.cache.tf.f = integrator.f
        integrator.u = x2
        fault_state = true
    end

    function regularState(integrator)
        sol = integrator.sol
        ode   = rhs(pg_postfault)
        index = PowerDynamics.variable_index(pg.nodes,"busv",:u_r)

        ic_tmp = deepcopy(integrator.sol.u[indexin(tfault[1],integrator.sol.t)[1]])
        #ic_tmp = getPreFaultVoltages(pg,ic_tmp,deepcopy(sol[end]))
        deleteat!(ic_tmp,index:index+1)
        op_prob = ODEProblem(ode, ic_tmp, (0.0, 1e-6), params, initializealg = BrownFullBasicInit())
        x3 = solve(op_prob,Rodas5())
        x3 = x3.u[end]

        resize!(integrator,length(ic0)-2)
        integrator.f = rhs(pg_postfault)
        integrator.cache.tf.f = integrator.f
        integrator.u = x3#sol2[end]

        postfault_state = true
        fault_state = false
        index_U_oltc = PowerDynamics.variable_index(pg_postfault.nodes,pg_postfault.lines[branch_oltc].to,:u_r)
        index_U_load = PowerDynamics.variable_index(pg_postfault.nodes,"bus3",:u_r)
    end

    function check_voltage(u,t,integrator)
            sqrt(u[index_U_load]*u[index_U_load] + u[index_U_load+1]*u[index_U_load+1]) < 0.2
    end

    function stop_integration(integrator)
        println("Terminated at $(integrator.t)")
        terminate!(integrator)
        #necessary, otherwise PowerGridSolution throws error
        integrator.sol = DiffEqBase.solution_new_retcode(integrator.sol, :Success)
    end

    cb1 = DiscreteCallback(voltage_deadband, timer_off)
    cb2 = DiscreteCallback(voltage_outside, timer_on)
    cb3 = DiscreteCallback(timer_hit, TapState)
    cb4 = DiscreteCallback(((u,t,integrator) -> t in tfault[1]), errorState)
    cb5 = DiscreteCallback(((u,t,integrator) -> t in tfault[2]), regularState)
    cb6 = DiscreteCallback(check_voltage, stop_integration)

    sol = solve(problem, Rodas4(), callback = CallbackSet(cb1,cb2,cb3,cb4,cb5,cb6), tstops=[tfault[1],tfault[2]], dtmax = 1e-2,progress =true) #
    sol = AddNaNsIntoSolution(pg,pg_postfault,deepcopy(sol))

    return PowerGridSolution(sol, pg)
end

nodes_postfault = deepcopy(pg.nodes)
branches_postfault = deepcopy(pg.lines)
delete!(nodes_postfault,"busv")
delete!(branches_postfault,"Line_1-v")
delete!(branches_postfault,"Line_v-2")
pg_postfault = PowerGrid(nodes_postfault,branches_postfault)

begin
    ic2 = deepcopy(ic0)
    deleteat!(ic2,3:4)
    f_test = rhs(pg_postfault)
    ode_test = ODEProblem(f_test,ic2,(0.0,10.0))
    testsol = solve(ode_test,Rodas4(),dtmax = 1e-3)
    pgtestsol = PowerGridSolution(testsol, pg_postfault)
    plot(pgtestsol,collect(keys(pg_postfault.nodes)), :v,size = (1000, 500),legend = (0.6, 0.75))
end

begin
    plot(pgsol,collect(keys(pg_postfault.nodes)),:v,legend = false)
    test = DataFrame(CSV.File("C:\\Users\\liemann\\Desktop\\u_pf.csv"; header=false, delim=';', type=Float64))
    plot!(test.Column1,test.Column2,label = "PF-bus1")
    plot!(test.Column1,test.Column3,label = "PF-bus2")
    plot!(test.Column1,test.Column4,label = "PF-bus3")
    plot!(test.Column1,test.Column5,label = "PF-gen")
end
begin
    plot(pgsol,["bus4"], :ifd,size = (1000, 500),legend = false)
    test = DataFrame(CSV.File("C:\\Users\\liemann\\Desktop\\u_ifd.csv"; header=false, delim=';', type=Float64))
    plot!(test.Column1,test.Column2,label = "PF-ifd")
end

begin
    plot(pgsol,["bus4"], :timer,size = (1000, 500),legend = false)
    test = DataFrame(CSV.File("C:\\Users\\liemann\\Desktop\\u_timer.csv"; header=false, delim=';', type=Float64))
    plot!(test.Column1,test.Column2,label = "PF-ifd")
end
xlims!((0.0,10.0))
ylims!((-18.2,-17))
