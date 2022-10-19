using PowerDynamics
using OrderedCollections: OrderedDict
using DifferentialEquations
using Plots

include("C:/Users/liemann/github/PowerDynamicsDevelopment/src/include_custom_nodes_lines_utilities.jl")
function getPG(;zf=0)
    zf = zf/(400e3^2/8000e6)
    Sbase = 8000e6
    buses=OrderedDict(
        "bus0" => SlackAlgebraic(U=1.0),
        "bus1" => VoltageDependentLoad(P=0.0, Q=0.0, U=1.0, A=1.0, B=0.0,Y_n = complex(0.0)),
        "bus_sm" => gentpj(Sbase=Sbase,Srated=5300e6,E_fd=1.0, H=6, P=4400e6/Sbase, D=0.0, Ω=50, R_a=0, T_d0s=7.0, T_q0s=1.5, T_d0ss=0.05, T_q0ss=0.05, X_d=2.2, X_q=2.0, X_ds=0.3, X_qs=0.4, X_dss=0.2, X_qss=0.2, X_l=0.15, S_10=0.1, S_12=0.3,K_is=0.0))

    
    Z = (0.01+1im*0.1)*400e3^2/5000e6
    Z = Z/(400e3^2/8000e6)                             
    branches= OrderedDict(
        "Line_0-1"=> PiModelLine(from= "bus0", to = "bus1",y=1.0/Z, y_shunt_km=0.0, y_shunt_mk=0.0),
        "Trafo_Netz"=> StaticPowerTransformer(from="bus1",to="bus_sm",Sbase=Sbase,Srated=5300e6,uk=0.15,XR_ratio=Inf,
                                           i0=0.0,Pv0=0.0,tap_side = "HV",tap_pos = 0,tap_inc = 1.0))
    pg=PowerGrid(buses, branches)
end
pg = getPG()
U1,δ1,ic0,cu = PowerFlowClassic(pg,iwamoto = false,max_tol = 1e-7)
display(U1.=> δ1)
pg, ic0 = InitializeInternalDynamics(pg,ic0)

#problem = ODEProblem(rhs(pg),ic0,(0,10))
#sol = solve(problem,Rodas4(),dtmax=1e-3)
#plot(sol)


pg_fault = deepcopy(pg)
zf = 10/(400e3^2/8000e6)
pg_fault.nodes["bus1"] = VoltageDependentLoad(P=0.0, Q=0.0, U=1.0, A=1.0, B=0.0,Y_n = complex(1.0/zf))

function errorState(integrator)
    sol1 = integrator.sol
    ode = rhs(pg_fault)
    op_prob = ODEProblem(ode, sol1[end], (0.0, 1e-6), initializealg = BrownFullBasicInit())
    x2 = solve(op_prob,Rodas5())
    x2 = x2.u[end]

    integrator.f = ode
    integrator.cache.tf.f = integrator.f
    integrator.u = x2
end

function regularState(integrator)
    sol = integrator.sol
    ode   = rhs(pg)
    ic_tmp = deepcopy(integrator.sol.u[indexin(5.0,integrator.sol.t)[1]])
    ic_tmp = getPreFaultVoltages(pg,ic_tmp,deepcopy(sol[end]))
    op_prob = ODEProblem(ode, ic_tmp, (0.0, 1e-6), initializealg = BrownFullBasicInit())
    x3 = solve(op_prob,Rodas5())
    x3 = x3.u[end]

    integrator.f = rhs(pg)
    integrator.cache.tf.f = integrator.f
    integrator.u = x3#sol2[end]
end

cb1 = DiscreteCallback(((u,t,integrator) -> t in [5]), errorState)
cb2 = DiscreteCallback(((u,t,integrator) -> t in [5.1]), regularState)

sol = solve(problem, Rodas4(), callback = CallbackSet(cb1,cb2), tstops=[5,5.1], dtmax = 1e-3,progress =true)
pgsol = PowerGridSolution(sol,pg)

plot(pgsol,["bus_sm"],:v)
plot(pgsol,["bus_sm"],:p)
plot(pgsol,["bus_sm"],:ω)

using CSV, DataFrames


upf = CSV.read("C:\\Users\\liemann\\Desktop\\u.csv", DataFrame)
um = CSV.read("C:\\Users\\liemann\\Desktop\\vgl dirty\\data.csv", DataFrame)
plot(pgsol,["bus_sm"],:v,label="julia")
plot!(upf[:,1],upf[:,2],label="PF")
xlims!((4.9,10))
#=
xlims!((4.9,5.3))

plot!(um[:,1],um[:,2],linewidth=1.5,label="MATLAB")

xlims!((4.9,10))
ylims!((0.95,1.01))



pq = CSV.read("C:\\Users\\liemann\\Desktop\\pq.csv", DataFrame)
plot(pgsol,["bus_sm"],:p,label="julia")
plot!(pq[:,1],pq[:,3]./0.8*0.53,label="PF")
xlims!((4.9,10))

plot(pgsol,["bus_sm"],:q,label="julia")
plot!(pq[:,1],pq[:,2]./0.8*0.53,label="PF")
xlims!((4.9,10))=#
pq = CSV.read("C:\\Users\\liemann\\Desktop\\pq.csv", DataFrame)
plot(pgsol,["bus_sm"],:p,label="julia")
plot!(pq[:,1],pq[:,3]./0.8*0.53,label="PF")
xlims!((4.9,10))
