using PowerDynamics, MAT
using OrderedCollections: OrderedDict
using Distributed, DataFrames
@everywhere using IfElse
import DiffEqBase: initialize_dae!
using PlotlyJS

include("C:/Users/liemann/github/PowerDynamicsDevelopment/src/include_custom_nodes_lines_utilities.jl")

Sbase = 100e6
Ubase = 110e3
Ibase = Sbase/Ubase/sqrt(3)
Zbase = Ubase^2/Sbase

buses=OrderedDict(
    "bus0" => SlackAlgebraic(U=1.0),# 1.05717
    "bus1" => VoltageDependentLoad(P=-6.0, Q=0.0, U=1.0, A=1.0, B=0.0,Y_n = complex(0.0)))
branches=OrderedDict(
    "OLTC"=> StaticPowerTransformerTapParam(from="bus0",to="bus1",Sbase=Sbase,Srated=100e6,uk=0.10,XR_ratio=Inf,
                i0=0.0,Pv0=0.0,tap_side = "LV",tap_pos = 0,tap_inc = 2.0,tap_max=20,tap_min=-20,p_ind=1))
pg = PowerGrid(buses, branches);

U1,δ1,ic0,cu = PowerFlowClassic(pg,iwamoto = false,max_tol = 1e-3,Ustart = [1.0, 0.85],δstart=[0,-30/180*pi])
pg, ic = InitializeInternalDynamics(pg,ic0)

pgsol = sim_tap(pg,ic,(0.0,20.0));
plot(plotallvoltages(pgsol))

dir = "\\\\fs0\\home\\liemann\\"
file = matopen(dir*"mattap.mat")
vltvs = read(file, "vmat");
close(file)

file = matopen(dir*"pftap.mat")
pftap = read(file, "vpf");
close(file)

p1 = plotallvoltages(pgsol);
append!(p1,[scatter(x=vltvs[:,1],y=vltvs[:,2],name="Matlab")])
append!(p1,[scatter(x=pftap[:,1],y=pftap[:,2],name="PowerFactory")])
plot(p1)

function sim_tap(pg::PowerGrid,ic::Vector{Float64},tspan::Tuple{Float64,Float64})
    problem= ODEProblem{true}(rhs(pg),ic,tspan,[0])
    timer_start = -1.0
    tap_dir = 1

    branch_oltc = "OLTC"
    index_U_oltc = PowerDynamics.variable_index(pg.nodes,pg.lines[branch_oltc].to,:u_r)

    function TapState(integrator)
        timer_start = integrator.t
        sol1 = integrator.sol
        integrator.p[1] += 1*tap_dir 
        initialize_dae!(integrator,BrownFullBasicInit())
        auto_dt_reset!(integrator)
    end

    function voltage_deadband(u,t,integrator)
         0.99 <= sqrt(u[index_U_oltc]*u[index_U_oltc] + u[index_U_oltc+1]*u[index_U_oltc+1]) <= 1.01
    end

    function timer_off(integrator)
        if timer_start != -1
            timer_start = -1
        end
    end

    function voltage_outside_low(u,t,integrator)
         sqrt(u[index_U_oltc]*u[index_U_oltc] + u[index_U_oltc+1]*u[index_U_oltc+1]) < 0.99
    end

    function voltage_outside_high(u,t,integrator)
        sqrt(u[index_U_oltc]*u[index_U_oltc] + u[index_U_oltc+1]*u[index_U_oltc+1]) > 1.01
   end

    function timer_on_low(integrator)
        tap_dir = 1
        if timer_start == -1
            timer_start = integrator.t
        end
    end

    function timer_on_high(integrator)
        tap_dir = -1
        if timer_start == -1
            timer_start = integrator.t
        end
    end

    function timer_hit(u,t,integrator)
        if timer_start == -1
            return false
        else
            return t-timer_start > 1.0
        end
    end

    function check_OLTC_voltage(u,t,integrator)
            sqrt(u[index_U_load]*u[index_U_load] + u[index_U_load+1]*u[index_U_load+1]) < 0.3 && t > 5.0
    end

    cb1 = DiscreteCallback(voltage_deadband, timer_off)
    cb2 = DiscreteCallback(voltage_outside_low, timer_on_low)
    cb3 = DiscreteCallback(voltage_outside_high, timer_on_high)
    cb4 = DiscreteCallback(timer_hit, TapState)

    stiff  = repeat([:stiff],length(ic))
    sol = solve(problem, Rodas5(autodiff=true), callback = CallbackSet(cb1,cb2,cb4,cb3), dtmax = 1e-3,force_dtmin=false,maxiters=1e6, initializealg = BrownFullBasicInit(),alg_hints=:stiff,abstol=1e-8,reltol=1e-8) 

    return PowerGridSolution(sol, pg)
end
