using PowerDynamics
using OrderedCollections: OrderedDict

function Testgrid_oPFC(;u_new = 1.0)
    buses = OrderedDict(
        "bus1" => SlackAlgebraic(U = u_new),
        "bus2" => oPFC(Cd = 1.0 /(100*pi*0.036)*2,
                      Pdc = -1.0,
                      Ulow = 0.3,
                      Qn= 0.21*0,
                      t0 = 0.0,
                      ϵ=0.01,
                      p_ind = collect(1:6),
                      ),
    )

    branches = OrderedDict(
        "branch1" => PiModelLine(
            from = "bus1",
            to = "bus2",
            y = 1.0/(0.00+1im*0.20),#1.97712-1im*19.78677,
            y_shunt_km = 0.0,
            y_shunt_mk = 0.0,),
        )
    pg = PowerGrid(buses, branches)
end

function oPFC_params()
    pg = Testgrid_oPFC()
    PFC = pg.nodes["bus2"]
    return [
        PFC.Cd,
        PFC.Pdc,
        PFC.Ulow,
        PFC.Qn,
        PFC.t0,
        PFC.ϵ,
    ]
end


tfault_on() = 0.1
tfault_off() = 0.2

function sim_oPFC(;u_new = 1.0)
    pg = Testgrid_oPFC()
    pg_fault = Testgrid_oPFC(u_new =u_new)
    tfault = [tfault_on(),tfault_off()]

    equal_grid(u,t,integrator) = abs(abs(u[3]+1im*u[4])-u[5]) <= integrator.p[6] && abs(u[3]+1im*u[4]) > integrator.p[3]
    function affect_equal(integrator)
       integrator.u[6] = 1.0
       integrator.u[5] = abs(integrator.u[3]+1im*integrator.u[4])
    end

    low_grid(u,t,integrator) = u[5] > abs(u[3]+1im*u[4]) + integrator.p[6] #&& u[1] > integrator.p[3]
    function affect_low(integrator)
      integrator.u[6] = 2.0
    end

    high_grid(u,t,integrator) = u[5] < abs(u[3]+1im*u[4]) - integrator.p[6] && abs(u[3]+1im*u[4])  > integrator.p[3] && integrator.u[3] != 3.0
    function affect_high(integrator)
      integrator.u[9] = 0.0
      integrator.u[6] = 3.0
      integrator.p[5] = asin(integrator.u[2]) / (2*pi*50)
    end

    disconnect(u,t,integrator) = abs(u[3]+1im*u[4])  <= integrator.p[3] && abs(u[3]+1im*u[4]) <= integrator.p[3]
    function affect_disconnect(integrator)
      integrator.u[6] = 4.0
    end

    cb1 = DiscreteCallback(equal_grid,affect_equal)
    cb2 = DiscreteCallback(low_grid,affect_low)
    cb3 = DiscreteCallback(high_grid,affect_high)
    cb4 = DiscreteCallback(disconnect,affect_disconnect)

    function affect_fault_on(integrator)
        new_f = rhs(pg_fault)
        op_prob = ODEProblem(new_f, integrator.sol[end], (0.0, 1e-6),integrator.p, initializealg = BrownFullBasicInit())
        ic_new = solve(op_prob,Rodas5())
        integrator.f = new_f
        integrator.cache.tf.f = integrator.f
        integrator.u = ic_new.u[end]
    end
    function affect_fault_off(integrator)
        new_f = rhs(pg)
        op_prob = ODEProblem(new_f, integrator.sol[end], (0.0, 1e-6),integrator.p, initializealg = BrownFullBasicInit())
        ic_new = solve(op_prob,Rodas5())
        integrator.f = new_f
        integrator.cache.tf.f = integrator.f
        integrator.u = ic_new.u[end]
    end

    cb5 = DiscreteCallback(((u,t,integrator) -> t in tfault[1]), affect_fault_on)
    cb6 = DiscreteCallback(((u,t,integrator) -> t in tfault[2]), affect_fault_off)


    U,δ,ic0 = PowerFlowClassic(pg, iwamoto = true, max_tol = 1e-7)
    pg1 ,ic = InitializeInternalDynamics(pg,ic0)
    params = oPFC_params()
    prob = ODEProblem(rhs(pg1),ic,(0.0,0.3),params, initializealg = BrownFullBasicInit())

    sol = solve(prob,dtmax = 1e-4,callback = CallbackSet(cb1,cb2,cb3,cb4,cb5,cb6), tstops=tfault)
    return PowerGridSolution(sol,pg)
end
