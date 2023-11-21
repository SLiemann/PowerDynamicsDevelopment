using ModelingToolkit, DifferentialEquations, Plots

function VoltageStep(;irr=false)
    @variables t
    @connector function Pin(;name)
        sts = @variables v(t)=1.0 i(t)=1.0 [connect = Flow] 
        ODESystem(Equation[], t, sts, []; name=name)
    end

    function Ground(;name)
        @named g = Pin()
        eqs = [g.v ~ 0]
        compose(ODESystem(eqs, t, [], []; name=name), g)
    end

    function OnePort(;name,irv = false, iri = false)
        @named p = Pin()
        @named n = Pin()
        @variables v(t)=1.0 [irreducible=irv] 
        @variables i(t)=1.0 [connect = Flow, irreducible=iri]
        eqs = [
            v ~ p.v - n.v
            0 ~ p.i + n.i
            i ~ p.i
            ]
        compose(ODESystem(eqs, t, [v,i], []; name=name), p, n)
    end

    function Resistor(;name, R = 1.0)
        @named oneport = OnePort()
        @unpack v, i = oneport
        ps = @parameters R=R
        eqs = [
            v ~ i * R
            ]
        extend(ODESystem(eqs, t, [], ps; name=name), oneport)
    end

    function Capacitor(;name, C = 1.0)
        @named oneport = OnePort()
        @unpack v, i = oneport
        ps = @parameters C=C
        D = Differential(t)
        eqs = [
            D(v) ~ i / C
            ]
        extend(ODESystem(eqs, t, [], ps; name=name), oneport)
    end

    function ACStepVoltage(;name, V = 1.0, freq = 1.0, phase = 0.0,irreducible = false)
        @named oneport = OnePort(irv=irreducible)
        @unpack v = oneport
        ps = @parameters V=V freq=freq phase=phase
        eqs = [
                v ~ V*sin(2*pi*freq*t+phase)
            ]
        extend(ODESystem(eqs, t, [], ps; name=name), oneport)
    end

    @named AC = ACStepVoltage(V=230*sqrt(2),freq = 50.0,irreducible=irr) #voltage state is either observed or not
    @named R = Resistor(R=10)
    @named C = Capacitor(C=1e-3)
    @named ground = Ground()

    rc_eqs = [
        connect(AC.p,R.p)
        connect(R.n,C.p)
        connect(C.n,ground.g,AC.n)
    ]

    step = [0.1] => [AC.V ~ 0.3*230*sqrt(2)]

    @named _rc_model = ODESystem(rc_eqs, t,[],[]; systems = [AC,R,C,ground],discrete_events = [step])
    sys = structural_simplify(_rc_model)

    u0 = zeros(length(states(sys)))

    prob = ODEProblem(sys, u0, (0.0, 0.2));
    sol = solve(prob,Rodas4(autodiff=false),dtmax = 1e-4,tstops=[0.04]); 
    plot!(sol,vars=[AC.v,R.v,C.v],layout=(1,3))
end

VoltageStep(irr=false) #not correct
VoltageStep(irr=true)  #correct