using ModelingToolkit, Plots, DifferentialEquations
using IfElse


@variables t
@connector function Pin(;name)
    sts = @variables v(t)=1.0 i(t)=1.0 [connect = Flow] #[irreducible=true]
    ODESystem(Equation[], t, sts, []; name=name)
end

function Ground(;name)
    @named g = Pin()
    eqs = [g.v ~ 0]
    compose(ODESystem(eqs, t, [], []; name=name), g)
end

function OnePort(;name)
    @named p = Pin()
    @named n = Pin()
    @variables v(t)=1.0  
    @variables i(t)=1.0 [connect = Flow]
    eqs = [
        v ~ p.v - n.v
        0 ~ p.i + n.i
        i ~ p.i
        ]
    compose(ODESystem(eqs, t, [v,i], []; name=name), p, n)
end

function ACVoltage(;name, V = 1.0)
    @named oneport = OnePort()
    @unpack v = oneport
    ps = @parameters V=V 
    eqs = [
        v ~ V*cos(2*pi*50*t)
        ]
    extend(ODESystem(eqs, t, [], ps; name=name), oneport)
end

affect!(integ,u,v,ctx) = integ.p[2] = 10

function Resistor1(;name, R = 1.0)
    @named oneport = OnePort()
    @unpack v, i = oneport
    ps = @parameters R=R
    eqs = [
        v ~ i * R
        ]
    extend(ODESystem(eqs, t, [], ps;discrete_events =[v<0.5 => (affect!, [], [R], nothing)], name=name), oneport)
end

function Resistor2(;name, R = 2.0)
    @named p = Pin()
    @named n = Pin()
    @variables v(t)=1.0 [irreducible =true]
    @variables i(t)=1.0 [connect = Flow]
    ps = @parameters R=R
    eqs = [
        v ~ p.v - n.v
        0 ~ p.i + n.i
        i ~ p.i
        v ~ i * R
        ]
    compose(ODESystem(eqs, t, [], ps; discrete_events  = [v<0.5 => (affect!, [], [R=>:r], nothing)], name=name), n,p)
end

@named Grid = ACVoltage()
#@named Rload = Resistor1(R=0.1) #ERROR: MethodError: no method matching iterate(::Base.Iterators.Reverse{Term{Bool, Nothing}})
@named Rload = Resistor2() 
@named ground = Ground()

rc_eqs = [
        connect(Grid.p, Rload.p)
        connect(Rload.n, Grid.n)
        connect(Grid.n, ground.g)
        ]
@named _rc_model = ODESystem(rc_eqs, t,[],[];systems = [Rload,Grid,ground])

sys = structural_simplify(_rc_model)

prob = ODEProblem(sys, [1.0],(0, 0.02))
sol = solve(prob,Rodas4(),dtmax = 1e-6)
plot(sol, vars=[Rload.v]) # for Resistor2: ERROR: ArgumentError: System Rload: variable v does not exist
plot(sol, vars=[Grid.v]) #works
