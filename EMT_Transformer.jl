using ModelingToolkit, Plots, DifferentialEquations
using IfElse

@variables t
@connector function Pin(;name)
    sts = @variables v(t)=1.0 i(t)=1.0 [connect = Flow]
    ODESystem(Equation[], t, sts, []; name=name, continuous_events = [1~0])
end

function Ground(;name)
    @named g = Pin()
    eqs = [g.v ~ 0]
    compose(ODESystem(eqs, t, [], []; name=name, continuous_events = [1~0]), g)
end

function OnePort(;name)
    @named p = Pin()
    @named n = Pin()
    sts = @variables v(t)=1.0 i(t)=1.0
    eqs = [
           v ~ p.v - n.v
           0 ~ p.i + n.i
           i ~ p.i
          ]
    compose(ODESystem(eqs, t, sts, []; name=name, continuous_events = [1~0]), p, n)
end

function Resistor(;name, R = 1.0)
    @named oneport = OnePort()
    @unpack v, i = oneport
    ps = @parameters R=R
    eqs = [
           v ~ i * R
          ]
    extend(ODESystem(eqs, t, [], ps; name=name, continuous_events = [1~0]), oneport)
end

function Capacitor(;name, C = 1.0)
    @named oneport = OnePort()
    @unpack v, i = oneport
    ps = @parameters C=C
    D = Differential(t)
    eqs = [
           D(v) ~ i / C
          ]
    extend(ODESystem(eqs, t, [], ps; name=name, continuous_events = [1~0]), oneport)
end

function ConstantVoltage(;name, V = 1.0)
    @named oneport = OnePort()
    @unpack v = oneport
    ps = @parameters V=V
    eqs = [
           V ~ v
          ]
    extend(ODESystem(eqs, t, [], ps; name=name), oneport)
end

function ACVoltage(;name, V = 1.0, freq = 1.0, phase = 0.0)
    @named oneport = OnePort()
    @unpack v = oneport
    ps = @parameters V=V freq=freq phase=phase
    eqs = [
           V*cos(2*pi*freq*t+phase) ~ v
          ]
    extend(ODESystem(eqs, t, [], ps; name=name, continuous_events = [1~0]), oneport)
end

function Inductance(;name, L = 1.0)
    @named oneport = OnePort()
    @unpack v, i = oneport
    ps = @parameters L=L
    D = Differential(t)
    eqs = [
           D(i) ~ v / L
          ]
    extend(ODESystem(eqs, t, [], ps; name=name), oneport)
end

Sb = 1e6/3
Ub = 10e3/sqrt(3)
Zb = Ub^2/Sb
I0 = 0.05
Pfe = 10000
uk = 0.1
XR = 2.0
ü = 2.0

fak = 1.0

rfe = Sb/Pfe
Rfe = rfe*Zb/fak
zm = 1.0/I0
Xm = 1.0/sqrt(1/(zm^2)-1/(rfe^2))
lh = Xm/(2*pi*50)/fak
ra = uk / sqrt(1 + XR^2)
xa = XR * ra
r1 = 0.5*ra/fak
l1 = 0.5*xa/(2*pi*50)/fak
R1 = 0.5*ra*Zb/fak
L1 = 0.5*xa*Zb/(2*pi*50)/fak
r2 = r1*(ü^2)
l2 = l1*(ü^2)
R2 = R1*(ü^2)
L2 = L1*(ü^2)

@named R1 = Resistor(R= r1)
@named R2 = Resistor(R= r2)
@named L1 = Inductance(L= l1)
@named L2 = Inductance(L= l2)
@named Rfe = Resistor(R= rfe)
@named Lh = Inductance(L= lh)
@named Grid = ACVoltage(V=1.0*sqrt(2),freq = 50.0)
@named load = Resistor(R=1.0)
@named ground = Ground()


 rc_eqs = [
           connect(Grid.p, R1.p)
           connect(R1.n, L1.p)
           connect(L1.n, Rfe.p)
           connect(Rfe.p, Lh.p)
           connect(Rfe.n, ground.g)
           connect(Lh.p, L2.p)
           connect(L2.n, R2.p)
           connect(R2.n, load.p)
           connect(load.n, ground.g)
           connect(Lh.n, load.n)
           connect(Grid.n, ground.g)
           ]

@named _rc_model = ODESystem(rc_eqs, t)
@named rc_model = compose(_rc_model,
                          [R1,R2,L1,L2,Lh,Rfe,Grid,ground,load])

#ModelingToolkit.continuous_events(sys)

sys = structural_simplify(rc_model)
u0 = [
   L1.i => 0.0
   Lh.i => 0.0
   L2.i => 0.0
  ]
prob = ODEProblem(sys, u0, (0, 0.1))

sol = solve(prob, Rodas4(), dtmax =1e-6, progress = true)
plot(sol,vars = [load.v,Grid.v], layout=(1,1))
plot(sol,vars = [Grid.i,load.i], layout=(1,1))
xlims!(0.96,1.0)
ylims!(1.24,1.25)
ylims!(1.,sqrt(2))


iload = sol[load.i]
uload = sol[load.v]
time = sol.t
plot(time,RMS(uload,time))
ylims!(0.8,1)
