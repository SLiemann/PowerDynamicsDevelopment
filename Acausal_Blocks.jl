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
           V*sin(2*pi*freq*t+phase) ~ v
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

function Diode(;name, Uf = 0.8,Ron = 0.001, Roff = 9999)
    @named oneport = OnePort()
    @unpack v, i = oneport
    ps = @parameters Ron=Ron Roff= Roff Uf=Uf # + Uf*(1.0/Roff - 1.0/Ron)
    cond_uf = [v ~ Uf]
    affect = [v~v]

    eqs = [
           i ~ IfElse.ifelse(v >= Uf, (v-Uf)/Ron + Uf/Roff,v/Roff)
          ]
    extend(ODESystem(eqs, t, [v], ps; name=name, continuous_events =[1~0]), oneport) # continuous_events = [v ~ Uf] => [v~v]
end

function ShockleyDiode(;name,Is = 10e-8,n=1.5,Ut=25e-3)
    @named oneport = OnePort()
    @unpack v, i = oneport
    ps = @parameters Is=Is n=n Ut=Ut
    eqs = [
           i ~ Is * exp(v/(n*Ut)-1.0)
          ]
    extend(ODESystem(eqs, t, [], ps; name=name), oneport)
end

function ConstPower(;name, P = 1.0)
    @named oneport = OnePort()
    @unpack v, i = oneport
    ps = @parameters P=P
    eqs = [
           P ~ v * i
          ]
    extend(ODESystem(eqs, t, [], ps; name=name), oneport)
end

@named Rg = Resistor(R= 1.0)
@named Rload = Resistor(R= (230*sqrt(2))^2/100)
@named Lg = Inductance(L= 4.817e-4)
@named Ld = Inductance(L= 0.024876)
@named D1 = Diode(Uf = 0.6,Ron=0.01,Roff = 1e8)
@named D2 = Diode(Uf = 0.7,Ron=0.01,Roff = 1e8)
@named D3 = Diode(Uf = 0.8,Ron=0.01,Roff = 1e8)
@named D4 = Diode(Uf = 0.9,Ron=0.01,Roff = 1e8)
@named Cd = Capacitor(C= 2.3066e-04)
@named Grid = ACVoltage(V=230*sqrt(2),freq = 50.0)
@named ground = Ground()
@named CPL = ConstPower(P = 230.0)

 rc_eqs = [
           connect(Grid.p, Rg.p)
           connect(Rg.n, D1.p)
           connect(D3.n, D1.p)
           connect(D1.n, D2.n)
           connect(Cd.p, D2.n)
           connect(D2.p, D4.n)
           connect(Rload.p, Cd.p)
           connect(Rload.n, D4.p)
           connect(Cd.n, D4.p)
           connect(D3.p, D4.p)
           connect(Grid.n, D4.n)
           connect(Grid.n, ground.g)
          ]

# => [D1.v~D1.v] => [D2.v~D2.v] => [D3.v~D3.v] => [D4.v~D4.v]
ce = [
    [D1.v ~ D1.Uf]
    [D2.v ~ D2.Uf]
    [D3.v ~ D3.Uf]
    [D4.v ~ D4.Uf]
]


@named _rc_model = ODESystem(rc_eqs, t,[D2.v,D3.v,D4.v],[];
                            systems = [Rg,D1,D2,D3,D4,Cd,Grid,ground,Rload],
                            continuous_events = ce)

@named _rc_model = ODESystem(rc_eqs, t,[D2.v,D3.v,D4.v])
@named rc_model = compose(_rc_model,
                          [Rg,D1,D2,D3,D4,Cd,Grid,ground,Rload])

#ModelingToolkit.continuous_events(sys)

sys = structural_simplify(_rc_model, simplify = false)
u0 = [
  # Lg.i => 0.0
   #Ld.i => 0.0
   Cd.v => 230*sqrt(2)
   D2.v => -230*sqrt(2)*0
   D1.v => -230*sqrt(2)*0
   D4.v => -0.0
   #Rload.i => 100/(230*sqrt(2))
  ]
prob = ODEProblem(sys, u0, (0, 0.1))
cond1(u,t,int) =

sol = solve(prob, Rodas4(), dtmax =1e-6, progress = true)
plot(sol,vars = [Rg.i], layout=(1,1))
ylims!((-0.5,1))


@named C = Capacitor(C=0.01)
@named D1 = Diode(Uf = 0.8,Ron=0.1,Roff = 10000)
#@named D1 = ShockleyDiode()
rc_eqs = [
          connect(Grid.p, D1.p)
          connect(D1.n, Rg.p)
          connect(Rg.n, Grid.n)
          connect(Grid.n, ground.g)
         ]

@named _rc_model = ODESystem(rc_eqs, t,continuous_events = [D1.v ~ D1.Uf])
@named rc_model = compose(_rc_model,
                         [Rg,D1,Grid,ground])

sys = structural_simplify(rc_model)
u0 = [
   #Ld.i => 0.0
   #Rg.i => 283/100
   #C.v => 283
   D1.v => 0.0
   ]
prob = ODEProblem(sys, u0, (0, 0.1))
sol = solve(prob, Rodas4(), dtmax =1e-6, progress = true)
plot(sol,vars = [D1.i], layout=(2,1))

ModelingToolkit.continuous_events(sys)






function math_diode(U; Uf = 0.08,Roff = 0.2,Ron = 0.1)
    It = Uf/Roff
    b2 = It - Uf/Ron
    if U <= Uf
        return U/Roff;
    else
        return U/Ron + Uf*(1.0/Roff - 1.0/Ron)
    end
end

function shockley(U;Is = 1e-8,n=1.5,Ut=25e-3)
    I = Is * exp(U/(n*Ut)-1.0)
end
U = 0:1e-3:0.8

plot!(U,shockley.(U))

plot(U,math_diode.(U))
