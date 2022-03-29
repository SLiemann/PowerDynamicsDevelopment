using ModelingToolkit, Plots, DifferentialEquations

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

function OnePort(;name)
    @named p = Pin()
    @named n = Pin()
    sts = @variables v(t)=1.0 i(t)=1.0
    eqs = [
           v ~ p.v - n.v
           0 ~ p.i + n.i
           i ~ p.i
          ]
    compose(ODESystem(eqs, t, sts, []; name=name), p, n)
end

function Resistor(;name, R = 1.0)
    @named oneport = OnePort()
    @unpack v, i = oneport
    ps = @parameters R=R
    eqs = [
           v ~ i * R
          ]
    extend(ODESystem(eqs, t, [], ps; name=name  ), oneport)
end

function Capacitor(;name, C = 1.0)
    @named oneport = OnePort()
    @unpack v, i = oneport
    ps = @parameters C=C
    D = Differential(t)
    eqs = [
           D(v) ~ i / C
          ]
    extend(ODESystem(eqs, t, [], ps; name=name  ), oneport)
end

function ACVoltage(;name, V = 1.0, freq = 1.0, phase = 0.0)
    @named oneport = OnePort()
    @unpack v,i = oneport
    ps = @parameters V=V freq=freq phase=phase
    eqs = [
           V*sin(2*pi*freq*t+phase) ~ v
          ]
    extend(ODESystem(eqs, t, [], ps; name=name ), oneport)
end

@named R = Resistor(R= 1)
@named C = Capacitor(C= 1e-3)
@named source = ACVoltage(V=1*sqrt(2),freq = 50.0)
@named ground = Ground()

 rc_eqs = [
            connect(source.p,R.p)
            connect(R.n,C.p)
            connect(C.n,source.n)
            connect(source.n,ground.g)
          ]


@named rc_model = ODESystem(rc_eqs, t,[R.v],[];
                            systems = [R,C,source,ground],
                            continuous_events = [C.v ~  1.0] => [source.freq ~ 1.0])

sys = structural_simplify(rc_model)

u0 = [
    C.v => 0.0
  ]
prob = ODEProblem(sys, u0, (0,2.0))

sol = solve(prob, Rodas4())
plot(sol,vars = [R.i], layout=(1,1))

ModelingToolkit.continuous_events(sys)


function dq(a,b,c,θ)
    θ += pi/2
    faktor = sqrt(2/3)
    d = faktor * (sin(θ)*a + sin(θ-2*pi/3)*b + sin(θ-4*pi/3)*c)
    q = faktor * (cos(θ)*a + cos(θ-2*pi/3)*b + cos(θ-4*pi/3)*c)
    z = faktor * (sqrt(0.5)*  (a+b+c))

    return d,q,z
end

theta = 0:1e-3:6*pi
ua = sqrt(1)*cos.(theta)
ub = sqrt(1)*cos.(theta.-2*pi/3)
uc = sqrt(1)*cos.(theta.+2*pi/3)

ia = 0.5*sqrt(2)*sin.(theta)
ib = 0.5*sqrt(2)*sin.(theta.-2*pi/3)
ic = 0.5*sqrt(2)*sin.(theta.+2*pi/3)


dqz = dq.(ua,ub,uc,theta)
idqz = dq.(ia,ib,ic,theta)

pdqz= dqz.*idqz

dqz = [collect(i) for i in dqz]
for i in enumerate(dqz)
    dqz

function tup2mat(tup)
    u = zeros(length(tup))
    v = zeros(length(tup))
    w = zeros(length(tup))

    for (ind,val) in enumerate(tup)
        u[ind] = val[1]
        v[ind] = val[2]
        w[ind] = val[3]
    end
    return [u v w]
end
function invdq(d,q,z,θ)
    θ += pi/2
    faktor = sqrt(2/3)
    a = faktor * (sin(θ)*d + cos(θ)*q + sqrt(0.5)*z)
    b = faktor * (sin(θ-2*pi/3)*d + cos(θ-2*pi/3)*q + sqrt(0.5)*z)
    c = faktor * (sin(θ-4*pi/3)*d + cos(θ-4*pi/3)*q + sqrt(0.5)*z)

    return a,b,c
end

abc = invdq.(d,q,z,theta)

tmp = tup2mat(abc)

plot(tmp)
plot!(ua)

t = 0:1e-4:0.04
ua = sqrt(1)*cos.(2*pi*50 .*t)
ub = sqrt(1)*cos.(2*pi*50 .*t .-2*pi/3)
uc = sqrt(1)*cos.(2*pi*50 .*t .+2*pi/3)

plot(t,ua)
plot!(t,ub)
plot!(t,uc)
