using ModelingToolkit, DifferentialEquations
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
    extend(ODESystem(eqs, t, [], ps; name=name), oneport)
end

function ACStepVoltage(;name, V = 1.0, freq = 1.0, phase = 0.0, dV = 1.0)
    @named oneport = OnePort(irv=false)
    @unpack v = oneport
    ps = @parameters V=V freq=freq phase=phase dV = dV switch=1
    eqs = [
        # v~IfElse.ifelse(switch<0,dV*V*sin(2*pi*freq*t+phase),V*sin(2*pi*freq*t+phase))
        v~ dV*V*sin(2*pi*freq*t+phase)
        ]
    extend(ODESystem(eqs, t, [], ps; name=name), oneport)
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

function affect!(integ,u,v,ctx)
    display("------------")
    display(integ.t)
    display(integ.u[u.v])
end

function Diode(;name, Uf = 0.8,Ron = 0.01, Roff = 0.5e4)
    @named oneport = OnePort(irv = true)
    @unpack v, i = oneport
    ps = @parameters Ron=Ron Roff= Roff Uf=Uf 
    eqs = [
        i ~ IfElse.ifelse(v<=Uf, v/Roff ,v/Ron + Uf*(1.0/Roff - 1.0/Ron))
        ]
    extend(ODESystem(eqs, t, [], ps;continuous_events  = [v-Uf~0.0], name=name), oneport) # =>(affect!, [v], [], nothing)]
    #extend(ODESystem(eqs, t, [], ps;continuous_events  = [[v~Uf],[v~Uf-0.01],[v~Uf+0.01]], name=name), oneport)
end

affectID!(integ,u,v,ctx) = integ.u[u.i] = 0

function IdealDiode(;name)
    @named p = Pin()
    @named n = Pin()
    @variables v(t)=1.0 [irreducible=true] 
    @variables i(t)=1.0 [connect = Flow, irreducible=true]
    eqs = [
        v ~ p.v - n.v
        0 ~ p.i + n.i
        i ~ p.i
        #0 ~ IfElse.ifelse(v<0, i ,v)
        i ~ IfElse.ifelse(v<=0, v/1e5 ,v/1e-3)
        ]
    #compose(ODESystem(eqs, t, [v,i], [];discrete_events  = [v<0.0 => (affectID!, [i], [], nothing)], name=name),p,n)
    compose(ODESystem(eqs, t, [v,i], [];continuous_events  = [v~0.0], name=name),p,n)
end

function DiodeModelica(;name,Ron = 1e-3, Roff = 0.5e4)
    @named oneport = OnePort()#irv = true,iri = true
    @unpack v, i = oneport
    sts = @variables s(t)=0.0
    ps = @parameters Ron=Ron Roff= Roff
    eqs = [
        v ~ s*IfElse.ifelse(s<0, 1,Ron) #+ Uf
        i ~ s*IfElse.ifelse(s<0, 1/Roff,1) #+ Uf/Roff
        ]
    extend(ODESystem(eqs, t, sts, ps;continuous_events  = [s~0], name=name), oneport)
end

function ShockleyDiode(;name,Is = 10e-8,n=1.5,Ut=25e-3, Steigung = 100)
    @named oneport = OnePort(irv=true)
    @unpack v, i = oneport
    ps = @parameters Is=Is n=n Ut=Ut
    Uknick = n*Ut*(log(Steigung*n*Ut/Is)+2)

    eqs = [
        i ~ IfElse.ifelse(v <= Uknick, Is * (exp(v/(n*Ut))-1.0),(v-Uknick)*Steigung + Is * (exp(Uknick/(n*Ut))-1.0)) 
        ]
    extend(ODESystem(eqs, t, [], ps;continuous_events  = [v~Uknick], name=name), oneport)
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
