using ModelingToolkit, Plots, DifferentialEquations
using IfElse
using DiffEqSensitivity

begin
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

    function Diode(;name, Uf = 0.8,Ron = 0.001, Roff = 9999)
        @named oneport = OnePort(irv = true)
        @unpack v, i = oneport
        ps = @parameters Ron=Ron Roff= Roff Uf=Uf 
        eqs = [
            i ~ IfElse.ifelse(v<=Uf, v/Roff ,v/Ron + Uf*(1.0/Roff - 1.0/Ron))
            ]
        extend(ODESystem(eqs, t, [], ps;continuous_events  = [[v-Uf~0.0] =>(affect!, [v], [], nothing)], name=name), oneport)
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

    function DiodeModelica(;name, Uf = 0.8,Ron = 1e-5, Roff = 1e5)
        @named oneport = OnePort()#irv = true,iri = true
        @unpack v, i = oneport
        sts = @variables s(t)=0.0
        ps = @parameters Ron=Ron Roff= Roff Uf=Uf 
        eqs = [
            v ~ s*IfElse.ifelse(s<0, 1,Ron) #+ Uf
            i ~ s*IfElse.ifelse(s<0, 1/Roff,1) #+ Uf/Roff
            ]
        extend(ODESystem(eqs, t, sts, ps;continuous_events  = [s~0], name=name), oneport)#continuous_events  = [[v-Uf~0.0],[i~0.0]]
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

    Rparam= 0.311
    @named Rg = Resistor(R= Rparam)
    @named Rload = Resistor(R= (230^2*2)/230)
    @named Lg = Inductance(L= 24.89183309e-3)
    @named Ld = Inductance(L= 0.024876)
    #@named D1 = Diode(Uf = 0.8,Ron=0.01,Roff = 1e4)
    #@named D2 = Diode(Uf = 0.8,Ron=0.01,Roff = 1e4)
    #@named D3 = Diode(Uf = 0.8,Ron=0.01,Roff = 1e4)
    #@named D4 = Diode(Uf = 0.8,Ron=0.01,Roff = 1e4)
    @named D1 = DiodeModelica()
    @named D2 = DiodeModelica()
    @named D3 = DiodeModelica()
    @named D4 = DiodeModelica()
    #@named D1 = IdealDiode()
    #@named D2 = IdealDiode()
    #@named D3 = IdealDiode()
    #@named D4 = IdealDiode()
    #@named D1 = ShockleyDiode()
    #@named D2 = ShockleyDiode()
    #@named D3 = ShockleyDiode()
    #@named D4 = ShockleyDiode()
    @named Cd = Capacitor(C= 0.2306593e-3)#2.3066e-04
    @named Grid = ACVoltage(V=230*sqrt(2),freq = 50.0)
    @named ground = Ground()
    @named CPL = ConstPower(P = 230.0)

    rc_eqs = [
            connect(Grid.p, Rg.p)
            connect(Rg.n, Lg.p)
            connect(Lg.n, D1.p)
            connect(D2.n, D1.n)
            connect(D2.p, D3.n)
            connect(Cd.p, D2.n)
            connect(D3.p, D4.p)
            connect(CPL.p, Cd.p)
            connect(CPL.n, D4.p)
            connect(Cd.n, D4.p)
            connect(D1.p, D4.n)
            connect(Grid.n, D3.n)
            connect(Grid.n, ground.g)
            ]
    @named _rc_model = ODESystem(rc_eqs, t,[],[];
                        systems = [Rg,Lg,D1,D2,D3,D4,Cd,Grid,ground,CPL])

    sys = structural_simplify(_rc_model)

    #=su0 = [
    #Lg.i => 0.0
    Cd.v => 230*sqrt(2)
    D1.v => -230*sqrt(2)/2
    D2.v => -230*sqrt(2)/2
    D3.v => -230*sqrt(2)/2
    D4.v => -230*sqrt(2)/2
    #Cd.i => 230*sqrt(2)/1000
    #D1.i => 0.0
    #D2.i => 0.0
    #D3.i => 0.0
    #D4.i => 0.0
    ]
    
    rc_eqs = [
            connect(Grid.p, D1.n)
            connect(D1.p, Rload.p)
            connect(Rload.n, Grid.n)
            connect(Grid.n, ground.g)
            ]
    @named _rc_model = ODESystem(rc_eqs, t,[],[];systems = [D1,Grid,ground,Rload])
    sys = structural_simplify(_rc_model) =#
end
begin    
    #u0 = [0.0]
    #u0= [Lg.i =>0,Cd.v =>230*sqrt(2),D1.s =>0,D2.s =>0,D3.s =>0,D4.s =>0]
    u0= [0,230*sqrt(2),0,0,0,0,sqrt(2)]
    prob = ODEProblem(sys, u0, (0.0, 0.1))
    #prob = ODEForwardSensitivityProblem(sys, u0, (0.0, 0.1),werte,ForwardDiffSensitivity())
    sol = solve(prob,Rodas4(autodiff=false),dtmax = 1e-6,maxiters = 1e6,force_dtmin=true)# ,dtmax = 1e-6 ,dtmin=1e-6,
    sol.t[end]
end

x,dp = extract_local_sensitivities(sol)
plot(sol.t,dp[5][1:1,:]')

plot(sol,vars = [D1.v,D1.i])
plot(sol,vars = [D1.v])
plot(sol,vars = [D1.i])
plot(sol,vars = [D1.s])
plot(sol,vars = [Rload.i])
plot(sol,vars = [D1.i,D2.i,D3.i,D4.i], layout = (2,2),linewidth = 1.5)
plot(sol,vars = [D1.v,D2.v,D3.v,D4.v], layout = (2,2),linewidth = 2,legend=:left)

plot(sol,vars = [Ds1.v,Ds4.v], layout = (2,2),linewidth = 2)
plot(sol,vars = [Grid.v])
plot(sol,vars = [Grid.i])
plot(sol,vars = [Cd.v])
plot(sol,vars = [Lg.i])
plot(sol)


st_p = collect(keys(ModelingToolkit.get_defaults(sys)))
indizes = setdiff(indexin(ModelingToolkit.get_ps(sys),st_p),[nothing])
werte = Float64.(collect(values(ModelingToolkit.get_defaults(sys)))[indizes])

irms = RMS(sol[Grid.i],sol.t)
urms = RMS(sol[Grid.v],sol.t)
plot(sol.t,irms.*urms)
plot(sol.t,-sol[Grid.i].*sol[Grid.v])

ModelingToolkit.continuous_events(sys)
ModelingToolkit.discrete_events(sys)


ylims!((-0.5,1))

@named ac = ACVoltage(V=230*sqrt(2))
@named C = Capacitor(C=0.01)
@named D1 = Diode(Uf = 0.8,Ron=0.1,Roff = 10000)
@named R = Resistor(R=10)
@named gr = Ground()
rc_eqs = [
        connect(ac.p, D1.n)
        connect(D1.p, R.p)
        connect(R.n, ac.n)
        connect(ac.n, gr.g)
        ]
    
@named _rc_model = ODESystem(rc_eqs, t)
@named rc_model = compose(_rc_model,[R,D1,ac,gr])
sys = structural_simplify(rc_model)

prob = ODEProblem(sys, [0.0], (0, 0.1))
sol = solve(prob, Rodas4(), dtmax =1e-6, progress = true,force_dtmin=true)
plot(sol,vars = [ac.i,R.v], layout=(2,1))





function math_diode(U; Uf = 0.8,Roff = 1e4,Ron = 0.001)
    if U <= Uf
        return (U)/Roff;
    else
        return (U)/Ron + Uf*(1.0/Roff - 1.0/Ron)
    end
end

function math_diode2(U; Uf = 0.8,Roff = 1e4,Ron = 0.001)
    if U <= Uf
        return (U-Uf)/Roff;
    else
        return (U-Uf)/Ron #+ Uf*(1.0/Roff - 1.0/Ron)
    end
end

u = collect(0.7999:0.00001:0.8001)
plot(u,math_diode.(u),legend=false)

function shockley(U;Is = 1e-8,n=1.5,Ut=25e-3)
    Steigung = 100
    Uknick = n*Ut*(log(Steigung*n*Ut/Is)+2)

    if U <= Uknick
        return Is * (exp(U/(n*Ut))-1.0)
    else
        return (U-Uknick)*Steigung + Is * (exp(Uknick/(n*Ut))-1.0)
    end
end
U = 0.0:1e-4:1

plot(U,shockley.(U),legend=false)

plot(U,shockley1.(U))
plot!(U,shockley2.(U))
ylims!(-0.1,0.1)
plot(U,math_diode2.(U))

Is = 1e-8
n=1.5
Ut=25e-3

x = n*Ut*(log(100*n*Ut/Is)+2)



@parameters t
@variables u(t)


D = Differential(t)

eqs = [D(u) ~ -u]

affect1!(integ, u, p, ctx) = integ.u[u.r] = 0

# set condition to be an equation instead of a timestop
@named sys = ODESystem(eqs, t, [u], [], discrete_events = [u < 1.0=> (affect1!, [u => :r], [], nothing)])
prob = ODEProblem(sys, [u => 10.0], (0, 10.0)) 
sol = solve(prob, Tsit5(); tstops = [4.0])
plot(sol)