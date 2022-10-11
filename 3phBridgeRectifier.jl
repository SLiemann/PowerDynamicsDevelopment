using Plots
include("C://Users//liemann//github//PowerDynamicsDevelopment//Acausal_models.jl")
#function Get3phBridgeRectifier()
Ubase = 230.0; #Volt
Pbase = 3000.0; #Watt
Zbase = Ubase^2 / Pbase;

@named Phase_a = ACStepVoltage(V=230 * sqrt(2), freq=50.0)
@named Phase_b = ACStepVoltage(V=230 * sqrt(2), freq=50.0, phase=-2 / 3 * pi)
@named Phase_c = ACStepVoltage(V=230 * sqrt(2), freq=50.0, phase=+2 / 3 * pi)
@named Rt1 = Resistor(R=0.08083)#.082
@named Rt2 = Resistor(R=0.08083)#0.079
@named Rt3 = Resistor(R=0.08083)#0.08
@named Lt1 = Inductance(L=0.223e-3)#0.221e-3
@named Lt2 = Inductance(L=0.223e-3)#0.226e-3)
@named Lt3 = Inductance(L=0.223e-3)#0.223e-3)

#@named D1 = DiodeModelica()
#@named D2 = DiodeModelica()
#@named D3 = DiodeModelica()
#@named D4 = DiodeModelica()
#@named D5 = DiodeModelica()
#@named D6 = DiodeModelica()

#@named D1 = ShockleyDiode()
#@named D2 = ShockleyDiode()
#@named D3 = ShockleyDiode()
#@named D4 = ShockleyDiode()
#@named D5 = ShockleyDiode()
#@named D6 = ShockleyDiode()

@named D1 = Diode(uf=0.8)
@named D2 = Diode(uf=0.8 - 0.1)
@named D3 = Diode(uf=0.8 + 0.1)
@named D4 = Diode(uf=0.8 - 0.2)
@named D5 = Diode(uf=0.8 + 0.2)
@named D6 = Diode(uf=0.8 - 0.3)

Ldv = 0.0125 / (100 * pi) * Zbase; #L = x_pu*Zbase/(2*50*ω0)
@named Ld = Inductance(L=Ldv)
Cdv = 1.0 / (0.361 * Zbase) / (100 * pi); #C = 1.0/(x_pu*Zbase)/(2*50*ω0)
@named Cd = Capacitor(C=Cdv)
@named CPL = ConstPower(P=3000.0)
@named ground = Ground()

rc_eqs = [
        connect(Phase_a.p, Rt1.p)
        connect(Phase_b.p, Rt2.p)
        connect(Phase_c.p, Rt3.p)connect(Rt1.n, Lt1.p)
        connect(Rt2.n, Lt2.p)
        connect(Rt3.n, Lt3.p)connect(Lt1.n, D1.p, D6.n)
        connect(Lt2.n, D2.p, D5.n)
        connect(Lt3.n, D3.p, D4.n)connect(D1.n, D2.n, D3.n, Ld.p)
        connect(Ld.n, Cd.p, CPL.p)
        connect(D6.p, D5.p, D4.p, CPL.n, Cd.n)connect(Phase_a.n, Phase_b.n, Phase_c.n, ground.g)]
ssa_start = [0.06] => [Phase_a.dV ~ 0.5]
ssb_start = [0.06] => [Phase_b.dV ~ 0.5]
ssc_start = [0.06] => [Phase_c.dV ~ 0.5]

ssa_end = [0.10] => [Phase_a.dV ~ 0.5]
ssb_end = [0.10] => [Phase_b.dV ~ 0.5]
ssc_end = [0.10] => [Phase_c.dV ~ 0.5]

@named _rc_model = ODESystem(rc_eqs, t, [], [];
        systems=[Phase_a, Phase_b, Phase_c, Rt1, Rt2, Rt3, Lt1, Lt2, Lt3, D1, D2, D3, D4, D5, D6, Ld, Cd, CPL, ground], discrete_events=[ssa_start, ssb_start, ssc_start])#,ssa_end,ssb_end,ssc_end
sys = structural_simplify(_rc_model)
ModelingToolkit.calculate_jacobian(sys)
u0_shock = [
        Ld.i => 9.45806283188738 * 0
        Cd.v => 542.976018704617
        Lt2.p.i => -9.596686636791368 * 0
        Lt3.p.i => 9.59668868756968 * 0
        Lt1.v => 0.00683609415944265 * 0
        D1.v => -277.25609798283506 * 0
        D2.v => -554.5012294190553 * 0
        D3.v => 9.54123856462757 * 0
        D4.v => -554.5012294211061 * 0
        D5.v => 9.541236513849462 * 0
        D6.v => -277.23559019970645 * 0
        CPL.i => 5.525105891706098
        D(Lt1.p.i) => 5.525105891706098 * 0
]
u0_mod = [
        Ld.i => 9.45806283188738 * 0
        Cd.v => 542.976018704617
        Lt2.p.i => -9.596686636791368 * 0
        Lt3.p.i => 9.59668868756968 * 0
        Lt1.v => 0.00683609415944265 * 0
        D1.s => -277.25609798283506 * 0
        D2.s => -554.5012294190553 * 0
        D3.s => 9.54123856462757 * 0
        D4.s => -554.5012294211061 * 0
        D5.s => 9.541236513849462 * 0
        D6.s => -277.23559019970645 * 0
        CPL.i => 5.525105891706098
        D(Lt1.p.i) => 5.525105891706098 * 0
]
u0_normal = [
        Ld.i => 9.45806283188738 * 0
        Cd.v => 545.0
        Lt2.p.i => -9.596686636791368 * 0
        Lt3.p.i => 9.59668868756968 * 0
        Phase_a.v => 230 * sqrt(2) * sin(100pi * 0)
        Phase_b.v => 230 * sqrt(2) * sin(100pi * 0 - 2 / 3 * pi)
        Phase_c.v => 230 * sqrt(2) * sin(100pi * 0 + 2 / 3 * pi)
        Lt1.v => 0.00683609415944265 * 0
        D1.v => -277.25609798283506 * 0
        D2.v => -554.5012294190553 * 0
        D3.v => 9.54123856462757 * 0
        D4.v => -554.5012294211061 * 0
        D5.v => 9.541236513849462 * 0
        D6.v => -277.23559019970645 * 0
        CPL.i => 3000 / 545.0
        D(Lt1.p.i) => 5.525105891706098 * 0
]
u0_normal_02 = [
        Ld.i => 9.550944953265073
        Cd.v => 541.0423908312898
        Lt2.p.i => -9.827611531668929
        Lt3.p.i => 9.827621734252663
        Phase_a.v => 1.5144965274442842e-12
        Phase_b.v => -281.69132042006663
        Phase_c.v => 281.69132042006396
        Lt1.v => 0.013671844908952688
        D1.v => -276.33231916116046
        D2.v => -553.5307928173809
        D3.v => 0.9071675557566079
        D4.v => -553.5207929594064
        D5.v => 0.9171674137310546
        D6.v => -276.2813062424894
        CPL.i => 5.544852031632163
        D(Lt1.p.i) => 61.308721564810455
]

#=u0 = [
  9.45806283188738
  542.976018704617
  -9.596686636791368  
  9.59668868756968
  0.00683609415944265 
  -277.25609798283506
  -554.5012294190553 
  9.54123856462757 
  -554.5012294211061 
  9.541236513849462
  -277.23559019970645
  5.525105891706098
  ] =#
u0_start = [
        9.45806283188738
        542.976018704617
        -9.596686636791368
        9.59668868756968
        #0
        #230*sqrt(2)*sin(-2/3*pi)
        #230*sqrt(2)*sin(+2/3*pi)
        0.00683609415944265
        277.25609798283506
        -554.5012294190553
        9.54123856462757
        -554.5012294211061
        9.541236513849462
        -277.23559019970645
        5.525105891706098
        0
]


#tmp = Get3phBridgeRectifier()
prob = ODEProblem(sys, u0_normal_02, (0.0, 0.2))
sol = solve(prob, Rodas4(autodiff=false), alg_hints=:stiff, dtmax=1e-6, maxiters=1e6, force_dtmin=true, tstops=[0.06, 0.1])
#,abstol=1e-12,reltol=1e-12
plot(sol, vars=[Phase_a.i, Phase_b.i, Phase_c.i], legend=:left)
#,ylims=(-0.5,0.5),xlims=(0.02,0.06))
plot(sol, vars=[Phase_a.i], ylims=(-17, 17), xlims=(0.02, 0.06))
plot(sol, vars=[Lt1.v, Lt2.v, Lt3.v], legend=:false, xlims=(0.05995, 0.06005))
plot(sol, vars=[Cd.v])

plot(sol, vars=[Phase_a.v, Phase_b.v, Phase_c.v], legend=:left)
plot(sol, vars=[CPL.i])
plot(sol, vars=[D(Lt1.p.i), Lt1.p.i], layout=(2, 1))
plot(sol, vars=[CPL.i])





u0_tmp = deepcopy(Float64.(sol.u[end]))
for i in ModelingToolkit.default_p(sys)
        display(i)
end
ModelingToolkit.discrete_events(sys)
ModelingToolkit.continuous_events(sys)
