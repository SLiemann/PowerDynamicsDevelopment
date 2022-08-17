using Plots
include("C://Users//liemann//github//PowerDynamicsDevelopment//Acausal_models.jl")
#function Get3phBridgeRectifier()
Ubase = 230.0; #Volt
Pbase = 3000.0; #Watt
Zbase = Ubase^2/Pbase;

@named Phase_a = ACStepVoltage(V=230*sqrt(2),freq = 50.0, dV = 0.8)
@named Phase_b = ACStepVoltage(V=230*sqrt(2),freq = 50.0, dV = 0.8, phase =-2/3*pi)
@named Phase_c = ACStepVoltage(V=230*sqrt(2),freq = 50.0, dV = 0.8, phase =+2/3*pi)
@named Rt1 = Resistor(R= 0.08083)#.082
@named Rt2 = Resistor(R= 0.08083)#0.079
@named Rt3 = Resistor(R= 0.08083)#0.08
@named Lt1 = Inductance(L= 0.223e-3)#0.221e-3
@named Lt2 = Inductance(L= 0.223e-3)#0.226e-3)
@named Lt3 = Inductance(L= 0.223e-3)#0.223e-3)

@named D1 = DiodeModelica()
@named D2 = DiodeModelica()
@named D3 = DiodeModelica()
@named D4 = DiodeModelica()
@named D5 = DiodeModelica()
@named D6 = DiodeModelica()

#@named D1 = Diode()
#@named D2 = Diode()
#@named D3 = Diode()
#@named D4 = Diode()
#@named D5 = Diode()
#@named D6 = Diode()

Ldv = 0.0125/(100*pi)*Zbase; #L = x_pu*Zbase/(2*50*ω0)
@named Ld = Inductance(L= Ldv)
Cdv = 1.0/(0.361*Zbase)/(100*pi); #C = 1.0/(x_pu*Zbase)/(2*50*ω0)
@named Cd = Capacitor(C= Cdv)
@named CPL = ConstPower(P = 3000.0)
@named ground = Ground()

rc_eqs = [
        connect(Phase_a.p,Rt1.p)
        connect(Phase_b.p,Rt2.p)
        connect(Phase_c.p,Rt3.p)

        connect(Rt1.n,Lt1.p)
        connect(Rt2.n,Lt2.p)
        connect(Rt3.n,Lt3.p)

        connect(Lt1.n,D1.p,D6.n)
        connect(Lt2.n,D2.p,D5.n)
        connect(Lt3.n,D3.p,D4.n)

        connect(D1.n,D2.n,D3.n,Ld.p)
        connect(Ld.n,Cd.p,CPL.p)
        connect(D6.p,D5.p,D4.p,CPL.n,Cd.n)

        connect(Phase_a.n,Phase_b.n,Phase_c.n,ground.g)

        ]
ssa_start = [0.0400] => [Phase_a.dV ~ 0.6]
ssb_start = [0.0401] => [Phase_b.dV ~ 0.6]    
ssc_start = [0.0402] => [Phase_c.dV ~ 0.6]    

#ssa_end =   [0.10] => [Phase_a.dV ~ 1.0]  
#ssb_end =   [0.10] => [Phase_b.dV ~ 1.0]  
#ssc_end =   [0.10] => [Phase_c.dV ~ 1.0]  

@named _rc_model = ODESystem(rc_eqs, t,[],[];
                    systems = [Phase_a,Phase_b,Phase_c,Rt1,Rt2,Rt3,Lt1,Lt2,Lt3,D1,D2,D3,D4,D5,D6,Ld,Cd,CPL,ground])
                    #discrete_events = [ssa_start,ssb_start,ssc_start])#,ssa_end,ssb_end,ssc_end])
sys = structural_simplify(_rc_model)

#=u0 = [
      Ld.i  => 9.45806283188738
      Cd.v  => 542.976018704617
      Lt2.i => -9.596686636791368  
      Lt3.i => 9.59668868756968
      Lt1.v => 0.00683609415944265 
      D1.s  => -277.25609798283506
      D2.s  => -554.5012294190553 
      D3.s  => 9.54123856462757 
      D4.s  => -554.5012294211061 
      D5.s  => 9.541236513849462
      D6.s  => -277.23559019970645
      CPL.i => 5.525105891706098
      ]=#
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
                0
                230*sqrt(2)*sin(-2/3*pi)
                230*sqrt(2)*sin(+2/3*pi)
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
prob = ODEProblem(sys, u0_tmp, (0.0, 0.1))
sol = solve(prob,Rodas4(autodiff=false),alg_hints=:stiff,dtmax = 1e-6,maxiters = 1e6,force_dtmin=true,
        abstol=1e-14,reltol=1e-14) #,tstops = [0.04,0.401,0.402]

plot!(sol,vars = [Phase_a.i,Phase_b.i,Phase_c.i])#,ylims=(-0.5,0.5),xlims=(0.02,0.06))
plot(sol,vars = [Phase_a.i],ylims=(-15,15),xlims=(0.02,0.06))
plot(sol,vars = [Lt3.i])
plot(sol,vars = [Cd.v])

plot(sol,vars = [Phase_a.v,Phase_b.v,Phase_c.v])
plot(sol,vars = [Lt1.i])


u0_tmp = deepcopy(Float64.(sol.u[end]))
ModelingToolkit.default_p(sys)
ModelingToolkit.discrete_events(sys)
ModelingToolkit.continuous_events(sys)