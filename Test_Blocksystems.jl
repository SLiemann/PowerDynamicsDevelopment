using Pkg
Pkg.activate(joinpath(@__DIR__, "../"))
using ModelingToolkit
using BlockSystems
using OrdinaryDiffEq
using Plots

@parameters t M D P_m(t) P_e(t)
@variables ω(t)
dt = Differential(t)

swing = IOBlock([dt(ω) ~ 1/M * (P_m - D*ω - P_e)], # equations
                 [P_m, P_e],                       # inputs
                 [ω])                              # outputs

swing = set_p(swing, :D=>0.5, :M=>1)

# system without pid controller, fixed pref
@variables P_m(t)
pfix = IOBlock([P_m ~ 1], [], [P_m])

wo_pid = IOSystem([pfix.P_m => swing.P_m], [pfix, swing]) |> connect_system

# system with pid controller for p_ref
@parameters input(t)
@variables int(t) out(t) pid(t)
pid = IOBlock([dt(int) ~ input,
               pid ~ input + int + dt(input),
               out ~ 1 - pid],
              [input],
              [out])

w_pid = IOSystem([pid.out => swing.P_m,
                  swing.ω => pid.input],
                 [swing, pid];
                 namespace_map=[pid.out => :P_m]
                 )
w_pid = connect_system(w_pid)

# connect step function
function simulate(sys)
    @variables P_e(t)
    step = IOBlock([P_e ~ 1 - 0.25/(1 + exp(-30*(t)))], [], [P_e])

    sys = IOSystem([step.P_e => sys.P_e], [step, sys]) |> connect_system

    gen = generate_io_function(sys)

    u0 = zeros(length(gen.states))
    odef = ODEFunction((du, u, p, t) -> gen.f_ip(du, u, nothing, p, t); mass_matrix=gen.massm, syms=Symbol.(gen.states))
    prob = ODEProblem(odef, u0, (-0.5, 7))

    sol = solve(prob, Rodas4())
end

tmp = generate_io_function(w_pid)
tmp.states

###################
function get_Res(R::Float64)
    @parameters t ϕ1(t) ϕ2(t)
    @variables i(t)
    Res = IOBlock([i ~ (ϕ1 - ϕ2)/R],[ϕ1,ϕ2],[i])
end

R1 = get_Res(5.0)

function get_Ind(L::Float64)
    @parameters t ϕ1(t) ϕ2(t)
    @variables i(t)
    dt = Differential(t)
    Ind = IOBlock([dt(i) ~ (ϕ1 - ϕ2)/L],[ϕ1,ϕ2],[i])
end
L1 = get_Ind(3.0)

function get_AC(Ud::Float64, freq::Float64)
    @parameters t
    @variables ϕ1(t) ϕ2(t)
    AC = IOBlock([ϕ1 ~ Ud*sin(freq*2*pi*t), ϕ2 ~ 0.0],[],[ϕ1,ϕ2])
end

AC1 = get_AC(230.0,50.0)

sys = IOSystem([AC1.ϕ1 => R1.ϕ1, AC1.ϕ2 => R1.ϕ2], [AC1,R1]) |> connect_system
equations(sys.system)



@parameters t L ϕ1(t) ϕ2(t) ϕ0
@variables i(t)
dt = Differential(t)
Res = IOBlock([dt(i) ~ (ϕ1 - ϕ2)/L],[ϕ1,ϕ2],[i])
AC  = IOBlock([])
