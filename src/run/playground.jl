using DifferentialEquations
using Plots

function Subsystem1(du,u,p,t)
    du[1] = -u[1] + p[1]*u[2]
    du[2] = p[2]*u[1] - u[2]
end

function Subsystem2(du,u,p,t)
    du[1] = -u[1] + 100*u[2]
    du[2] = -10*u[1] - u[2]
end

function condition1(u,t,integrator) # Event when event_f(u,t) == 0
    u[1]*u[2] >= 0.0
end

function condition2(u,t,integrator) # Event when event_f(u,t) == 0
    u[1]*u[2] < 0.0
end

function affect1(integrator)
    integrator.p[1] = 10.0
    integrator.p[2] = -100.0
end

function affect2(integrator)
    integrator.p[1] = 100.0
    integrator.p[2] = -10.0
end

cb1 = DiscreteCallback(condition1,affect1)
cb2 = DiscreteCallback(condition2,affect2)

u0 = [0.7;0.7]
#u0 = [0.7;0.0]
tspan = (0.0,1.0)
p = [10.0,-100.0]
prob = ODEProblem(Subsystem1,u0,tspan,p)
sol = solve(prob,Rodas4(),dtmax=1e-3, callback = CallbackSet(cb1,cb2))
plot(sol,vars =(1,2))

fig = plot(sol,vars=(1,2),size = (400, 300), legend = false,
            framestyle = :box,
            linewidth = 1.5,
            xlabel = "\$x_1\$",
            labelfontsize = 15,
            #xlims = (-2.5,2.5),
            xlims = (-0.7,0.7),
            xtickfont = font(8, "Akkurat"),
            ylabel = "\$x_2\$",
            #ylims = (-2.5,2.5),
            ylims = (-0.7,0.7),
            #linecolor = RGB(146/255,208/255,80/255),
            #linecolor = RGB(244/255,177/255,131/255),
            ytickfont = font(8, "Akkurat"),
            grid = true,
            gridalpha = 0.5,
            gridstyle = :dash
            )

savefig(fig,"C:\\Users\\liemann\\Desktop\\Präsi_DFG\\switched_stable.svg")
