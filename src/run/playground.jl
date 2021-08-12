using DifferentialEquations
using Plots
using DataFrames

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

plot()
fig2 = []
u0 = -[0.7;0.7]
for i in 1:7
    tspan = (0.0,1.0)
    p = [10.0,-100.0]
    prob = ODEProblem(Subsystem1,u0,tspan,p)
    sol = solve(prob,Rodas4(),dtmax=1e-4, callback = CallbackSet(cb1,cb2))

    fig2 = plot!(sol,vars=(1,2),size = (400, 300), legend = false,
                framestyle = :box,
                linewidth = 1.5,
                xlabel = "\$x_1\$",
                labelfontsize = 15,
                linestyle = :dot,
                #xlims = (-2.5,2.5),
                xlims = (-1,1),
                xtickfont = font(8, "Akkurat"),
                ylabel = "\$x_2\$",
                #ylims = (-2.5,2.5),
                ylims = (-1,1),
                #linecolor = RGB(146/255,208/255,80/255),
                #linecolor = RGB(244/255,177/255,131/255),
                linecolor = RGB(0/255,0/255,255/255),
                ytickfont = font(8, "Akkurat"),
                grid = true,
                gridalpha = 0.5,
                gridstyle = :dash
                )
                u0[1] = u0[1] +0.1;
end
df = DataFrame(sol')
df.t = sol.t


savefig(fig2,"C:\\Users\\liemann\\Desktop\\Präsi_DFG\\example_hybrid_trajectory2.svg")


styles = filter((s->begin
                s in Plots.supported_styles()
            end), [:solid, :dash, :dot, :dashdot, :dashdotdot])
styles = reshape(styles, 1, length(styles))
n = length(styles)
y = cumsum(randn(20, n), dims = 1)
plot(y, line = (5, styles), label = map(string, styles), legendtitle = "linestyle")
