using Plots
using DifferentialEquations
using ModelingToolkit
using LinearAlgebra
using IfElse

function example1(du,u,p,t)
    #display(t)
    du[1] = u[1] + u[4]*u[2]
    du[2] = u[5]*u[1] + u[2]
    if u[3] < 0.0
        du[3] = p[1]*u[1] -u[2] -u[6]*u[3]
    else
        du[3] = u[2] -0.36*u[1] -u[6]*u[3]
    end
    du[4] = 0.0
    du[5] = 0.0
    du[6] = 0.0
end

function Hybrid()
    function example11(du, u, p, t)
        du[1] = u[1] + p[2] * u[2]
        du[2] = p[3] * u[1] + u[2]
        du[3] = IfElse.ifelse(
            p[4] > 0.0,
            p[1] * u[1] - u[2] - p[4] * u[3],
            u[2] - 0.36 * u[1] - p[4] * u[3],
        )
    end

    #condition_up(u,t,integrator) = integrator.p[1]*u[1]-u[2]
    #condition_down(u,t,integrator) = u[2]-0.36*u[1]
    x0 = [0.0,1.0,-1.0]
    p  = [2.75,-100.0,10.,1.0]
    speicher = hcat(0.0, p')

    condition(u,t,integrator) = u[3]
    function affect(integrator)
        tmp = deepcopy(integrator.p[2])
        integrator.p[2] = integrator.p[3]
        integrator.p[3] = tmp
        integrator.p[4] *= -1
        if integrator.p[4] > 0.0
            integrator.u[3] =
                -(0.36 * integrator.u[1] - integrator.u[2]) / integrator.p[4]
        else
            integrator.u[3] = -(integrator.u[2] - integrator.u[1] * integrator.p[1]) / integrator.p[4]
        end
        speicher = vcat(speicher,[integrator.t integrator.p[1] integrator.p[2] integrator.p[3] integrator.p[4]])
    end

    #cb1 = ContinuousCallback(condition_down,affect,nothing)
    #cb2 = ContinuousCallback(condition_up,affect,nothing)
    cb = ContinuousCallback(condition,affect)

    fun = ODEFunction(example11,mass_matrix=Diagonal([1,1,0]))
    prob =  ODEProblem(fun,x0,(0.0,0.2),p)
    sol = solve(prob,Rodas4(),callback=cb,dtmax=1e-3)
    return sol, prob,speicher #CallbackSet(cb1,cb2)
end
sol,prob,speicher = Hybrid()
plot(sol,vars=[(0,1),(0,2),(0,3)],legend=(0.7,0.5))
mtkprob = modelingtoolkitize(prob)
include("C:/Users/liemann/github/PowerDynamicsDevelopment/src/sensitivity_analyses/Local_Sensitivity.jl")


#Aufteilung in den ersten Ber
ind = indexin(speicher[:,1],sol.t)
sol_part = sol[ind[1]:ind[2]]
time_interval = (sol.t[ind[1]],sol.t[ind[2]])
sensis = TimeDomainSensitivies(mtkprob,time_interval,[0.0,1.0,-1.0],[2.75,-100.0,10.,1.0],[1,2],[1],sol_part)

eqs = equations(mtkprob)
st  = states(mtkprob)
p   = parameters(mtkprob)
@variables t
mtsys = ODESystem(eqs[1:2],t,st[1:2],p[2:3])
s = eqs[3]

condition(u,t,integrator) = p[1]*u[1]-u[2]

fun_ = ODEFunction(condition)
modelingtoolkitize(condition)

sx = GetJacobian([s],st[1:2])
fx,fy,gx,gy = GetSymbolicFactorizedJacobian(mtsys)
eqs, aeqs, D_states, A_states =GetSymbolicEquationsAndStates(mtsys)
params = parameters(mtsys)
sub_args = Pair.([D_states[1:2];params],[sol_part.u[end][1:2];[-100.0,10.]])
f_minus = Substitute(Num.(my_rhs.(equations(mtsys))),sub_args)
sub_args = Pair.([D_states[1:2];params],[sol.u[ind[2]+1][1:2];[10.0,-100.]])
f_plus = Substitute(Num.(my_rhs.(equations(mtsys))),sub_args)
sub_args = Pair.(p,speicher[2,2:5])
sx_minus = Substitute(Num.(sx),sub_args)
τx0 = sx_minus*sensis[3][1:2,end]./(sx_minus*f_minus)

sensis_neu = sensis[3][1:2,end] -(f_plus-f_minus).*τx0

plot(sol.t[1:end-1],sensis[3]')

###########new Test
include("C:/Users/liemann/github/PowerDynamicsDevelopment/src/grids/OLTC_Hybrid_Sensis.jl")
pg = OLTC_Hybrid_Sensi()

U,δ,ic = PowerFlowClassic(pg)
