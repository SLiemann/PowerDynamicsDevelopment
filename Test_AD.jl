using DifferentialEquations
using Plots
using DiffEqSensitivity, OrdinaryDiffEq, ForwardDiff, Zygote

ForwardDiff.can_dual(::Type{ComplexF64{Int64}}) = true
ForwardDiff.can_dual(::Type{Any}) = true

function f(du,u,p,t)
  du[1] = dx = p[1]*u[1] - p[2]*u[1]*u[2]
  du[2] = dy = -p[3]*u[2] + u[1]*u[2]
end


f(x) = sum(sin, x) + prod(tan, x) * sum(sqrt, x);
f2(x) = sin(x) + cos(2*x)

xi = [0.1+1im*0.2,0.1+1im*0.2*8]
x0 = collect(0:1e-3:0.5)
xd = Complex.(ForwardDiff.Dual.(rand(5)))

g = x -> ForwardDiff.derivative(f2, x)

dx = g.(x0)

plot(x0,f2.(x0))
plot!(x0,dx)


f(x0)


function abc(y,x)
    y[1] = abs(3.0*x[1] + 1im*x[2])
    y[2] = x[2]
    nothing
end

y = [0.0,0.0]
dy = [0 0.0; 0 0]
abc(y,[3,4])
ForwardDiff.jacobian!(dy,abc,y,[3.0,1])




function fiip(du,u,p,t)
  du[1] = dx = p[1]*u[1] - p[2]*u[1]*u[2]
  du[2] = dy = -p[3]*u[2] + p[4]*u[1]*u[2]
end
p = [1.5,1.0,3.0,1.0]; u0 = [1.0;1.0]
prob = ODEProblem(fiip,u0,(0.0,10.0),p)
sol = solve(prob,Tsit5(),reltol=1e-6,abstol=1e-6)
plot(sol)
sum(sol)

du01,dp1 = Zygote.gradient((u0,p)->sum(solve(prob,Tsit5(),u0=u0,p=p,saveat=0.1,sensealg=QuadratureAdjoint())),u0,p)



function sum_of_solution(x)
    _prob = remake(prob,u0=x[1:2],p=x[3:end])
    sol = solve(_prob,Tsit5(),reltol=1e-6,abstol=1e-6,saveat=0.1)
    sol.u
end
dv = rand(101,6)

dx = ForwardDiff.jacobian(sum_of_solution,[u0;p])

tmp = values(dx[1])
tmp2 = values(tmp)
tmp3 = values(tmp3)

tmp_f = x -> ForwardDiff.value(x[2])

dx1 = tmp_f.(dx)
plot(dx1)


function RL(du,u,p,t)
    ω  = p[1]
    U1 = p[2]
    U2 = p[3]
    K = p[4]
    L = p[5]
    R = p[6]

    du[1] = u[1] - U1*sin(ω*t)
    du[2] = (u[1] - u[2]*R - U2*sin(ω*t+pi/6)) / L
    du[3] = abs2(u[1]+1im*u[2])
end
using LinearAlgebra
mm = Diagonal([false,true,true])
fun = ODEFunction(RL,mass_matrix =  mm)
u0 = [0.0,0.0,0.0]
p = [100*pi,230*sqrt(2),230*sqrt(2),0.1,1/(100*pi),1]

prob = ODEForwardSensitivityProblem(fun,u0,(0.0,0.1),p,ForwardDiffSensitivity())
sol = solve(prob,Rodas4(),dtmax=1e-5)
#ForwardDiff.value.(sol.u)
#plot(sol)

x,dp = extract_local_sensitivities(sol)
da = dp[6]

plot(sol.t,x[3,:],lw=3)


include("C:/Users/liemann/github/PowerDynamicsDevelopment/src/include_custom_nodes_lines_utilities.jl")
include("C:/Users/liemann/github/PowerDynamicsDevelopment/src/grids/OLTC_Hybrid_Sensis.jl")

pg, ic0 = GetInitializedOLTCHisken();
rhs_pg = rhs(pg);
p = GetParametersOLTCHisken(5.0);
prob = ODEForwardSensitivityProblem(rhs_pg,ic0,(0.0,0.1),p,ForwardDiffSensitivity())
sol = solve(prob,Rodas4(),dtmax=1e-3)

x,dp = extract_local_sensitivities(sol)



a,b = SimulateOLTCHIsken()
x,dp = extract_local_sensitivities(a)

u4 = sqrt.(x[8,:].^2 .+x[7,:].^2)
plot(a.t,u4)



plot(a.t,dp[2][:,:]')
plot(a,"bus4",:v)


begin
    include("C:/Users/liemann/github/PowerDynamicsDevelopment/src/grids/LTVS_Test_System_TapParam.jl")
    pg = GFC_LTVS_Test_SystemTapParam(nTap = 5)
    Qmax   = [Inf, Inf, Inf,Inf, Inf,Inf*sqrt(1-0.95^2)]
    Qmin   = -Qmax
    U,δ,ic0 = PowerFlowClassic(pg,iwamoto = true, Qmax = Qmax, Qmin = Qmin, Qlimit_iter_check = 2,max_tol = 1e-6)
    pg, ic0 = InitializeInternalDynamics(pg,ic0)
    sol,evr  = run_LTVS_simulationTapParam(pg,ic0,(0.0,70.0))
end

plot(pgsol_per,"bus4",:i_abs, legend=:bottomleft)
plot(pgsol_per,"bus4",:v)
#display(plot!(pgsol_stkvq,"bus4",:v,label = "Enh kvq"))
#display(plot(pgsol,"bus4",:v,label = "U-Original",xlims=(5,70),ylims=(0.9,1.01), legend = (0.5,0.1)))
#display(plot!(pgsol_per,"bus4",:v,label = "real perturbed"))
#display(plot!(pgsol_per,"bus4",:v, label ="real perturbed")) #linestyle = :dash
