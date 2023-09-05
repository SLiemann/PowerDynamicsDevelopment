using DifferentialEquations
using Plots
using LinearAlgebra
using Roots

function f(du,u,p,t)
    Ud,w,Pdc,Cd = p
    du[1] = u[1] - Ud*sin(w*t)#*exp(-t)
    du[2] = 0.0 # u_off
    du[3] = 0.0 # t_sum
    du[4] = 0.0 # t_on
    du[5] = 0.0 # t_off
    du[6] = 0.0 # half_cycle
    du[7] = 0.0 # a
    du[8] = 0.0 # b
    du[9] = u[9] - Ud#*exp(-t)
end

u0 = [0.0, 0.7*230*sqrt(2),0.0,0.0025,0.00511210289264389,0.0,0.0,0.0,230*sqrt(2)]
p =  [230*sqrt(2),100*pi,1000.0,700e-6]

M = Diagonal([0, 1, 1, 1, 1, 1, 1, 1,0])

s0(u, t, integrator) = true
function h0(integrator) 
    integrator.u[3] = integrator.uprev[3] + integrator.t -integrator.tprev   
end

s1(u, t, integrator) = u[3] >= 0.01
function h1(integrator) 
    Ud,w,Pdc,Cd = integrator.p
    integrator.u[6] = 1.0
    integrator.u[3] = mod(integrator.t,0.01)
    if integrator.u[4] > 0
        u_off = integrator.u[9]*sin(w*integrator.u[5])
        integrator.u[2] = sqrt(u_off^2 - 2*Pdc*(0.01-integrator.u[5])/Cd)
    else
        u_off = integrator.u[9]*sin(w*integrator.u[5])
        integrator.u[2] = sqrt(integrator.u[2]^2 - 2*Pdc*(0.01)/Cd)
    end
end    

s2(u, t, integrator) = u[3] < 0.01 && u[6] == 1.0
function h2(integrator) 
    integrator.u[6] = 0.0
end 

s3(u, t, integrator) = u[6] == 1.0 || mod(t,0.01) < u[5]
function h3(integrator) 
    Ud,w,Pdc,Cd = integrator.p
     integrator.u[5] = (pi + asin(2*Pdc/(w*Cd*integrator.u[9]^2)))/(2*w)
end 

s4(u, t, integrator) = u[6] == 1.0 || mod(t,0.01) < u[4]
function h4(integrator)
    Ud,w,Pdc,Cd = integrator.p
    uofft2 = integrator.u[2]
    u0 = integrator.u[9]
    x1 = 0.004517042542168
    x2 = −0.084973441092720
    x3 =  1.367157627046003
    x4 = −0.580239811988436
    function fton(ton)
        ton^3*(u0^2*w^3*x4) + ton^2*(u0^2*w^2*x3) + ton*(u0^2*w*x2 + 2*Pdc/Cd) + u0^2*x1 - uofft2^2  
    end
     integrator.u[4] = find_zero(fton,0.004)
end 

cb0 = DiscreteCallback(s0,h0)
cb1 = DiscreteCallback(s1,h1)
cb2 = DiscreteCallback(s2,h2)
cb3 = DiscreteCallback(s3,h3)
cb4 = DiscreteCallback(s4,h4)

ode_fun = ODEFunction(f,mass_matrix = M)
prob = ODEProblem(ode_fun, u0, (0.0, 1), p)
sol = solve(prob, Rodas4P(), callback = CallbackSet(cb0,cb1,cb2,cb3,cb4),dtmax=1e-4);

plot(sol,idxs=[4,5])

plot(sol.t.-sol[9,:])



f(t) = t^3 -5

find_zero(f,1)
