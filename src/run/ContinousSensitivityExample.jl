using DifferentialEquations
using Plots
using SciMLSensitivity
using ForwardDiff


function f2(dx,x,p,t)
    Ud,w,L,C,R = p
    dx[1] = (Ud*cos(w*t) - x[2])/L
    dx[2] = (x[1] - x[2]/R)/C
    dx[3] = x[3] - (x[1] - C*dx[2]) 
end

M = [1. 0  0
     0  1. 0
     0  0  0];
u0 = [0.0,0.0,0.0]
p = [sqrt(2),100*pi,5e-3,1e-3,50]
odef = ODEFunction(f2,mass_matrix = M)
prob = ODEProblem(odef, u0, (0.0, 0.5), p)


sol = solve(prob, Rodas4(),dtmax=1e-3)[1,:]
plot(sol,idxs=(3))

###############################
sens_prob = ODEForwardSensitivityProblem(odef, u0, (0.0, 0.5), p,ForwardDiffSensitivity();)
sol = solve(sens_prob, Rodas4(),dtmax=1e-3);
x, dp = extract_local_sensitivities(sol);
da = dp[3];
plot(da[3,:], lw = 3)

#############################

function f3(theta)
    _prob = remake(prob,p=theta) #theta[6:end])
    sol_ = Array(solve(_prob,Rodas4(autodiff = false),dtmax=1e-3))'
end
res = ForwardDiff.jacobian(f3,p)
plot(res[:,1])
