using DifferentialEquations
using Plots
using SciMLSensitivity
using ForwardDiff
include("C:/Users/liemann/github/PowerDynamicsDevelopment/src/sensitivity_analyses/Local_Sensitivity.jl")

function Example1(;p=[2.75])
    function f(dx,x,p,t)
        dx[1] = x[1] + x[3]*x[2] 
        dx[2] = x[4]*x[1] + x[2]
        dx[3] = 0.0
        dx[4] = 0.0
        dx[5] = 0.0
        dx[6] = x[6]
    end
    M = [1. 0  0  0  0  0
         0  1. 0  0  0  0
         0  0  1. 0  0  0
         0  0  0  1. 0  0
         0  0  0  0  1. 0
         0  0  0  0  0  0];
    evr = Array{Float64}(undef,0,4)
    s1(u, t, integrator) = (integrator.p[1]*u[1] - u[2])/u[5] + u[5] -1
    s2(u, t, integrator) = (0.36*u[1] - u[2])/u[5] + u[5] + 1

    function h1(integrator)
        tmp =  deepcopy(integrator.u[3])
        integrator.u[3] = integrator.u[4] 
        integrator.u[4] = tmp
        integrator.u[5] = -integrator.u[5]
        evr = vcat(evr, [integrator.t 1 1 1])
    end
    function h2(integrator) 
        tmp =  deepcopy(integrator.u[3])
        integrator.u[3] = integrator.u[4] 
        integrator.u[4] = tmp
        integrator.u[5] = -integrator.u[5]
        evr = vcat(evr, [integrator.t 1 2 1])
    end
    cb1 = ContinuousCallback(s1,h1)
    cb2 = ContinuousCallback(s2,h2)

    u0 = [0.0,1,-100,10,1,0]
    
    ode_fun = ODEFunction(f,mass_matrix = M)#
    prob = ODEProblem(ode_fun, u0, (0.0, 0.2), p)
    sol = solve(prob, Rodas4(), callback = CallbackSet(cb1,cb2),dtmax=1e-4)
    return sol, evr
end

sol,evr = Example1();
plot(sol.t,sol[1,:])
plot!(sol.t,sol[2,:])
# plot(sol.t,sol[5,:])
# plot(sol.t,sol[6,:])

# plot(sol,idxs=(1,2),xlim=(-4,2),ylim=(-2,2))
# x = collect(-4:0.1:2)
# y = 0.36.*x
# plot!(x,y)
# plot!(x,2.75/0.36*y)

#######
mtk = modelingtoolkitize(sol.prob)
sym_states = states(mtk)
symp = parameters(mtk)

hs = Vector{Vector{Num}}(undef,0)
ident = Num.(vcat(sym_states[1:5],symp))
ident[3] = sym_states[4]
ident[4] = sym_states[3]
ident[5] = -sym_states[5]
push!(hs,ident)

@parameters t
s = Vector{Num}(undef,0)
push!(s,(symp[1]*sym_states[1]-sym_states[2])/sym_states[5]+sym_states[5]-1)
push!(s,(0.36*sym_states[1]-sym_states[2])/sym_states[5]+sym_states[5]+1)

hybrid_sen,Δτ = CalcHybridTrajectorySensitivity([mtk],sol,evr,s,hs);
# plot(sol.t,hybrid_sen[6][1,:])
# plot!(sol.t,hybrid_sen[6][2,:])

############
Δp = 3.0 - 2.75
sol_appr = ApproximatedTrajectory(sol,hybrid_sen[6],Δp)
sol_per,evr = Example1(p=[3.0]);

scatter(sol.t,sol[1,:])
scatter!(sol.t,sol_appr[1,:])
scatter!(sol_per.t,sol_per[1,:])
scatter!(sol.t,sol_refin[1,:])
#,xlims=(0.11,0.12),ylims=(-0.55,-0.49)
scatter(sol.t,hybrid_sen[6][1,:],xlims=(0.11,0.12))

plot(sol.t,sol[1,:],xlims=(0.18,0.2),ylims=(-0.05,-0.))
plot!(sol.t,sol_appr[1,:],xlims=(0.09,0.2),ylims=(-0.55,-0.))
plot(sol_per.t,sol_per[1,:],xlims=(0.18,0.2),ylims=(-0.05,-0.))
plot!(sol.t,sol_refin[1,:],xlims=(0.09,0.2),ylims=(-0.55,0))
plot!(sol.t,sol_refin[1,:],xlims=(0.18,0.2),ylims=(-0.05,-0.))


#Input: [mtk],sol,evr, sensis,Δτ x0_ind, Δx0
x0_ind = 6
Δx0 = 0.25
mtk = [mtk]
sensis = deepcopy(hybrid_sen)
# ind1 = Anfang der Indizes für die kontinuierlichen Sensis
ind1 = vcat(1,setdiff(indexin(evr[:,1],sol.t),[nothing]).+1)
# ind2 = Ende der Indizes für die kontinuierlichen Sensis
ind2 = vcat(setdiff(indexin(evr[:,1],sol.t),[nothing]),length(sol.t))
# zwischen den Anfang und End Indizes liegt das Event mit Zustandsänderung
ind_sol = Int.(hcat(ind1,ind2,zeros(size(ind1)[1],2)))

for i=1:size(ind_sol)[1]-1
    if Δτ[i,x0_ind] >= 0.0
        ind = findall(>(0),sol.t[ind_sol[i,2]+1] .<= sol.t .< sol.t[ind_sol[i,2]+1] + Δτ[i,x0_ind]*Δx0)
    else
        ind = findall(>(0),sol.t[ind_sol[i,2]+1] + Δτ[i,x0_ind]*Δx0 .< sol.t .< sol.t[ind_sol[i,2]+1])
    end
    ind_sol[i,3] = ind[1]
    ind_sol[i,4] = ind[end]
end


f_all = Vector{Array{Num,1}}(undef,length(mtk))
Fx_all = Vector{Array{Num,2}}(undef,length(mtk))
Fy_all = similar(Fx_all)
Gx_all = similar(Fx_all)
Gy_all = similar(Fx_all)
for (ind,val) in enumerate(mtk)
    fulleqs = equations(val)
    symstates = states(val)
    sym_params = parameters(val)
    eqs, aeqs, x, y = GetSymbolicEquationsAndStates(fulleqs, symstates)

    Fx_all[ind],Fy_all[ind],Gx_all[ind],Gy_all[ind] = GetSymbolicFactorizedJacobian(eqs, aeqs, x, y)
    f_all[ind] = my_rhs.(eqs)
end
sol_refin = zeros(size(sol))

#get indices
sym_states = states(mtk[1]);
sym_params = parameters(mtk[1]);

D0_states, A_states = GetFactorisedSymbolicStates(mtk[1]);

D0_indices= Int64.(setdiff(indexin(D0_states, sym_states), [nothing]))
A_indices = Int64.(setdiff(indexin(A_states, sym_states), [nothing]))

pk = sym_params .=> sol.prob.p
for i=1:size(ind_sol)[1]
    if i !=size(ind_sol)[1]
        if Δτ[i,x0_ind] >= 0.0
            
            f = f_all[Int(evr[i,2])]
            Gx = Gx_all[Int(evr[i,2])]
            Gy = Gy_all[Int(evr[i,2])]
            
            #Approximation until jump
            if i == 1
                tmp_sol = sol[:,ind_sol[i,1]:ind_sol[i,2]]
                tmp_sens = sensis[x0_ind][:,ind_sol[i,1]:ind_sol[i,2]]
                #display(length(tmp_sol))
                sol_refin[:,ind_sol[i,1]:ind_sol[i,2]] = ApproximatedTrajectory(tmp_sol,tmp_sens,Δx0)
            else # the solution and sensis after the new jumps have to be used
                tmp_sol = sol[:,ind_sol[i-1,4]+1:ind_sol[i,2]]
                tmp_sens = sensis[x0_ind][:,ind_sol[i-1,4]+1:ind_sol[i,2]]
                sol_refin[:,ind_sol[i-1,4]+1:ind_sol[i,2]] = ApproximatedTrajectory(tmp_sol,tmp_sens,Δx0)
            end

            xk = sym_states .=> sol[:,ind_sol[i,2]]  # the last values BEFORE the original jump -> are used
            f_f  = Float64.(Substitute(f,[xk;pk]))   # for all approximations within t+Δτ
            Gx_f = Float64.(Substitute(Gx,[xk;pk]))
            Gy_f = Float64.(Substitute(Gy,[xk;pk]))

            gygx = inv(Gy_f)*Gx_f*f_f

            x_part = sol_refin[D0_indices,ind_sol[i,2]]
            y_part = sol_refin[A_indices, ind_sol[i,2]]
 
            for (ind,j) in enumerate(sol.t[ind_sol[i,3]:ind_sol[i,4]])
                dt = j .- sol.t[ind_sol[i,2]]
                sol_refin[D0_indices, ind_sol[i,3]+ind-1] = x_part + f_f  * dt
                sol_refin[A_indices,ind_sol[i,3]+ind-1] = y_part + gygx * dt
            end
        else
            error("not implemented yet")
        end
    else
        tmp_sol = sol[:,ind_sol[i-1,4]+1:ind_sol[i,2]]
        tmp_sens = sensis[x0_ind][:,ind_sol[i-1,4]+1:ind_sol[i,2]]
        sol_refin[:,ind_sol[i-1,4]+1:ind_sol[i,2]] = ApproximatedTrajectory(tmp_sol,tmp_sens,Δx0)
    end
end

x1 = sol.t[ind_sol[2,2]]+Δτ[2,6]*0.25

plot!([x1,x1],[-1,1])

hybrid_sen[6][1,ind_sol[2,2]+1]

setdiff(indexin(sol.t[ind_sol[2,2]+1:end] .< sol.t[ind_sol[2,2]+1] +Δτ[2,6]*0.25,[1]),[nothing])

findall(>(0),sol.t[ind_sol[2,2]+1:end] .< sol.t[ind_sol[2,2]+1] +Δτ[2,6]*0.25)



