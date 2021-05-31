using PowerDynamics
using DifferentialEquations
using Plots
using ModelingToolkit
using CSV
using DataFrames

#data1 = sqrt.((sensi[2][7, :] .^ 2) .+ (sensi[2][8, :] .^ 2))
#plot(pgsol.dqsol.t[1:end-1], sensi[2][3, :])

begin
    include("C:/Users/liemann/github/PowerDynamicsDevelopment/src/grids/OLTC_Hybrid_Sensis.jl")
    include("C:/Users/liemann/github/PowerDynamicsDevelopment/src/sensitivity_analyses\\Local_Sensitivity.jl")
end

mtk = GetMTKOLTCSystem()
pg, ic = GetInitializedOLTCHisken()
tspan = tspanOLTCHisken()
params = GetParametersOLTCHisken()
pgsol,evr = SimulateOLTCHIsken()
sol = pgsol.dqsol
#sensi = TimeDomainSensitivies(mtk, tspan, ic, params, [1], [3], pgsol.dqsol)
begin
    plot(pgsol,["bus4"],:v,label="PowerDynamics")
    test = DataFrame(CSV.File("C:\\Users\\liemann\\Desktop\\uPF.csv"; header=false, delim=';', type=Float64))
    plot!(test.Column1,test.Column2,label = "PF-bus3")
    ylims!((0.82,1.02))
    xlims!((0.0,200.0))
    yticks!(collect(0.82:0.02:1.02))
    xticks!(collect(0.0:20.:200.))
end
toll = CalcHybridTrajectorySensitivity(mtk,sol,evr,[1],[3])
plot(sol.t[1:end-2],toll[2][7,1:end-1])
ylims!((-0.008,0.024)) #.+toll[2][8,1:end-1]
xlims!((0.0,200.0))
yticks!(collect(-0.008:0.002:0.024))
xticks!(collect(0.0:20.:200.))

function CalcHybridTrajectorySensitivity(mtk,sol,evr,u0_sensi,p_sensi)
    ic = sol.prob.u0
    p  = sol.prob.p
    xx0_k, yx0_k, sym_states,sym_params, A_states, D_states, M, N, O, symp, Δt,len_sens, f, g, J =
        InitTrajectorySensitivity(mtk, ic, p, u0_sensi, p_sensi)
    xx0 = [i[1] for i in xx0_k]
    yx0 = [i[1] for i in yx0_k]
    fx,fy,gx,gy = J
    hx,hy,sx,sy = CalcTriggerAndStateResetJacobians(mtk)
    sensis = Vector{Array{Float64}}(undef, len_sens)
    for i = 1:length(sensis)
      sensis[i] = Array{Float64}(
        undef,
        size(D_states)[1] + size(A_states)[1],
        size(sol)[2] - 1,
      )
    end

    ind_sol = vcat(1,setdiff(indexin(evr[:,1],sol.t),[nothing]),indexin(sol.t[end],sol.t))
    for i = 1:length(ind_sol)-1
        sol_part = sol[ind_sol[i]:ind_sol[i+1]]
        sensi_part,xx0_k,yx0_k = ContinuousSensitivity(
                                    sol_part,
                                    xx0_k,
                                    yx0_k,
                                    sym_states,
                                    A_states,
                                    D_states,
                                    M,
                                    N,
                                    O,
                                    symp,
                                    Δt,
                                    len_sens,
                                )
        for j = 1:length(sensi_part)
            sensis[j][:,ind_sol[i]:ind_sol[i+1]-1] = sensi_part[j]
        end
        if ind_sol[i+1] ≠ ind_sol[end]
            display("Event at t = $(evr[i,1])")

            xx0_pre = [i[2] for i in xx0_k]
            yx0_pre = [i[2] for i in yx0_k]
            x0_pre = sol.u[ind_sol[i+1]]    # ind_sol[i+1] is before jump
            x0_post = sol.u[ind_sol[i+1]+1] # ind_sol[i+1]+1 is after jump
            p_pre = p
            p_post = evr[i,2:end-2]
            hx_tmp = hx[Int(evr[i,end])]
            hy_tmp = hy[Int(evr[i,end])]
            sx_tmp = sx[Int(evr[i,end-1])]
            sy_tmp = sy[Int(evr[i,end-1])]

            #display("Befor Jump")
            #display(xx0_pre)
            #display(yx0_pre)


            xx0_post, yx0_post = CalcSensitivityAfterJump(
                                    sym_states,
                                    sym_params,
                                    xx0_pre,
                                    yx0_pre,
                                    x0_pre,
                                    x0_post,
                                    p_pre,
                                    p_post,
                                    f,
                                    g,
                                    J,
                                    hx_tmp,
                                    hy_tmp,
                                    sx_tmp,
                                    sy_tmp,
                                )
            #=
            display("After Jump")
            display([yx0_post./yx0_pre;xx0_post./xx0_pre])

            xx0_post, yx0_post = xx0_pre, yx0_pre
            =#
        else
            return sensis
        end

        xx0_k = xx0 .=> xx0_post
        yx0_k = yx0 .=> yx0_post
        symp = sym_params .=> p_post
    end
end
