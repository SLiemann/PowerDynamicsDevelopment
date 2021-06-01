using PowerDynamics
using DifferentialEquations
using Plots
using ModelingToolkit
using CSV
using DataFrames

begin
    include("C:/Users/liemann/github/PowerDynamicsDevelopment/src/include_costum_nodes_lines_utilities.jl")
    include("C:/Users/liemann/github/PowerDynamicsDevelopment/src/grids/OLTC_Hybrid_Sensis.jl")
    include("C:/Users/liemann/github/PowerDynamicsDevelopment/src/sensitivity_analyses\\Local_Sensitivity.jl")
end

mtk = GetMTKOLTCSystem()
pgsol,evr = SimulateOLTCHIsken()
sol = pgsol.dqsol
begin
    plot(pgsol,["bus4"],:v,label="PowerDynamics")
    #test = DataFrame(CSV.File("C:\\Users\\liemann\\Desktop\\uPF.csv"; header=false, delim=';', type=Float64))
    #plot!(test.Column1,test.Column2,label = "PF-bus3")
    ylims!((0.84,0.855))
    xlims!((30.0,30.05))
    #yticks!(collect(0.82:0.02:1.02))
    xticks!(collect(30.0:0.01:30.05))
end
toll = CalcHybridTrajectorySensitivity(mtk,sol,evr,[1],[3])
plot(sol.t[1:end-1],toll[2][7,1:end]) #
ylims!((-0.008,0.024)) #.+toll[2][8,1:end-1]
xlims!((30.0099,30.0101))
yticks!(collect(-0.008:0.002:0.024))
xticks!(collect(0.0:20.:200.))

function CalcHybridTrajectorySensitivity(mtk,sol,evr,u0_sensi,p_sensi)
    ic = sol.prob.u0
    p_pre  = GetParametersOLTCHisken()
    xx0_k, yx0_k, sym_states,sym_params, A_states, D_states, M, N, O, symp, Δt,len_sens, f, g, J =
        InitTrajectorySensitivity(mtk, ic, p_pre, u0_sensi, p_sensi)
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
    ind_sol = vcat(1,setdiff(indexin(evr[:,1],sol.t).+1,[nothing]),length(sol.t))
    #ind_sol = [1]
    #for i in evr[:,1] # DifferentialEquations.jl has multiple time points
    #    ind_sol = vcat(ind_sol,findall(x->x==i,sol.t)[end])
    #end
    #ind_sol = vcat(ind_sol,length(sol.t))

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
            x0_pre = sol.u[ind_sol[i+1]]    # ind_sol[i+1]-1 is before jump
            x0_post = sol.u[ind_sol[i+1]+1]    # ind_sol[i+1] is after jump
            p_post = evr[i,2:end-2]
            hx_tmp = hx[Int(evr[i,end])]
            hy_tmp = hy[Int(evr[i,end])]
            sx_tmp = sx[Int(evr[i,end-1])]
            sy_tmp = sy[Int(evr[i,end-1])]

            display("p_pre")
            display(p_pre)

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
        else
            return sensis
        end
        display("p_post")
        display(p_post)

        xx0_k = xx0 .=> xx0_post
        yx0_k = yx0 .=> yx0_post
        symp = sym_params .=> p_post
        p_pre = p_post
    end
end
