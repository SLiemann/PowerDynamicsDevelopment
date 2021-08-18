using PowerDynamics
using DifferentialEquations
using Plots
using JLD
#using ModelingToolkit
#using CSV
#using DataFrames

begin
    include("C:/Users/liemann/github/PowerDynamicsDevelopment/src/include_costum_nodes_lines_utilities.jl")
    include("C:/Users/liemann/github/PowerDynamicsDevelopment/src/grids/GFC_Test_Grid.jl")
    include("C:/Users/liemann/github/PowerDynamicsDevelopment/src/sensitivity_analyses/Local_Sensitivity.jl")

    pg = GFC_Test_Grid()
    U,δ,ic0 = PowerFlowClassic(pg, iwamoto = true, max_tol = 1e-7)

    pg1 ,ic = InitializeInternalDynamics(pg,ic0)
    params = GFC_params()
    prob = ODEProblem(rhs(pg1),ic,(0.0,5.0),params)
    pgsol,evr = simGFC(prob)
end

SaveComparingTrajectoryPlots(pg1,ic,pgsol,sens,labels_p,0.1)



params = GFC_params()
params[8] += 0.1*params[8]

display(plot(pgsol,["bus3"],:v, title = "Iabs (ohne CSA)"))
plot(pgsol,collect(keys(pg.nodes))[2:3],:v, title = "Spannung (ohne CSA)")

mtk_normal = GetMTKSystem(pg,(0.0,1.0),params)
mtk_fault = GetMTKSystem(GFC_Test_Grid(y_new = yfault()),(0.0,1.0),params)
mtk = [mtk_normal,mtk_fault]
s = GetTriggCondsGFCTest(mtk_normal)
h = GetStateResFunGFCTest(mtk_normal)
@time sens = CalcHybridTrajectorySensitivity(mtk,pgsol.dqsol,params,evr,s,h,[],collect(1:16))    #collect(1:15)
save("C:/Users/liemann/Desktop/Sens_GFC_Test/sens_all_CSA.jld", "sens", sens,"ic0",ic0,"p_pre",params,"evr",evr,"sensis_p",collect(1:16))

plot(pgsol,["bus3"],:i_abs, label = "Kvi = " * string(pg.nodes["bus3"].Kvi) * ", K_vq = " *string(pg.nodes["bus3"].K_vq))
plot(pgsol,["bus3"],:i_abs, label = "Strom VSC")
plot!(pgsol,["bus3"],:i_setabs1, label = "Regelungsstrom VSC")
ylims!((0.2,1.6))
xlims!((0.99,2.5))
plot!(pgsol,collect(keys(pg.nodes))[2],:v, label = "U_3", legend = (0.9,0.8))
plot!(xlabel = "\$t\$")
plot!(                grid = true,
                gridalpha = 0.25,
                gridstyle = :dash)

plot(pgsol,["bus3"],:i_abs, label = "mit CSA, imax = "* string(pg.nodes["bus3"].imax_csa))
plot(pgsol,["bus3"],:p, label = "imax_csa = " * string(pg.nodes["bus3"].imax_csa), title = "P_VSC (AS)")
plot(pgsol,["bus3"],:q)
plot!(pgsol,["bus3"],:θ)
plot(pgsol,["bus3"],:φ)
plot(pgsol,["bus3"],:ω)
plot(pgsol,["bus3"],:Qm)
xlims!((4.9,5.50))
plot(pgsol,["bus3"],:e_uq)
plot(pgsol,["bus3"],:e_ud)
plot(pgsol,["bus3"],:e_id)
plot!(pgsol,["bus3"],:e_iq)

sol = ExtractResult(pgsol,:ω_3)

f, freq = DFT(sol,pgsol.dqsol.t)
plot(freq,abs.(f), xlim=(-10, 10))
#Calculating approximated trajectory
#sol_original = deepcopy(pgsol.dqsol)
#pgsol_or = deepcopy(pgsol)
sol_appr = deepcopy(sol_original)
for (ind,val) in enumerate(collect(eachcol(sens[15])))
    sol_appr.u[ind+1] .+= val*(0.01)
end
pgsol_tmp = PowerGridSolution(sol_appr,pg)
plot!(pgsol_tmp,"bus3",:i_abs, label = "approximated")
plot!(pgsol,"bus3",:i_abs, label = "real perturbed", xlims = (0.0,1.0))
plot(pgsol_or,"bus3",:i_abs, label = "Original")
plot!(pgsol_tmp,collect(keys(pg.nodes)),:v,legend = fals


#Plot all parameter sensitivities
labels_p = [
    "Kp_droop", #1
    "Kq_droop", #2
    "ωf_P", #3
    "ωf_Q", #4
    "xlf", #5
    "rf", #6
    "xcf", #7
    "Kp_u", #8
    "Ki_u", #9
    "Kp_i", #10
    "Ki_i", #11
    "imax", #12
    "Kvi", #13
    "σXR", #14
    "K_vq", #15
    "imax_csa" #16
    ]
syms = rhs(pg).syms
look_on = 14
plot(pgsol.dqsol.t[1:end-1],sens[1][look_on,1:end], title = "Sensis of $(String(syms[look_on]))",
    label = labels_p[1],
    legend = :outertopright,
    size = (1000,750))
for i in collect(16:16)
    display(plot!(pgsol.dqsol.t[1:end-1],sens[i][look_on,1:end], label = labels_p[i]))
    #sleep(3.0)
end

#Calculating "per unit" sensis
sens_pu = deepcopy(sens)
for i in collect(1:length(sens))
    sens_pu[i] = sens_pu[i]*params[i]
end
look_on = 14
plot(pgsol.dqsol.t[1:end-1],sens_pu[2][look_on,1:end], title = "Sensis of $(String(syms[look_on]))",
    label = labels_p[2],
    legend = :outertopright,
    size = (1000,750))
for i in collect(2:15)
    display(plot!(pgsol.dqsol.t[1:end-1],sens_pu[i][look_on,1:end], label = labels_p[i]))
    #sleep(3.0)
end

#Calculating
Imax = 1.0
I_tresh = 0.95
ΔXR = 10.0
V0 = 1.0
a = (Imax - I_tresh)^2 * (1 + ΔXR^2)
b = 2.0*(Imax - I_tresh) * (ΔXR)
c = - V0^2/Imax^2

kp = (-b + sqrt(b^2 - 4.0*a*c))/(2*a)


plot(ExtractResult(pgsol,"bus2",:u_r))
plot!(ExtractResult(pgsol,"bus2",:u_i))

i = ExtractResult(pgsol,"bus3",:i_abs)
t = pgsol.dqsol.t
using CSV
using DataFrames

time_Series = DataFrame(zeit = t, strom = i)
CSV.write("C:\\Users\\liemann\\Desktop\\test.csv",time_Series, delim = ";")

function SavePlotApprTrajectories(pg,sol_or,sol_per,sens,par_num,or_value,new_value,labels_p; state_sym = :i_abs)
    sol_appr = deepcopy(sol_or.dqsol)
    for (ind,val) in enumerate(collect(eachcol(sens[par_num])))
        sol_appr.u[ind+1] .+= val*(new_value-or_value)
    end
    pgsol_tmp = PowerGridSolution(sol_appr,pg)
    title_str  = labels_p[par_num] * ": from " * string(or_value) * " to " * string(new_value)
    plot(sol_or,"bus3",state_sym, label = "Original - " * string(state_sym), title = title_str)
    plot!(sol_per,"bus3",state_sym, label = "Real perturbed - " * string(state_sym))
    display(plot!(pgsol_tmp,"bus3",state_sym, label = "Approximated - " * string(state_sym),linestyle = :dash))
    savefig("C:\\Users\\liemann\\Desktop\\Julia_figures\\"*labels_p[par_num] *".svg")
end


function SaveComparingTrajectoryPlots(pg,ic_or,pgsol_or,sens,labels_p,delta)
    params = GFC_params()
    for (ind,val) in enumerate(params)
        params_new = deepcopy(params)
        params_new[ind] += delta*params_new[ind]
        prob = ODEProblem(rhs(pg),ic,(0.0,3.0),params_new)
        try
            pgsol,evr = simGFC(prob)
            SavePlotApprTrajectories(pg,pgsol_or,pgsol,sens,ind,params[ind],params_new[ind],labels_p)
        catch
            @warn "Could not calc: " *labels_p[ind]
        end
    end
end
