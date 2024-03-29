using PowerDynamics
using DifferentialEquations
using Plots
using MAT
#using DiffEqSensitivity
#using JLD
#using ModelingToolkit
#using CSV
#using DataFrames

begin
    include("C:/Users/liemann/github/PowerDynamicsDevelopment/src/include_custom_nodes_lines_utilities.jl")
    include("C:/Users/liemann/github/PowerDynamicsDevelopment/src/operationpoint/InitializeInternals.jl")
    include("C:/Users/liemann/github/PowerDynamicsDevelopment/src/grids/GFC_Test_Grid.jl")
    #include("C:/Users/liemann/github/PowerDynamicsDevelopment/src/sensitivity_analyses/Local_Sensitivity.jl")
end
begin
    pg = GFC_Test_Grid()
    U,δ,ic0 = PowerFlowClassic(pg, iwamoto = true, max_tol = 1e-7)
    pg1 ,ic = InitializeInternalDynamics(pg,ic0)
    params = GFC_params()
    prob = ODEProblem(rhs(pg1),ic,(0.0,2.0),params, initializealg = BrownFullBasicInit())
    #prob_new = ODEForwardSensitivityProblem(rhs(pg1),ic,(0.0,0.1),params)
    pgsol,evr = simGFC(prob)
end
plot(pgsol,collect(keys(pg.nodes))[2:end],:v)
plot(pgsol,["bus3"],:i_abs, label = "Iabs", legend = false)#, ylims=(1.1,1.4)









plot(pgsol,["bus3"],:θ, label ="Droop-Winkel VSC")
theta = ExtractResult(pgsol,:θ_3)*180/pi
plot(pgsol,["bus3"],:ω)
plot(pgsol,["bus3"],:Pout)
plot(pgsol,["bus3"],:UQmeas)

file = matopen("C:\\Users\\liemann\\Desktop\\ucf.mat")
tmp2 = read(file, "ucf")'
close(file)
ucf = tmp2
plot(ucf[:,1],ucf[:,2],label = "MATLAB - Ucf")
plot!(pgsol,collect(keys(pg.nodes))[3],:v, label = "PD - Ucf", legend = false,xlims=(0.0,1.5))
plot(ucf[:,1],ucf[:,3],label = "MATLAB - Iabs")
plot!(pgsol,["bus3"],:i_abs, label = "PD - Iabs", legend = (0.8,0.5),xlims=(0.0,1.5)) #, xlims =(5.0,5.5)
plot(ucf[:,1],ucf[:,4].+ic[7],label = "MATLAB - θ")
plot!(pgsol,["bus3"],:θ, label = "PD - θ", legend = (0.8,0.5),xlims=(0.6,1.0),ylims=(6.1,7))

plot(ucf[:,1],ucf[:,5],label = "MATLAB - id")
plot!(pgsol,["bus3"],:i_abs, label = "PD - id", legend = (0.8,0.5),xlims=(-5,0.0))

plot(ucf[:,1],ucf[:,6],label = "MATLAB - iq")
plot!(pgsol,["bus3"],:i_absq, label = "PD - iq", legend = (0.8,0.5))


file2 = matopen("C:\\Users\\liemann\\Desktop\\ibf.mat")
tmp3 = read(file2, "ibf")'
close(file2)
ibf = tmp3
plot(ibf[:,1],ibf[:,2],label = "MATLAB - id (global)")
plot!(pgsol,["bus3"],:i_absd, label = "PD - id", legend = (0.8,0.5))

plot(ibf[:,1],ibf[:,3],label = "MATLAB - iq(global)")
plot!(pgsol,["bus3"],:i_absd, label = "PD - iq", legend = (0.3,0.2),xlims=(0.6,1.2))

id_tmp = ExtractResult(pgsol,:i_absd_3)
plot(id_tmp)
plot(iq_tmp)
iq_tmp = ExtractResult(pgsol,:i_absd_3)
plot(sqrt.(id_tmp[1:end-1].^2+iq_tmp.^2))
plot(pgsol.dqsol.t,abs.(id_tmp[1:end-1].^2+iq_tmp.^2))
sol_try, evr = simGFC(prob_new)
x,dp = extract_local_sensitivities(sol_try)

mtk_normal = GetMTKSystem(pg1,(0.0,1.0),params)
mtk_fault = GetMTKSystem(GFC_Test_Grid(y_new = yfault()),(0.0,1.0),params)
mtk = [mtk_normal,mtk_fault]
s = GetTriggCondsGFCTest(mtk_normal)
h = GetStateResFunGFCTest(mtk_normal)
@time sens = CalcHybridTrajectorySensitivity(
    mtk,
    pgsol.dqsol,
    params,
    evr,
    s,
    h,
    [],
    collect(1:16),
)
save("C:/Users/liemann/Desktop/Sens_GFC_Test/sens_REAL_Param.jld", "sens", sens,"ic0",ic0,"p_pre",params,"evr",evr,"sensis_p",collect(1:15))
#tmp = JLD.load("C:/Users/liemann/Desktop/Sens_GFC_Test/sens_all_CSA2.jld")

plot(pgsol,["bus3"],:i_abs, label = "Iabs", ylims=(0.0,1.4))
plot(pgsol,collect(keys(pg.nodes))[3],:v, label = "Spannung VSC", legend = (0.8,0.75))
plot(pgsol,["bus3"],:θ, label ="Droop-Winkel VSC",xlims=(0.0,0.5))
plot(pgsol,["bus3"],:p)
plot(pgsol,["bus3"],:Pout)
plot(pgsol,collect(keys(pg.nodes))[2:3],:φ)

theta = plot(pgsol,["bus3"],:θ)[1][1][:y]
p_vsc = plot(pgsol,["bus3"],:Pout)[1][1][:y]
plot(theta,p_vsc, legend =false)

#theta = ExtractResult(pgsol,:θ_3)*180/pi
#p_vsc = ExtractResult(pgsol,:Pout_3)
theta = plot(pgsol,["bus3"],:θ)[1][1][:y]*180/pi
p_vsc = plot(pgsol,["bus3"],:p)[1][1][:y]
plot(theta,p_vsc, label = "P(δ) - zeitlicher Verlauf", xlims=(0.0,180.), ylims = (0.0,10.0))
winkel = collect(0:0.01:pi)
P_δ = 1.0/(0.1+0.15/12)*sin.(winkel) #+0.0177*(380/320)^2+0.00059095*(380/320)^2+
plot!([0.0,180.0],[3.0,3.0], linestyle = :dash, color = "black", label = "Psoll")
plot!(winkel*180/pi,P_δ, color = "grey", label = "P(δ) - statisch")
xlabel!("δ_VSC in °")
ylabel!("P in p.u.")

plot(pgsol,["bus3"],:q)
plot!(pgsol,["bus3"],:θ)
plot(pgsol,["bus3"],:φ)
plot(pgsol,["bus3"],:ω)
plot(pgsol,["bus3"],:Qm)
plot(pgsol,["bus3"],:e_uq)
plot(pgsol,["bus3"],:e_ud)
plot(pgsol,["bus3"],:e_id)
plot!(pgsol,["bus3"],:e_iq)

plot(pgsol,collect(keys(pg.nodes))[3],:v, label = "Spannung VSC", legend = (0.6,0.85))
plot!([0.0,0.2],[0.9,0.9], linestyle = :dash, color = "red", label= "FRT- Limiting curve")
plot!([0.2,0.2],[0.9,0.05], linestyle = :dash, color = "red", label= nothing)
plot!([0.2,0.45],[0.05,0.05], linestyle = :dash, color = "red", label= nothing)
plot!([0.45,0.65],[0.05,0.60], linestyle = :dash, color = "red", label= nothing)
plot!([0.65,3.0],[0.60,0.9], linestyle = :dash, color = "red", label= nothing)
plot!([3.0,5.0],[0.90,0.9], linestyle = :dash, color = "red", label= nothing)


uqmeas3 = ExtractResult(pgsol,:u_r_3)
plot(pgsol.dqsol.t,uqmeas3)

#Calculating approximated trajectory
sol_original = deepcopy(pgsol.dqsol)
pgsol_or = deepcopy(pgsol)
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
    #"imax_csa" #16
    ]
syms = rhs(pg).syms
look_on = 15
plot(pgsol.dqsol.t[1:end-1],sens[1][look_on,1:end], title = "Sensis of $(String(syms[look_on]))",
    label = labels_p[1],
    legend = :outertopright,
    size = (1000,750))
for i in setdiff(collect(2:15),[2,5,6,7])
    display(plot!(pgsol.dqsol.t[1:end-1],sens[i][look_on,1:end], label = labels_p[i]))
    sleep(3.0)
end

#Calculating "per unit" sensis
sens_pu = deepcopy(sens)
for i in collect(1:length(sens))
    sens_pu[i] = sens_pu[i]*params[i]
end
look_on = 15
plot(pgsol.dqsol.t[1:end-1],sens_pu[1][look_on,1:end], title = "Sensis of $(String(syms[look_on])) in p.u.",
    label = labels_p[1],
    legend = :outertopright,
    size = (1000,750))
for i in setdiff(collect(2:15),[2,5,6,7])
    display(plot!(pgsol.dqsol.t[1:end-1],sens_pu[i][look_on,1:end], label = labels_p[i]))
    sleep(3.0)
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

function SavePlotApprTrajectories(pg_tmp,sol_or,sol_per,sens,par_num,or_value,new_value,labels_p; state_sym = :i_abs)
    sol_appr = deepcopy(sol_or.dqsol)
    for (ind,val) in enumerate(collect(eachcol(sens[par_num])))
        sol_appr.u[ind+1] .+= val*(new_value-or_value)
    end
    pgsol_tmp = PowerGridSolution(sol_appr,pg_tmp)
    title_str  = labels_p[par_num] * ": from " * string(or_value) * " to " * string(new_value)
    plot(sol_or,"bus3",state_sym, label = "Original - " * string(state_sym), title = title_str)
    plot!(sol_per,"bus3",state_sym, label = "Real perturbed - " * string(state_sym))
    display(plot!(pgsol_tmp,"bus3",state_sym, label = "Approximated - " * string(state_sym),linestyle = :dash))
    savefig("C:\\Users\\liemann\\Desktop\\Julia_figures\\"*labels_p[par_num] *".svg")
end


function SaveComparingTrajectoryPlots(pg_tmp,ic_or,pgsol_or,sens,labels_p,delta)
    params_tmp = GFC_params()
    for (ind,val) in enumerate(params_tmp)
        params_new = deepcopy(params_tmp)
        params_new[ind] += delta*params_new[ind]
        prob_tmp = ODEProblem(rhs(pg_tmp),ic_or,(0.0,3.0),params_new)
        try
            pgsol_per,evr = simGFC(prob_tmp)
            SavePlotApprTrajectories(pg_tmp,pgsol_or,pgsol_per,sens,ind,params[ind],params_new[ind],labels_p)
        catch
            @warn "Could not calc: " *labels_p[ind]
        end
    end
end
