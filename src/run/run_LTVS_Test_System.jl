using PowerDynamics
#using OrderedCollections: OrderedDict
using Plots
#import PowerDynamics: PiModel
using DifferentialEquations
#using CSV #read PF DataFrames
#using DataFrames #for CSV
#using Distributed
using JLD

Ubase = 380e3
Sbase = 100e6
Zbase = (Ubase^2)/Sbase

begin
    #include("C:/Users/liemann/github/PowerDynamicsDevelopment/src/include_costum_nodes_lines_utilities.jl")
    #include("C:/Users/liemann/github/PowerDynamicsDevelopment/src/grids/LTVS_Test_System.jl")
    include("C:/Users/liemann/github/PowerDynamicsDevelopment/src/grids/LTVS_Test_System_TapParam.jl")
    #include("C:/Users/liemann/github/PowerDynamicsDevelopment/src/sensitivity_analyses/Local_Sensitivity.jl")
end
begin
    pg = GFC_LTVS_Test_SystemTapParam(nTap = 5)
    Qmax   = [Inf, Inf, Inf,Inf, Inf,Inf*sqrt(1-0.95^2)]
    Qmin   = -Qmax
    U,δ,ic0 = PowerFlowClassic(pg,iwamoto = true, Qmax = Qmax, Qmin = Qmin, Qlimit_iter_check = 2,max_tol = 1e-6)
    Uc = U.*exp.(1im*δ/180*pi)
    Ykk = NodalAdmittanceMatrice(pg)
    Ic = abs.(Ykk*Uc./5.5)
    S  = round.(Uc.*(conj.(Ykk)*conj.(Uc)),digits=3)
    pg, ic0 = InitializeInternalDynamics(pg,ic0)

    #pg, ic0 = GetInitializedLTVSSystem()
    pgsol,evr  = run_LTVS_simulationTapParam(pg,ic0,(0.0,120.0))
    display(plot!(pgsol,collect(keys(pg.nodes))[2:end-1],:v,legend = false, linestyle = :dash))
end
display(plot!(pgsol,collect(keys(pg.nodes))[2:end-1],:v,ylims=(0.7,1.01),xlims=(0.0,90), label =false))
vline!([67.5],linestyle = :dash,color = "black",label =false)
pg= GFC_LTVS_Test_System()
dimension(pg.nodes["bus3"])
rhs(pg).syms[10]

plot(pgsol,"bus4",:i_abs)
plo!t(pgsol,collect(keys(pg.nodes)),:v,legend = false,xlim=(60.0,110.0))
plot(pgsol,"bus4",:i_abs, legend = false, ylims =(0.95,1.1))
plot(pgsol,"bus4",:ω, legend = (0.8,0.1))
plot(pgsol,"bus4",:θ, legend = (0.8,0.8),ylims=(-0.1,0.15))
plot!(pgsol_droop,"bus4",:θ, legend = (0.8,0.8),ylims=(-0.05,0.05))
plot(pgsol,"bus4",:Pout, legend = (0.8,0.5))
p = ExtractResult(pgsol,"bus4",:Pout)
q = ExtractResult(pgsol,"bus4",:Qout)
s = hypot.(p,q)
plot(pgsol.dqsol.t,s, xlims=(0,2))
plot(0.0,4,:int,3)

variable_index(pgsol.powergrid.nodes, "bus4", 1)

mtk_normal = GetMTKLTVSSystemTapParam(pg_state = "gfc_normal")
mtk_fault = GetMTKLTVSSystemTapParam(pg_state = "gfc_fault")
mtk_postfault = GetMTKLTVSSystemTapParam(pg_state = "gfc_postfault")
mtk = [mtk_normal; mtk_fault; mtk_postfault]

s = GetTriggCondsLTVS(mtk_normal)
h = GetStateResFunLTVS(mtk_normal)
p_pre = GFC_LTVS_params_TapParam()
sensis_p = collect(1:16)
@time toll_tap = CalcHybridTrajectorySensitivity(mtk,pgsol.dqsol,p_pre,evr,s,h,[],sensis_p)

save("C:/Users/liemann/Desktop/Sens_LTVS/sens_pscc_kq_1em1_dt_1em2_tap_param.jld", "sens", toll_tap,"ic0",ic0,"p_pre",p_pre,"evr",evr,"sensis_p",sensis_p)


toll = load("C:/Users/liemann/Desktop/Sens_LTVS/sens_kq_1em3_t_90_dt_1em2.jld")
toll = toll["sens"]
plot(1:length(toll[1][16,1:end]),toll[1][16,1:end])
plot!(pgsol.dqsol.t[1:end-1],toll[1][16,1:end])
xlims!((1,900))
xlims!((1,20))
ylims!((0.8,1.005))
plot!(pgsol,"bus4",:v, label = "Imax = 1.0") #, linestyle = :dash)
a = [rhs(pg).syms sol.u[end] ic0]
pgsol = PowerGridSolution(sol,pg)
collect(1:15)

plot!(pgsol,collect(keys(pg.nodes)),:v,legend = false) #, linestyle = :dash
xlims!((0,90))#, linestyle = :dash
plot(pgsol,"bus4",:i_abs, label = "original")
xlims!((0.9,5.3))
ylims!((0.99,1.01))
plot(pgsol,"bus4",:Pout)
plot!(pgsol,"bus4",:Qout)
plot(pgsol,"bus4",:v)

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
    "imax_csa",#16
    ]
syms = rhs(pg).syms
look_on = 16
t_cut = 1200*0
t_cut2 = 300*0
plot(pgsol.dqsol.t[1:end-1],toll[1][look_on,1:end], title = "Sensis of $(String(syms[look_on]))",
    label = labels_p[1],
    legend = :outertopright,
    size = (1000,750),
    #xlims = (1.5,65.0),
    xlims = (0.8,1.5),
    ylims = (-20,100))
    #ylims = (-1.5,1.5))
for i in sensis_p[2:end]
    display(plot!(pgsol.dqsol.t[1:end-1],toll[i][look_on,1:end], label = labels_p[i]))
    #sleep(3.0)
end

xlims!((0,10))
ylims!((-1 ,2))
toll = toll_new
begin
    indi = 14 #15 ist am interessantesten!!!!
    display(plot(pgsol.dqsol.t[1:end-1],toll[indi][look_on,1:end], label = labels_p[indi], title = "Sensis of $(String(syms[look_on]))"))
    #xlims!((0.0,5.0))
end

#Calculation of approximated solution
sol_sensi_per = deepcopy(pgsol.dqsol)
for (ind,val) in enumerate(collect(eachcol(toll[1])))# - 39019
    sol_sensi_per.u[ind] .+= val*(0.004)
    #display(val)
end
pgsol_tmp = PowerGridSolution(sol_sensi_per,pg)
plot!(pgsol_tmp,"bus4",:i_abs, label = "approximated", legend = (0.2,0.5))
plot!(pgsol,"bus4",:i_abs, label = "real perturbed")
plot!(pgsol_tmp,collect(keys(pg.nodes)),:v,legend = false, linestyle = :dash) #


# Plotting influence on voltage
sol_sensi_per = deepcopy(pgsol.dqsol)
for (ind,val) in enumerate(collect(eachcol(toll[1])))# - 39019
    sol_sensi_per.u[ind] .+= val*(0.009)
    #display(val)
end
pgsol_tmp = PowerGridSolution(sol_sensi_per,pg)
plot(pgsol_tmp,"bus2",:v, label = labels_p[1], title = "Influence on V2 Betrag", legend = :outertopright, size = (1000,1000))
plot!(pgsol_tmp,collect(keys(pg.nodes)),:v,legend = false, linestyle = :dash)

for i in 3:length(toll)
    sol_sensi_per = deepcopy(pgsol.dqsol)
    for (ind,val) in enumerate(collect(eachcol(toll[i])))# - 39019
        sol_sensi_per.u[ind] .+= val*(0.01)
        #display(val)
    end
    pgsol_tmp = PowerGridSolution(sol_sensi_per,pg)
    display(plot!(pgsol_tmp,"bus2",:v, label = labels_p[i]))
end

xlims!((10,90))
ylims!((0.86 ,0.94))
ylims!((0.79 ,0.82))

display(plot(pgsol.dqsol.t[1:end-1],toll[1][1,1:end]))
for i in 2:18
    display(plot!(pgsol.dqsol.t[1:end-1],toll[1][i,1:end])) #, label = String(rhs(pg).syms[i])
end
