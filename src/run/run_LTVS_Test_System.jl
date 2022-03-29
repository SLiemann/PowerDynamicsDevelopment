using PowerDynamics
#using OrderedCollections: OrderedDict
using Plots
#import PowerDynamics: PiModel
using DifferentialEquations
#using CSV #read PF DataFrames
#using DataFrames #for CSV
#using Distributed
using JLD

#Ubase = 380e3
#Sbase = 100e6
#Zbase = (Ubase^2)/Sbase

begin
    include("C:/Users/liemann/github/PowerDynamicsDevelopment/src/include_costum_nodes_lines_utilities.jl")
    #include("C:/Users/liemann/github/PowerDynamicsDevelopment/src/grids/LTVS_Test_System.jl")
    #include("C:/Users/liemann/github/PowerDynamicsDevelopment/src/grids/LTVS_Test_System_TapParam.jl")
    include("C:/Users/liemann/github/PowerDynamicsDevelopment/src/sensitivity_analyses/Local_Sensitivity.jl")
    include("C:/Users/liemann/github/PowerDynamicsDevelopment/src/utility/utility_functions.jl")
end
begin
    include("C:/Users/liemann/github/PowerDynamicsDevelopment/src/grids/LTVS_Test_System_TapParam.jl")
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
    pgsol_per,evr  = run_LTVS_simulationTapParam(pg,ic0,(0.9,120.0))
    #plot(pgsol,"bus4",:i_abs,label = "I-Original",xlims=(5,70),ylims=(1.02,1.11), legend = (0.5,0.1))
    #display(plot!(pgsol_per,"bus4",:i_abs, label ="real perturbed"))
    #display(plot(pgsol_stkvi,"bus4",:v,label = "Enhanced"))
    #display(plot!(pgsol_stkvq,"bus4",:v,label = "Enh kvq"))
    #display(plot(pgsol,"bus4",:v,label = "U-Original",xlims=(5,70),ylims=(0.9,1.01), legend = (0.5,0.1)))
    #display(plot!(pgsol_per,"bus4",:v,label = "real perturbed"))
    #display(plot!(pgsol_per,"bus4",:v, label ="real perturbed")) #linestyle = :dash
end
x = CalcEigenValues(pg,GFC_LTVS_params_TapParam(),output=true, plot=true)


sensi = load("C:/Users/liemann/Desktop/DIESE_sens_pscc_kq_1em1_dt_1em2_tap_param.jld")
toll_tap = sensi["sens"]
PlotApproTrajectories(pg,pgsol,pgsol_per,sensi["sens"],15,0.1,0.15,labels_p,:i_abs)
ylims!(0.95,0.99)
xlims!(-5,0)
xlims!(0,3)
ylims!(-2,2)
ylims!(0,280)

plot(pgsol,"bus4",:i_abs)
plot!(pgsol_per,"bus4",:i_abs, linestyle = :dash)
plot!(pgsol_per,"bus4",:v, label="voltage", title="Long-term voltage stability", grid=true,gridstyle = :dash,gridalpha = 0.5,linewidth = 2)
plot(pgsol,"bus4",:i_abs, legend = false, ylims =(0.95,1.1))
plot(pgsol,"bus4",:ω, legend = (0.8,0.1))
plot(pgsol,"bus4",:θ, legend = (0.8,0.8),ylims=(-0.1,0.15))
plot!(pgsol_droop,"bus4",:θ, legend = (0.8,0.8),ylims=(-0.05,0.05))
plot(pgsol,"bus4",:Pout, legend = (0.8,0.5))
plot(pgsol,"bus4",:Qout, legend = (0.8,0.5))
vi = ExtractResult(pgsol,:u_i_4)
vr = ExtractResult(pgsol,:u_r_4)
v_abs = sqrt.(vi.^2+vr.^2)
i_abs = ExtractResult(pgsol,:i_abs_4)
time = pgsol.dqsol.t
f = plot(time,[v_abs i_abs],layout=(1,2),label=["\$|u_f|\$" "\$|i_c|\$"], color = "green", size=(800,350),legend=:bottomright,
    ytickfont = font(11, "Times"),
    xtickfont = font(11, "Times"),
    legendfont= font(11,"Times"),
    grid = true,
    gridalpha = 0.5,
    gridstyle = :dash,
    framestyle = :box,
    linewidth = 1.5,
    ylims = (0.6,1.2),)
savefig(f,"Test.png")


plot!(pgsol_per,"bus4",:Pout, legend = (0.8,0.5))

mtk_normal = GetMTKLTVSSystemTapParam(pg_state = "gfc_normal")
mtk_fault = GetMTKLTVSSystemTapParam(pg_state = "gfc_fault")
mtk_postfault = GetMTKLTVSSystemTapParam(pg_state = "gfc_postfault")
mtk = [mtk_normal; mtk_fault; mtk_postfault]

s = GetTriggCondsLTVS(mtk_normal)
h = GetStateResFunLTVS(mtk_normal)
p_pre = GFC_LTVS_params_TapParam()
sensis_p = collect(1:16)
@time toll_tap = CalcHybridTrajectorySensitivity(mtk,pgsol.dqsol,p_pre,evr,s,h,[],sensis_p)

save("C:/Users/liemann/Desktop/Sens_LTVS/sens_short_pscc_kq_1em1_dt_1em3_tap_param.jld", "sens", toll_tap,"ic0",ic0,"p_pre",p_pre,"evr",evr,"sensis_p",sensis_p)

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

plot(pgsol.dqsol.t[1:end-1],toll_tap[1][look_on,1:end], title = "Sensis of $(String(syms[look_on]))",
    label = labels_p[1],
    legend = :outertopright,
    size = (1000,750),
    #xlims = (1.5,65.0),
    #xlims = (0.9,1.5),
    ylims = (-30,100))
    #ylims = (-1.5,1.5))
for i in sensis_p[2:end]
    display(plot!(pgsol.dqsol.t[1:end-1],toll_tap[i][look_on,1:end], label = labels_p[i]))
    #sleep(3.0)
end

xlims!((0.9,1.2))
ylims!((-3 ,1.5))
toll = toll_new
begin
    indi = 14 #15 ist am interessantesten!!!!
    display(plot(pgsol.dqsol.t[1:end-1],toll[indi][look_on,1:end], label = labels_p[indi], title = "Sensis of $(String(syms[look_on]))"))
end


sens_ur4,sens_ui4 = GetVoltageSensis(toll_tap,7,8);

params = GFC_LTVS_params_TapParam()
plot(pgsol.dqsol.t[1:end-1],sens_ur4[:,1].*params[1])
for i in collect(2:16)
    display(plot!(pgsol.dqsol.t[1:end-1],sens_ur4[:,i].*params[i]))
end

plot(pgsol.dqsol.t[1:end-1],sens_ui4[:,1].*params[1])
for i in collect(2:16)
    display(plot!(pgsol.dqsol.t[1:end-1],sens_ui4[:,i].*params[i]))
end

plot(pgsol.dqsol.t[1:end-1],sqrt.(sens_ui4[:,1])

for i in collect(2:16)
    display(plot!(pgsol.dqsol.t[1:end-1],sens_ui4[:,i].*params[i]))
end
