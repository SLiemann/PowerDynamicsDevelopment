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
    include("C:/Users/liemann/github/PowerDynamicsDevelopment/src/include_costum_nodes_lines_utilities.jl")
    include("C:/Users/liemann/github/PowerDynamicsDevelopment/src/grids/LTVS_Test_System.jl")
    include("C:/Users/liemann/github/PowerDynamicsDevelopment/src/sensitivity_analyses/Local_Sensitivity.jl")
    #include("C:/Users/liemann/github/PowerDynamicsDevelopment/src/sensitivity_analyses/LS_old.jl")
end

pg, ic0 = GetInitializedLTVSSystem(gfc = "gfc_normal")
pgsol,evr  = run_LTVS_simulation(pg,ic0,(0.0,5.0))

mtk_normal = GetMTKLTVSSystem(pg_state = "gfc_normal")
mtk_fault = GetMTKLTVSSystem(pg_state = "gfc_fault")
mtk_postfault = GetMTKLTVSSystem(pg_state = "gfc_postfault")
mtk = [mtk_normal; mtk_fault; mtk_postfault]
mtk = [mtk_normal; mtk_normal; mtk_normal]

s = GetTriggCondsLTVS(mtk_normal)
h = GetStateResFunLTVS(mtk_normal)
p_pre = GFC_LTVS_params()
sensis_p = collect(1:15)
@time toll_new = CalcHybridTrajectorySensitivity(mtk,pgsol.dqsol,p_pre,evr,s,h,[],sensis_p)

save("C:/Users/liemann/Desktop/Sens_LTVS/sens_kq_on_90_1em2.jld", "sens", toll,"ic0",ic0,"p_pre",p_pre,"evr",evr,"sensis_p",sensis_p)



toll2 = load("/tmp/myfile.jld")
plot(pgsol.dqsol.t[1:end-1],toll[1][16,1:end])
xlims!((1.99,2.2))
xlims!((2.1,3.0))
ylims!((-0.3,0.60))

a = [rhs(pg).syms sol.u[end] ic0]
pgsol = PowerGridSolution(sol,pg)
collect(1:15)
plot(pgsol,collect(keys(pg.nodes)),:v,legend = (0.8, 0.5))
plot!(pgsol,"bus4",:i_abs)
xlims!((1.9,2.3))
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
    ]
syms = rhs(pg).syms
look_on = 16
plot(pgsol.dqsol.t[1:end-1],toll[1][look_on,1:end], title = "Sensis of $(String(syms[look_on]))",
    label = labels_p[1],
    legend = :outertopright,
    size = (1000,750))
    xlims!((1.9,2.3))
for i in sensis_p[3:end]
    display(plot!(pgsol.dqsol.t[1:end-1],toll[i][look_on,1:end], label = labels_p[i]))
    #sleep(3.0)
end
xlims!((1.9,2.3))
ylims!((-40.0,10.0))
begin
    indi = 15
    display(plot(pgsol.dqsol.t[1:end-1],toll[indi][look_on,1:end], label = labels_p[indi], title = "Sensis of $(String(syms[look_on]))"))
    xlims!((0.0,5.0))
end


sol_sensi_per = deepcopy(pgsol.dqsol)
for (ind,val) in enumerate(collect(eachcol(toll[15])))# - 39019
    sol_sensi_per.u[ind] .+= val*(0.001-0.01)
    #display(val)
end
pgsol_tmp = PowerGridSolution(sol_sensi_per,pg)
plot!(pgsol,"bus4",:i_abs)
plot!(pgsol_tmp,"bus4",:i_abs)



plot(pgsol,"bus4", :Um,size = (1000, 500),label = "PD-Um")
begin
    plot(pgsol,collect(keys(pg.nodes))[1:end-1], :v,size = (1000, 500))
    #plot(pgsol,"bus4", :Um,size = (1000, 500),label = "PD-Um")
    #plot(pgsol,"bus4", :v,label = "PD-V0", legend= false)
    #ylims!((0.90,1.01))
    xlims!((1.9,2.3))
    test = DataFrame(CSV.File("C:\\Users\\liemann\\Desktop\\Basisszenario\\data.csv"; header=false, delim=';', type=Float64))
    plot!(test.Column1,test.Column2,label = "PF-bus1")
    plot!(test.Column1,test.Column3,label = "PF-bus2")
    plot!(test.Column1,test.Column4,label = "PF-bus3")
    plot!(test.Column1,test.Column5,label = "Matlab-V0")
    #plot!(test.Column1,test.Column6,label = "Matlab-Um")
end

begin
    #plot(pgsol,"bus4", :iabs,size = (1000, 500),label = "PD-I0")
    plot(pgsol,"bus4", :i_abs,label = "PD-Idq")
    #ylims!((0.95,1.05))
    xlims!((1.9,2.2))
    test = DataFrame(CSV.File("C:\\Users\\liemann\\Desktop\\Basisszenario\\data_c.csv"; header=false, delim=';', type=Float64))
    plot!(test.Column1,test.Column2,label = "Matlab-Idq")
    plot!(test.Column1,test.Column3,label = "Matlab-I0")
end

begin
    plot(pgsol,"bus4", :Ixcf,size = (1000, 500),label = "PD-Ixcf")
    ylims!((0.06,0.07))
    #xlims!((0,2.0))
    test = DataFrame(CSV.File("C:\\Users\\liemann\\Desktop\\Basisszenario\\data_c.csv"; header=false, delim=';', type=Float64))
    plot!(test.Column1,test.Column4,label = "Matlab-Icf")
end
