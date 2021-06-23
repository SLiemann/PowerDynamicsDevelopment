using PowerDynamics
#using OrderedCollections: OrderedDict
using Plots
#import PowerDynamics: PiModel
using DifferentialEquations
#using CSV #read PF DataFrames
#using DataFrames #for CSV
#using Distributed
#using JLD

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

s = GetTriggCondsLTVS(mtk_normal)
h = GetStateResFunLTVS(mtk_normal)
p_pre = GFC_LTVS_params()
@time toll = CalcHybridTrajectorySensitivity(mtk,pgsol.dqsol,p_pre,evr,s,h,[1],[13,14])

toll2 = load("/tmp/myfile.jld")
plot(pgsol.dqsol.t[1:end-1],toll[3][16,1:end])

a = [rhs(pg).syms sol.u[end] ic0]
pgsol = PowerGridSolution(sol,pg)

plot(pgsol,collect(keys(pg.nodes)),:v,legend = (0.8, 0.5))
plot(pgsol,"bus4",:i_abs)
xlims!((1.9,2.3))
ylims!((0.99,1.01))
plot(pgsol,"bus4",:Pout)
plot!(pgsol,"bus4",:Qout)
plot(pgsol,"bus4",:v)

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
