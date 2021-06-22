using PowerDynamics
using DifferentialEquations
using Plots
using ModelingToolkit
using CSV
using DataFrames

begin
    include("C:/Users/liemann/github/PowerDynamicsDevelopment/src/include_costum_nodes_lines_utilities.jl")
    include("C:/Users/liemann/github/PowerDynamicsDevelopment/src/grids/OLTC_Hybrid_Sensis.jl")
    include("C:/Users/liemann/github/PowerDynamicsDevelopment/src/sensitivity_analyses/Local_Sensitivity.jl")
end

#simulate nominal trajectory
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

#Calc hybrid sensis
mtk = GetMTKOLTCSystem()
p_pre = GetParametersOLTCHisken(5.0)
s = GetTriggCondsOLTCHisken(mtk)
h = GetStateResFunOLTCHisken(mtk)
@time toll = CalcHybridTrajectorySensitivity(mtk,sol,p_pre,evr,s,h,[1],[3])
plot(pgsol.dqsol.t[1:end-1],toll[2][8,1:end]) #
ylims!((-0.008,0.024)) #.+toll[2][8,1:end-1]
xlims!((30.0099,30.0101))
yticks!(collect(-0.008:0.002:0.024))
xticks!(collect(0.0:20.:200.))

#Add sensis
sol_sensi_per = deepcopy(sol)
for (ind,val) in enumerate(collect(eachcol(toll[2])))# - 39019
    sol_sensi_per.u[ind] .+= val*(-1.5)
    #display(val)
end
pgsol_sensi_per = PowerGridSolution(sol_sensi_per,pgsol.powergrid)
plot(pgsol,["bus4"],:v,label="nominal")
plot!(pgsol_sensi_per,["bus4"],:v,label="approximated", linestyle = :dash)

pgsol_per,evr = SimulateOLTCHIsken(Tp = 3.5)
plot!(pgsol_per,["bus4"],:v,label="perturbed")
