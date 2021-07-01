using Plots
using PowerDynamics
using DifferentialEquations
using ModelingToolkit


begin
    include("C:/Users/liemann/github/PowerDynamicsDevelopment/src/include_costum_nodes_lines_utilities.jl")
    include("C:/Users/liemann/github/PowerDynamicsDevelopment/src/grids/LTVS_Test_System.jl")
    include("C:/Users/liemann/github/PowerDynamicsDevelopment/src/sensitivity_analyses/Local_Sensitivity.jl")
    #include("C:/Users/liemann/github/PowerDynamicsDevelopment/src/sensitivity_analyses/LS_old.jl")
end

pg = GFC_LTVS_Test_System()
pgpf = GetPostFaultLTVSPG(pg)

p = GFC_LTVS_params()
ew = CalcEigenValues(pgpf,p; output = true, plot =true)
xlims!((-2.5,0.25))
ylims!((-10,10))
