using PlotlyJS, DataFrames
using ModelingToolkit   

Sbase = 8000e6
Ubase = 400e3
Zbase = Ubase^2/Sbase

begin
    include("C:/Users/liemann/github/PowerDynamicsDevelopment/src/include_custom_nodes_lines_utilities.jl")
    include("C:/Users/liemann/github/PowerDynamicsDevelopment/src/grids/LTVS_Test_System_N32 GFM_PEL_TS.jl")
    nothing
end

pg0,ic0 = Initialize_N32_GFM_PEL_TS();
@time pgsol0_per, evr_sol = simulate_LTVS_N32_simulation_N32_GFM_PEL_TS(pg0,ic0,(0.0,0.6),(100.0)/Zbase);
plotallvoltages(pgsol0_per);
plot(myplot(pgsol0,"bus_gfm",:q_imax))
plot(myplot(pgsol0,"bus_gfm",:q_idcmax))
plot(myplot(pgsol0,"bus_gfm",:idc0))
plot([myplot(pgsol0,"bus_gfm",:θ),myplot(pgsol0_per,"bus_gfm",:θ)])

plot([myplot(pgsol0,"bus_gfm",:P0),myplot(pgsol0_per,"bus_gfm",:P0)])

