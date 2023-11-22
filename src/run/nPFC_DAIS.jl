using PlotlyJS, Plots, DataFrames
using ModelingToolkit   

Sbase = 8000e6
Ubase = 400e3
Zbase = Ubase^2/Sbase

begin
    include("C:/Users/liemann/github/PowerDynamicsDevelopment/src/include_custom_nodes_lines_utilities.jl")
    include("C:/Users/liemann/github/PowerDynamicsDevelopment/src/grids/LTVS_Test_System_N32 PEL_TS.jl")
    nothing
end

pg0,ic0 = Initialize_N32_PEL_TS();
rhs(pg0).syms .=> ic0

@time pgsol0, evr_sol = simulate_LTVS_N32_simulation_PEL_TS(pg0,ic0,(0.0,0.3),(80.0+1im*0)/Zbase);
plotallvoltages(pgsol0);
myplot(pgsol0,"bus_load",:q1);
myplot(pgsol0,"bus_load",:p1);
myplot(pgsol0,"bus_load",:tsum);
myplot(pgsol0,"bus_load",:ton);
myplot(pgsol0,"bus_load",:toff);
myplot(pgsol0,"bus_load",:vofft2);

pgsol0.dqsol.t[end]
length(pgsol0.dqsol)
plotv(pgsol0,"bus_load")



