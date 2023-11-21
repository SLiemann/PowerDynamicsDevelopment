using PlotlyJS, DataFrames
using ModelingToolkit   

Sbase = 8000e6
Ubase = s
Zbase = Ubase^2/Sbase

begin
    include("C:/Users/liemann/github/PowerDynamicsDevelopment/src/include_custom_nodes_lines_utilities.jl")
    include("C:/Users/liemann/github/PowerDynamicsDevelopment/src/grids/LTVS_Test_System_N32 PEL_TS.jl")
    nothing
end

pg0,ic0 = Initialize_N32_PEL_TS();
@time pgsol0, evr_sol = simulate_LTVS_N32_simulation_TS(pg0,ic0,(0.0,0.6),(20.0+1im*20)/Zbase);
plotallvoltages(pgsol0);
plot(myplot(pgsol0,"bus_gfm",:q_imax))
plot(myplot(pgsol0,"bus_gfm",:q_idcmax))
plot(myplot(pgsol0,"bus_gfm",:idc0))
plot(myplot(pgsol0,"bus_gfm",:udc))
