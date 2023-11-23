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

pg0,ic0 = Initialize_N32_PEL_TS(share_pe= 0.10);
rhs(pg0).syms .=> ic0

@time pgsol0, evr_sol = simulate_LTVS_N32_simulation_PEL_TS(pg0,ic0,(0.0,0.201),(20.0+1im*20)/Zbase);
ppel2 = plotallvoltages(pgsol0);
PlotlyJS.plot([ppel;ppel2])

myplot(pgsol0,"bus_load",:q1);
myplot(pgsol0,"bus_load",:p1);
myplot(pgsol0,"bus_load",:tsum);
myplot(pgsol0,"bus_load",:ton);
myplot(pgsol0,"bus_load",:toff);
myplot(pgsol0,"bus_load",:vofft2);

myplot(pgsol0,"bus_load",[:tsum,:ton,:toff]);

p1 = myplot(pgsol0,"bus_load",:ton,y_norm=1/50);
p2 = myplot(pgsol0,"bus_load",:toff,y_norm=1/500,y_bias=-2.32);
p3 = myplot(pgsol0,"bus_load",:p1,y_norm=-1);
p4= myplot(pgsol0,"bus_load",:vofft2,y_norm=10);
p5 = plotv(pgsol0,"bus_load",y_norm=5)
PlotlyJS.plot([p1,p2,p3,p4,p5])


pgsol0.dqsol.t[end]
length(pgsol0.dqsol)
p5 = plotv(pgsol0,"bus_load")



