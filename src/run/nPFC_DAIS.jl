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
# unstable case: share_pe = 0.3; Rf = 10, Xf = 0
pg0,ic0 = Initialize_N32_PEL_TS(share_pe= 0.30);
rhs(pg0).syms .=> ic0

@time pgsol0, evr_sol = simulate_LTVS_N32_simulation_PEL_TS(pg0,ic0,(0.0,0.6),(20.0+1im*20)/Zbase);
myplot(pgsol0,"bus_load",:q_on);
ppel3 = plotallvoltages(pgsol0);
PlotlyJS.plot([ppel2;ppel3])

pp2 = myplot(pgsol0,"bus_load",:q1)
pq2 = myplot(pgsol0,"bus_load",:p1);
p1 = myplot(pgsol0,"bus_load",:tsum);
myplot(pgsol0,"bus_load",:ton);
myplot(pgsol0,"bus_load",:toff);
myplot(pgsol0,"bus_load",:vofft2);
myplot(pgsol0,"bus_load",:Vabstoff);
myplot(pgsol0,"bus_load",:q_on);

myplot(pgsol0,"bus_load",[:tsum,:ton,:toff]);

p1 = myplot(pgsol0,"bus_load",:ton,y_norm=1/50);
p2 = myplot(pgsol0,"bus_load",:toff,y_norm=1/500,y_bias=-2.32);
p3 = myplot(pgsol0,"bus_load",:p1,y_norm=-1);
p4= myplot(pgsol0,"bus_load",:vofft2,y_norm=10);
p5 = plotv(pgsol0,"bus_load",y_norm=5)
PlotlyJS.plot([p1,p2,p3,p4,p5])
