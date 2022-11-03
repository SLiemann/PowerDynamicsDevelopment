using PlotlyJS, DataFrames
begin
    include("C:/Users/liemann/github/PowerDynamicsDevelopment/src/include_custom_nodes_lines_utilities.jl")
    include("C:/Users/liemann/github/PowerDynamicsDevelopment/src/grids/LTVS_Test_System_N32 GFM.jl")
end

pgsol = run_LTVS_N32_simulation(1,(0.0,15.0));

rhs(pgsol.powergrid).syms
plotallvoltages(pgsol);

tmp = myplot(pgsol,"bus_gfm",:u_r)
tmp2 = myplot(pgsol,"bus_gfm",:u_i)
plot([tmp,tmp2])

plot([plot(pgsol,:Um_7),plotv(pgsol,"bus_gfm")])
plotv(pgsol,["bus_gfm","bus_load"])


