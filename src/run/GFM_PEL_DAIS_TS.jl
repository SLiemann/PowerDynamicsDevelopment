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
@time pgsol0, evr_sol = simulate_LTVS_N32_simulation_N32_GFM_PEL_TS(pg0,ic0,(0.0,10.0),(20.0+1im*20)/Zbase);
jp = plotallvoltages(pgsol0);
plot(myplot(pgsol0,"bus_load",:toff))
plot(myplot(pgsol0,"bus_load",:ton))
plot(myplot(pgsol0,"bus_load",:q_on))
plot(myplot(pgsol0,"bus_load",:vofft2))
plot(myplot(pgsol0,"bus_load",:Vabstoff))
plot(myplot(pgsol0,"bus_load",:p1))
plot(myplot(pgsol0,"bus_load",:q1))

collect((rhs(pg0).syms))[33]

plot(myplot(pgsol0,"bus_gfm",:q_imax))
plot(myplot(pgsol0,"bus_gfm",:q_idcmax))
plot(myplot(pgsol0,"bus_gfm",:i_abs))
plot(myplot(pgsol0,"bus_gfm",:idc0))
myplot(pgsol0,"bus_gfm",:P0)
plot([myplot(pgsol0,"bus_gfm",:θ),myplot(pgsol0_per,"bus_gfm",:θ)])

plot([myplot(pgsol0,"bus_gfm",:P0),myplot(pgsol0_per,"bus_gfm",:P0)])

using MAT
datapath = "\\\\fs0\\home\\liemann\\Diss\\Results\\623_GFM_PEL\\"
mdata = matread(datapath*"PEL_0_PEL_share_0.3_Zf_20.mat")
vpel = mdata["Vpel"]
mp = [
        PlotlyJS.scatter(x=vpel[:,1],y=vpel[:,2],name="V0 (EMT)"),
        PlotlyJS.scatter(x=vpel[:,1],y=vpel[:,3],name="V1 (EMT)"),
        PlotlyJS.scatter(x=vpel[:,1],y=vpel[:,4],name="V2 (EMT)"),
        PlotlyJS.scatter(x=vpel[:,1],y=vpel[:,5],name="V_GFM (EMT)"),
        PlotlyJS.scatter(x=vpel[:,1],y=vpel[:,6],name="V_load (EMT)")
    ]
append!(mp,jp)
PlotlyJS.plot(mp)