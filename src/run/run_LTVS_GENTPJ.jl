using PlotlyJS, DataFrames

Sbase = 8000e6
Ubase = 400e3
Zbase = Ubase^2/Sbase

begin
    include("C:/Users/liemann/github/PowerDynamicsDevelopment/src/include_custom_nodes_lines_utilities.jl")
    include("C:/Users/liemann/github/PowerDynamicsDevelopment/src/grids/LTVS_Test_System_N32 GENTPJ.jl")
    nothing
end

pg0,ic0 = Initialize_N32_GENTPJ();
@time pgsol0, suc0,FRT0 = simulate_LTVS_N32_simulation_GENTPJ(pg0,ic0,(0.0,3.0),(20.0+1im*20)/Zbase);
plot(plotallvoltages(pgsol0))
plot([myplot(pgsol0,"bus_gfm",:LVRT),plotv(pgsol0,["bus_gfm"])[1]])

plot(myplot(pgsol0,"bus_gfm",:i_abs))
plot(myplot(pgsol0,"bus_gfm",:udc))
plot(myplot(pgsol0,"bus_gfm",:Q0))
plot(myplot(pgsol0,"bus_gfm",:Pf))
plotallvoltages(pgsol0)
plot(myplot(pgsol0,"bus_gfm",:θ))
plot(myplot(pgsol0,"bus_gfm",:θ,y_norm=pi/180,y_bias=18.5807))


using MAT
dir = "C:\\Users\\liemann\\Desktop\\GENTPJ_AVR_OEL\\"
file = matopen(dir*"Base_Case_3s.mat")
vltvs = read(file, "Vsimemt_pe");
close(file)

begin
    p1 = plotallvoltages(pgsol0)
    for i = 1:5
        p2 = scatter(x=vltvs[:,6],y=vltvs[:,i])
        push!(p1,p2)
    end
end
plot(p1)
