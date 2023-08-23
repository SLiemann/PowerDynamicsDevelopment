using PlotlyJS, DataFrames, MAT

Sbase = 8000e6
Ubase = 400e3
Zbase = Ubase^2/Sbase

begin
    include("C:/Users/liemann/github/PowerDynamicsDevelopment/src/include_custom_nodes_lines_utilities.jl")
    include("C:/Users/liemann/github/PowerDynamicsDevelopment/src/grids/LTVS_Test_System_N32 GFM.jl")
    nothing
end

pg0,ic0 = Initialize_N32_GFM(1,1);
@time pgsol0, suc0,FRT0 = simulate_LTVS_N32_simulation(pg0,ic0,(0.0,160.0),(1e6)/Zbase);
plot(plotallvoltages(pgsol0))
plot(myplot(pgsol0,"bus_gfm",:P0))
plot(myplot(pgsol0,"bus_gfm",:Q0))
plot(myplot(pgsol0,"bus_gfm",:Pdelta))
plot(myplot(pgsol0,"bus_gfm",:i_abs))


dir = "\\\\fs0\\home\\liemann\\"
file = matopen(dir*"vltvs_200s.mat")
file = matopen(dir*"vltvs_short.mat")
file = matopen(dir*"v_droop_secm4.mat")
#file = matopen(dir*"1em6_Rf_1e6_vltvs2.mat")
#file = matopen(dir*"awu1_100s_vltvs.mat")
vltvs = read(file, "Vltvs");
close(file)

p1 = plotallvoltages(pgsol0);
append!(p1,[scatter(x=vltvs[:,1],y=vltvs[:,2],name="Mat_Bus1")])
append!(p1,[scatter(x=vltvs[:,1],y=vltvs[:,3],name="Mat_EHV")])
append!(p1,[scatter(x=vltvs[:,1],y=vltvs[:,4],name="Mat_HV")])
append!(p1,[scatter(x=vltvs[:,1],y=vltvs[:,5],name="Mat_load")])
append!(p1,[scatter(x=vltvs[:,1],y=vltvs[:,6],name="Mat_GFM")])
append!(p1,[scatter(x=[0.0,32.0],y=[0.99,0.99],line="dash")])
plot(p1)

dir = "\\\\fs0\\home\\liemann\\"
file = matopen(dir*"pmat.mat")
pq = read(file, "Pgfm");
close(file)

p1 = [myplot(pgsol0,"bus_gfm",:P0)]
append!(p1,[scatter(x=pq[:,1],y=pq[:,2],name="Mat_P")])
plot(p1)

file = matopen(dir*"qmat.mat")
pq = read(file, "Qgfm");
close(file)
p1 = [myplot(pgsol0,"bus_gfm",:Q0)]
append!(p1,[scatter(x=pq[:,1],y=pq[:,2],name="Mat_Q")])
plot(p1)

p1 = [myplot(pgsol0,"bus_load",:P0)]
append!(p1,[scatter(x=pq[:,1],y=-pq[:,4],name="Mat_P")])
plot(p1)
p1 = [myplot(pgsol0,"bus_load",:Q0)]
append!(p1,[scatter(x=pq[:,1],y=-pq[:,5],name="Mat_Q")])
plot(p1)

dir = "\\\\fs0\\home\\liemann\\"
file = matopen(dir*"iabs.mat")
pq = read(file, "iabs");
close(file)
p1 = [myplot(pgsol0,"bus_gfm",:i_abs)]
append!(p1,[scatter(x=pq[:,1],y=pq[:,2],name="Mat_Iabs")])
plot(p1)

begin
    pg = LTVS_Test_System_N32_GFM(gfm=1,awu=0)
    s0 = pg.nodes["bus_gfm"].Srated/8000e6
    ps = pg.nodes["bus_gfm"].p0set/s0
    Qmax   = [Inf,Inf, Inf, Inf,Inf,Inf,5300/8000*sqrt(1-ps^2)]
    Qmin   = -Qmax
    U1,δ1,ic0,cu = PowerFlowClassic(pg,iwamoto = false,max_tol = 1e-4,iter_max = 100,Qmax = Qmax, Qmin = Qmin,Qlimit_iter_check=80)
end