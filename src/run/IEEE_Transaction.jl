using PlotlyJS, DataFrames

begin
    include("C:/Users/liemann/github/PowerDynamicsDevelopment/src/include_custom_nodes_lines_utilities.jl")
    include("C:/Users/liemann/github/PowerDynamicsDevelopment/src/grids/LTVS_Test_System_N32 GFM.jl")
    nothing
end

pgsol = run_LTVS_N32_simulation(4,0,(0.0,2.0));

plot([myplot(pgsol,"bus_gfm",:Ps),myplot(pgsol,"bus_gfm",:P0)])
plot([myplot(pgsol,"bus_gfm",:Qs),myplot(pgsol,"bus_gfm",:Q0)])
plot(plotv(pgsol,["bus_gfm"]))
plot(myplot(spgsol,"bus_gfm",:θ,y_norm=1/(180/pi)))
plot(myplot(pgsol,"bus_gfm",:w,y_norm=50/2/pi,y_bias=50.0))
plot(myplot(pgsol,"bus_gfm",:i_abs))
plot(myplot(pgsol,"bus_gfm",:udc))
plot(myplot(pgsol,"bus_gfm",:idc0_lim))
plot(plotv(sol[2:2:8],"bus_gfm"))


sol = Vector{PowerGridSolution}()
for i=1:4
    for j=1:-1:0
        pgsol = run_LTVS_N32_simulation(i,j,(0.0,2.0));
        push!(sol,pgsol);
    end
end

typeof.(sol)
tmp = Vector{GenericTrace}()
for i in sol
    push!(tmp,plotv(i,"bus_gfm"))
end
plot(tmp)


## VGL mit MATLAB
using MAT
dir = "\\\\fs0\\home\\liemann\\"
file = matopen(dir*"vm.mat")
vm = read(file, "vm")
close(file)
file = matopen(dir*"tmat.mat")
tm = read(file, "tout")
close(file)
file = matopen(dir*"theta.mat")
θm = read(file, "theta")
close(file)
file = matopen(dir*"freq.mat")
freq = read(file, "freq")
close(file)
file = matopen(dir*"PQ.mat")
pq = read(file, "pq")
close(file)
file = matopen(dir*"PQs.mat")
pqs = read(file, "pqs")
close(file)
file = matopen(dir*"iabs.mat")
iabs = read(file, "i_abs")
close(file)
file = matopen(dir*"udc.mat")
udc = read(file, "udc")
close(file)
file = matopen(dir*"dw.mat")
dw = read(file, "dw")
close(file)


pm = scatter(x=tm[:,1],y=vm[:,4],label="Matlab")
jm = plotv(pgsol,["bus_gfm"])
plot([pm,jm[1]])

pm = scatter(x=tm[:,1],y=θm[:,1].-18.882651941871828,label="Matlab")# 19.095781601625703
jm = myplot(pgsol,"bus_gfm",:θ,y_norm=1/180*pi)
plot([pm,jm])

pm = scatter(x=tm[:,1],y=freq[:,1],label="Matlab")
jm = myplot(pgsol,"bus_gfm",:w,y_norm=50/2/pi,y_bias=50.0)
plot([pm,jm])

pm = scatter(x=tm[:,1],y=pq[:,1],label="Matlab")
jm = myplot(pgsol,"bus_gfm",:P0)
plot([pm,jm])

pm = scatter(x=tm[:,1],y=pq[:,2],label="Matlab")
jm = myplot(pgsol,"bus_gfm",:Q0)
plot([pm,jm])

pm = scatter(x=tm[:,1],y=pqs[:,1],label="Matlab")
jm = myplot(pgsol,"bus_gfm",:Ps)
plot([pm,jm])

pm = scatter(x=tm[:,1],y=pqs[:,2],label="Matlab")
jm = myplot(pgsol,"bus_gfm",:Qs)
plot([pm,jm])

pm = scatter(x=tm[:,1],y=iabs[:,1],label="Matlab")
jm = myplot(pgsol,"bus_gfm",:i_abs)
plot([pm,jm])

pm = scatter(x=tm[:,1],y=udc[:,1],label="Matlab")
jm = myplot(pgsol,"bus_gfm",:udc,y_bias=1.0)
plot([pm,jm])

pm = scatter(x=tm[:,1],y=-dw[:,1],label="Matlab")
jm = myplot(pgsol,"bus_gfm",:w)


p4 = Vector{PlotlyJS.SyncPlot}() 
for i=2:2:8
    pm  = scatter(x=tm[:,1],y=vm[:,4+(i-1)*5],label="Matlab")
    jm = plotv(sol[i],"bus_gfm")
    push!(p4,plot([pm,jm]))
end
[p4[1] p4[2];p4[3] p4[4]]
