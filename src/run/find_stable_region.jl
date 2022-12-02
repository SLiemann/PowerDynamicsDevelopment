using PlotlyJS, DataFrames

Sbase = 8000e6
Ubase = 400e3
Zbase = Ubase^2/Sbase

begin
    include("C:/Users/liemann/github/PowerDynamicsDevelopment/src/include_custom_nodes_lines_utilities.jl")
    include("C:/Users/liemann/github/PowerDynamicsDevelopment/src/grids/LTVS_Test_System_N32 GFM.jl")
    nothing
end

pg0, ic0 =  Initialize_N32_GFM(1,0);

@time pgsol0, suc0,FRT0 = simulate_LTVS_N32_simulation(pg0,ic0,(0.0,30.0),(4+3.0im)/Zbase);
plot([myplot(pgsol0,"bus_gfm",:LVRT),plotv(pgsol0,["bus_gfm"])[1]])
plot(plotv(pgsol0,["bus_gfm"])[1])
plot(myplot(pgsol0,"bus_gfm",:iset_abs))
plot(myplot(pgsol0,"bus_gfm",:idc0))
plot(myplot(pgsol0,"bus_gfm",:udc))
plot(myplot(pgsol0,"bus_gfm",:θ))
plot(myplot(pgsol0,"bus_gfm",:w))

plot(myplot(pgsol0,"bus_gfm",:Pf))
plot(myplot(pgsol0,"bus_gfm",:θ,y_norm=pi/180,y_bias=18.5807))

using MAT
dir = "\\\\fs0\\home\\liemann\\"
dir = "C:\\Users\\liemann\\Desktop\\"
file = matopen(dir*"matlab.mat")
XRm = read(file, "vmt");
close(file)

p1 = plotv(pgsol0,["bus_gfm"])[1]
p2 = scatter(x=XRm[:,1],y=XRm[:,2],name="MATLAB")
plot([p1,p2])

p1 = myplot(pgsol0,"bus_gfm",:θ,y_norm=pi/180,y_bias=18.5807)
p2 = scatter(x=XRm[:,1],y=XRm[:,4],name="MATLAB")
plot([p1,p2])

p1 = myplot(pgsol0,"bus_gfm",:P0)
p2 = scatter(x=XRm[:,1],y=XRm[:,3],name="MATLAB")
plot([p1,p2])

p1 = myplot(pgsol0,"bus_gfm",:w)
p2 = scatter(x=XRm[:,1],y=XRm[:,5],name="MATLAB")
plot([p1,p2])

p1 = myplot(pgsol0,"bus_gfm",:Q0)
p2 = scatter(x=XRm[:,1],y=XRm[:,6],name="MATLAB")
plot([p1,p2])




function CalcXRMap(Rrange, Xrange)
    pg, ic =  Initialize_N32_GFM(1,0);

    length_dr = length(Rrange)
    length_dx = length(Xrange)
    XR = zeros(length_dr,length_dx)
    XR_tend = similar(XR)
    cache_xrkrit = Matrix{Float64}(undef,1,2)

    for (indR,R) = enumerate(Rrange)
        println(R)
        for (indX,X) = enumerate(Xrange)
            pgsol, suc,FRT_tmp = simulate_LTVS_N32_simulation(pg,ic,(0.0,30.0),(R+1im*X)/Zbase);
            if  suc == :DtLessThanMin 
                XR[indR,indX] = -3;
            elseif suc == :Unstable
                XR[indR,indX] = -2;
            elseif suc != :Success && FRT_tmp == 1 
                XR[indR,indX] = -4;
            else
                XR[indR,indX] = FRT_tmp;
            end
            XR_tend[indR,indX] = pgsol.dqsol.t[end]
        end
    end
    return XR, XR_tend
end

Rverlauf = 12:-1:0.0
Xverlauf = 8:-1:0.0

@time xr, xrt = CalcXRMap(Rverlauf,Xverlauf);

#droop = R=103, x=66
#mathing R=135, X = 87 ODER R=150, X = 98
#dVOC R=100, X = 65 (not finished)
plot(xr[:,1],xr[:,2])
plot(Rverlauf,xr)
plot(surface(x=Rverlauf,y=Xverlauf,z=xrt))
plot(surface(x=Rverlauf,y=Xverlauf,z=winkel))
plot(surface(x=Rverlauf,y=Xverlauf,z=betrag))

winkel = zeros(length(Rverlauf),length(Xverlauf))
betrag = similar(winkel)

for (indR,R) in enumerate(Rverlauf)
    for (indX,X) in enumerate(Xverlauf)
        winkel[indR,indX] = atan(X,R)*180/pi
        betrag[indR,indX] = hypot(R,X)
    end
end


using FileIO
save("droop_I00_R_12_1_0_X_08_1_0.jld","Rverlauf",Rverlauf,"Xverlauf",Xverlauf,"XR",xr,"XR_t","xrt")

function plotxkrit(vl_,xdata,ydata)
    x = Vector{Float64}()
    y = Vector{Float64}()

    for i=1:size(vl_)[1]-1
        ind1 = findfirst(x->x==1,vl_[i,:])
        ind2 = findfirst(x->x==1,vl_[i+1,:])
        #print(ind1)
        #print(ind2)
        if i==1 
            push!(x,xdata[i])
            push!(y,ydata[ind1])
        end
        if ind2 > ind1
            push!(x,xdata[i])
            push!(y,ydata[ind2])
            push!(x,xdata[i+1])
            push!(y,ydata[ind2])
        end
        if ind2 < ind1
            push!(x,xdata[i+1])
            push!(y,ydata[ind1])
            push!(x,xdata[i+1])
            push!(y,ydata[ind2])
        end
        if ind2 == ind1
            push!(x,xdata[i+1])
            push!(y,ydata[ind1])
        end
        if i == size(vl_)[2]-1
            push!(x,xdata[i+1])
            push!(y,ydata[ind2])
        end
    end
    return x,y
end
px,py = plotxkrit(reverse(xr),collect(reverse(Rverlauf)),collect(reverse(Xverlauf)));
plot(px,py)


tmp = load("droop_R_110_1_0_X_70_1_0.jld")
tmp = load("droop_I090_R_64_1_0_X_40_1_0.jld")
tmp = load("matching_R_136_1_0_X_88_1_0.jld")

Rv= tmp["Rverlauf"];
Xv = tmp["Xverlauf"];
XRv= tmp["XR"];
plot(surface(x=Rv,y=Xv,z=XRv))
px,py = plotxkrit(reverse(XRv),collect(reverse(Rv)),collect(reverse(Xv)));


using MAT
dir = "C:\\Users\\liemann\\Desktop\\"
file = matopen(dir*"matlab2.mat")
XRm = read(file, "res")
close(file)

p1 = scatter(x=px,y=py,name="PD")
p2= scatter(x=XRm[:,1],y=XRm[:,2],name="MATLAB")
plot([p1,p2])