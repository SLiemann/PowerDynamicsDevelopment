using PlotlyJS, DataFrames

Sbase = 8000e6
Ubase = 400e3
Zbase = Ubase^2/Sbase

begin
    include("C:/Users/liemann/github/PowerDynamicsDevelopment/src/include_custom_nodes_lines_utilities.jl")
    include("C:/Users/liemann/github/PowerDynamicsDevelopment/src/grids/LTVS_Test_System_N32 GFM.jl")
    nothing
end

pg0,ic0 = Initialize_N32_GFM(1,1);
@time pgsol0, suc0,FRT0 = simulate_LTVS_N32_simulation(pg0,ic0,(0.0,5.0),(8.0+8.0im)/Zbase);
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
            pgsol, suc,FRT_tmp = simulate_LTVS_N32_simulation(pg,ic,(0.0,5.0),(R+1im*X)/Zbase);
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

Rverlauf = 20:-5:0.0
Xverlauf = 10:-5:0.0

@time xr, xrt = CalcXRMap(Rverlauf,Xverlauf);

plot(surface(x=Rverlauf,y=Xverlauf,z=xr))
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
save("droop_I100_ownred_Tf_300ms_R_20_1_0_X_20_1_0.jld","Rverlauf",Rverlauf,"Xverlauf",Xverlauf,"XR",xr,"XR_t",xrt)

function plotxkrit(vl_,xdata,ydata)
    x = Vector{Float64}()
    y = Vector{Float64}()

    len1 = size(vl_)[1]
    len1 = size(vl_)[2]
    
    for i=1:size(vl_)[1]-1
        ind1 = findfirst(x->x==1,vl_[i,:])
        ind2 = findfirst(x->x==1,vl_[i+1,:])

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

px1,py1 = plotxkrit(reverse(xr),collect(reverse(Rverlauf)),collect(reverse(Xverlauf)));
plot(px1,py1)

function DetermineBoundary(XR,Rv,Xv)
    r = Vector{Float64}()
    x = Vector{Float64}()
    flag_dir = "bottom"
    flag_end = false
    ind = [(0,-1);(-1,0);(0,1);(1,0)]

    Rv = reverse(Rv)
    Xv = reverse(Xv)

    lenx = size(XR)[1]
    lenr = size(XR)[2]

    #starting indices
    indr = deepcopy(lenr)
    indx = findlast(k->k==1,XR[indr,:])
    #saving first entry
    push!(r,Rv[indr])
    push!(x,Xv[indx])

    function next_ind(start_ind)
        flag_end == false
        for i=start_ind:start_ind+3
            if mod(i,4) == 0
                dr = ind[4][1]
                dx = ind[4][2]
            else
                dr = ind[mod(i,4)][1]
                dx = ind[mod(i,4)][2]
            end

            indr = indr + dr
            indx = indx + dx

            println(indr)
            println(indx)
            if indr == 0 || indr > lenr || indx == 0 || indx > lenx 
                #skip boundarys
            elseif indx == lenx 
                println("END")
                flag_end = true
                break;
            elseif XR[indr,indx] == -1
                #skip unstable case
            elseif XR[indr,indx] == 1
                push!(r,Rv[indr])
                push!(x,Xv[indx])
                if mod(i,4) == 1
                   flag_dir = "right"
                   break;
                elseif mod(i,4) == 2
                   flag_dir = "bottom"
                   break;
                elseif mod(i,4) == 3
                    flag_dir = "left"
                    break;
                elseif mod(i,4) == 0
                    flag_dir = "top"
                    break;
                else
                    error("error2")
                end
            end 
        end
        return nothing
    end    
    counter = 0;

    while flag_end == false && counter < 10
        println(flag_dir)
        if flag_dir == "top"
            next_ind(1)
        elseif flag_dir == "left"
            next_ind(2)
        elseif flag_dir == "bottom"
            next_ind(3)
        elseif flag_dir == "right"
            next_ind(4)
        else
            error("error")
        end
        counter = counter + 1
    end
    return r,x
end

r1,x1 = DetermineBoundary(xr,collect(Rverlauf),collect(Xverlauf));
plot(r1,x1)

tmp = load("droop_I070_R_18_1_0_X_08_1_0.jld")
tmp = load("droop_I00_R_10_01_0_X_55_01_0.jld")
tmp = load("matching_R_136_1_0_X_88_1_0.jld")

Rv= tmp["Rverlauf"];
Xv = tmp["Xverlauf"];
XRv= tmp["XR"];
XRt= tmp["XR_t"];

plot(surface(x=Rv,y=Xv,z=XRv))
px,py = plotxkrit(reverse(XRv),collect(reverse(Rv)),collect(reverse(Xv)));
plot([scatter(x=px1,y=py1),scatter(x=px,y=py)])

using MAT
dir = "C:\\Users\\liemann\\Desktop\\"
file = matopen(dir*"matlab.mat")
XRm = read(file, "res")
close(file)

p4 = scatter(x=px,y=py,name="PD")
p2= scatter(x=XRm[:,1],y=XRm[:,2],name="MATLAB")
plot([p1,p2,p3,p4])