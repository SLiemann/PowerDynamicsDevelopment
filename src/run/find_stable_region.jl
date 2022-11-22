using PlotlyJS, DataFrames

Sbase = 8000e6
Ubase = 400e3
Zbase = Ubase^2/Sbase

begin
    include("C:/Users/liemann/github/PowerDynamicsDevelopment/src/include_custom_nodes_lines_utilities.jl")
    include("C:/Users/liemann/github/PowerDynamicsDevelopment/src/grids/LTVS_Test_System_N32 GFM.jl")
    nothing
end

pg0, ic0 =  Initialize_N32_GFM(2,0);

@time pgsol0, suc0,FRT0 = simulate_LTVS_N32_simulation(pg0,ic0,(0.0,4.0),(134+0im)/Zbase);
plot([myplot(pgsol0,"bus_gfm",:LVRT),plotv(pgsol0,["bus_gfm"])[1]])
plot(myplot(pgsol0,"bus_gfm",:iset_abs))
plot(myplot(pgsol0,"bus_gfm",:idc0))
plot(myplot(pgsol0,"bus_gfm",:udc))


function CalcXRMap(Rrange, Xrange)
    pg, ic =  Initialize_N32_GFM(2,0);

    length_dr = length(Rrange)
    length_dx = length(Xrange)
    XR = zeros(length_dr,length_dx)
    XR_tend = similar(XR)
    cache_xrkrit = Matrix{Float64}(undef,1,2)

    for (indR,R) = enumerate(Rrange)
        println(R)
        for (indX,X) = enumerate(Xrange)
            pgsol, suc,FRT_tmp = simulate_LTVS_N32_simulation(pg,ic,(0.0,4.0),(R+1im*X)/Zbase);
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

Rverlauf = 90:-1:45
Xverlauf = 98:-1:95.0

@time xr, xrt = CalcXRMap(Rverlauf,Xverlauf);

#droop = R=103, x=66
#mathing R=135, X = 87 ODER R=150, X = 98
#dVOC R=100, X = 65 (not finished)
plot(xr[:,1],xr[:,2])
plot(Rverlauf,xr)
plot(surface(x=Rverlauf,y=Xverlauf,z=xr))

droop_save = ["Rverlauf" => Rverlauf, "Xverlauf"=>Xverlauf,"XR"=>xr]

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
        if i == size(vl)[2]-1
            push!(x,xdata[i+1])
            push!(y,ydata[ind2])
        end
    end
    return x,y
end
px,py = plotxkrit(vl,collect(xdata),collect(ydata));
plot(px,py)

vl_match = deepcopy(vl)
px_match = deepcopy(px)
py_match = deepcopy(py)

function CalcRkrit(Rstart,dRstart;nenner=2.0,limit=1e-2)
    pg, ic =  Initialize_N32_GFM(1,0);
    R = Rstart
    Rkrit = R
    dR = dRstart
    Verlauf = Matrix{ComplexF64}(undef,1,3)
    while abs(dR) > limit
        params = GetParamsGFM(pg)
        prob = ODEProblem(rhs(pg),ic,(0.0,2.0),params)
        pgsol, suc = simulate_LTVS_N32_simulation(prob,(R)/Zbase);
        Verlauf = vcat(Verlauf,[imag(R) real(R) suc])
        if suc == 0
            R = R + dR/nenner
            dR = dR/nenner
        else
            Rkrit = R
            R = R - dR 
        end
    end
    return Rkrit,dR*nenner,Verlauf
end



function CalcXRMapNEW(Rstart,Xstart)
    pg, ic =  Initialize_N32_GFM(1,0);

    cache_xrkrit = Matrix{Float64}(undef,1,2)

    R = Rstart
    Rtmp = Inf
    X = Xstart
    Xtmp = 0
    #while abs(R-Rtmp) > 1
    for R=Rstart:-1.0:0.0
        println(R)
        while abs(X-Xtmp) > 1
            pgsol, suc = simulate_LTVS_N32_simulation(pg,ic,(0.0,2.0),(R+1im*X)/Zbase);
            if Xtmp == 0
                dX = X/5
            else
                dX = abs(X-Xtmp)/5
            end
            Xtmp = deepcopy(X)
            if suc == :Success
                X = X - dX
            else
                X = X + dX
            end         
        end
        cache_xrkrit = vcat(cache_xrkrit,[R X])
        X= Xstart
        Xtmp = 0
    end
    return cache_xrkrit[2:end,:]
end
tmp = 4
@time tmp = CalcXRMapNEW(103,66);
p1 = scatter(x=tmp[:,1],y=tmp[:,2])
p2 = scatter(x=px_droop,y=py_droop)
plot([p1,p2])

