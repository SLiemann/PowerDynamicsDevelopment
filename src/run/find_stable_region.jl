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
pg_post0 = GetPostFaultLTVSPG(pg0);


params0 = GetParamsGFM(pg0)
prob0 = ODEProblem(rhs(pg0),ic0,(0.0,2.0),params0)
@time pgsol0, suc0 = simulate_LTVS_N32_simulation(pg0,ic0,(0.0,5.0),(60+1im*75.0)/Zbase);
plot(plotv(pgsol0,["bus_gfm"]))


function CalcXRMap(Rstart,dRstart,Xstart,dXstart;nenner=2.0,limit=1e-2)
    pg, ic =  Initialize_N32_GFM(3,0);
    R = Rstart
    X = Xstart
    dR = dRstart
    dX = dXstart

    length_dr = length(Rstart:-dR:0.0)
    length_dx = length(Xstart:-dX:0.0)
    Verlauf = zeros(length_dr,length_dx)
    cache_xrkrit = Matrix{Float64}(undef,1,2)

    ind_r = length_dr
    ind_x = length_dx

    for R=Rstart:-dR:0.0
        ind_x  = length_dx
        suc_tmp = :Success

        for X=Xstart:-dX:0.0

            pgsol, suc = simulate_LTVS_N32_simulation(pg,ic,(0.0,2.0),(R+1im*X)/Zbase);
            if suc == :DtLessThanMin
                println(R," ",X)
                Verlauf[ind_r,ind_x] = -1
                break;
            elseif suc == :Unstable
                Verlauf[ind_r,ind_x] = 0
                if suc_tmp == :Success
                    cache_xrkrit = vcat(cache_xrkrit,[R X+dX])
                end
                break
            elseif suc == :Success
                Verlauf[ind_r,ind_x] = 1
                if X == 0.0
                    cache_xrkrit = vcat(cache_xrkrit,[R X])
                end
            else
                Verlauf[ind_r,ind_x] = -10
            end
            ind_x  = ind_x -1
            suc_tmp = deepcopy(suc)
        end
        ind_r  = ind_r -1
    end
    return Verlauf, cache_xrkrit[2:end,:]
end

@time vl, xr = CalcXRMap(100.0,5.0,80.0,5.0);

#droop = R=103, x=66
#mathing R=136, X = 88
#dVOC R=100, X = 65
plot(xr[:,1],xr[:,2])

xdata = 0:5:100
ydata = 0:5:80

plot(surface(x=xdata,y=ydata,z=vl))

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

