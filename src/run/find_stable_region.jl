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
    pg, ic =  Initialize_N32_GFM(1,0);
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
        x_tmp = Xstart
        r_tmp = Rstart
        for X=Xstart:-dX:0.0
            #params = GetParamsGFM(pg)
            #prob = ODEProblem(rhs(pg),ic,(0.0,2.0),params)
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

@time vl, xr = CalcXRMap(105.0,5.0,75.0,5.0);
plot(xr[:,1],xr[:,2])

xdata = 0:5:105
ydata = 0:5:75

plot(surface(x=xdata,y=ydata,z=vl))

function plotxkrit(vl,xr)
    for i in vl[:]
end


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

rkrit,dr,vl = CalcRkrit(1im*75,1im*5.0)
