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
@time pgsol0, suc0 = simulate_LTVS_N32_simulation(pg0,ic0,(0.0,5.0),(76+1im*56.0)/Zbase);
plot(plotv(pgsol0,["bus_gfm"]))


function CalcXRMap(Rstart,dRstart,Xstart,dXstart;nenner=2.0,limit=1e-2)
    pg, ic =  Initialize_N32_GFM(1,0);
    pg_post = GetPostFaultLTVSPG(pg);
    R = Rstart
    X = Xstart
    Zkrit = R + 1im*X
    dR = dRstart
    dX = dXstart

    length_dr = length(Rstart:-dR:0.0)
    length_dx = length(Xstart:-dX:0.0)
    Verlauf = zeros(length_dr,length_dx)

    ind_r = length_dr
    ind_x = length_dx

    for R=Rstart:-dR:0.0
        ind_x  = length_dx
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
                break
            elseif suc == :Success
                Verlauf[ind_r,ind_x] = 1
            else
                Verlauf[ind_r,ind_x] = -10
            end
            ind_x  = ind_x -1
        end
        ind_r  = ind_r -1
    end
    return Verlauf
end

@time vl = CalcXRMap(103.0,1.0,66.0,1.0);

xdata = 0:1:103
ydata = 0:1:66

plot(surface(x=xdata,y=ydata,z=vl))

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
