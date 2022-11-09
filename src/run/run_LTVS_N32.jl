using PowerDynamics
#using Plots
using IfElse

begin
    include("C:/Users/liemann/github/PowerDynamicsDevelopment/src/include_custom_nodes_lines_utilities.jl")
    include("C:/Users/liemann/github/PowerDynamicsDevelopment/src/grids/LTVS_Test_System_N32.jl")
    #include("C:/Users/liemann/github/PowerDynamicsDevelopment/src/sensitivity_analyses/Local_Sensitivity.jl")
end
begin
    pg = LTVS_Test_System_N32()
    Qmax   = [Inf,Inf, Inf, Inf,Inf,5300/8000*sqrt(1-0.8378^2),Inf]
    Qmin   = -Qmax
    U1,δ1,ic0,cu = PowerFlowClassic(pg,iwamoto = false,max_tol = 1e-4,iter_max = 100,Qmax = Qmax, Qmin = Qmin,Qlimit_iter_check=80)
    #display(U1.=> [1.0;U])
    #display(U.=> δ)
    pg, ic0 = InitializeInternalDynamics(pg,ic0)
    #display(rhs(pg).syms .=> ic0)
    pgsol  = run_LTVS_N32_simulation(pg,ic0,(0.0,165.0));
    nothing
end
plot(cu')
println.(rhs(pg).syms .=> ic0)

CalcEigenValues(pg,[],plot =true, output = true)
xlims!((-3,0))
plot(pgsol,["bus_ehv","bus_hv","bus_load","bus_sm"],:v, legend = legend=:bottomright)
plot(pgsol,"bus_sm",:timer)
plot(pgsol,"bus_sm",:ifd)
plot(pgsol,"bus_sm",:e_qs)
plot(pgsol,"bus_sm",:e_ds)
plot(pgsol,"bus_sm",:e_qss)
plot(pgsol,"bus_sm",:e_dss)
plot!(pgsol,"bus_load",:v,ylims=(0.98,1),xlims=(0,150))


display(U1.=> δ1)

Uc = U.*exp.(1im*δ/180*pi)
Ykk = NodalAdmittanceMatrice(pg)
S  = round.(Uc.*(conj.(Ykk)*conj.(Uc)),digits=5)*8000e6
abs.(S)


using CSV, DataFrames

upf = CSV.read("\\\\fs0\\home\\liemann\\u.csv", DataFrame)
plot(upf[:,1])
plot(upf[:,1],upf[:,5],label="PF",linestyle=:dash)
plot(pgsol,["bus1","bus_ehv","bus_hv","bus_load","bus_sm"],:v,label="julia",legend=:bottomleft)


ur = ExtractResult(pgsol,"bus_ehv",:u_r)
ui = ExtractResult(pgsol,"bus_ehv",:u_i)
u = sqrt.(ur.^2+ui.^2)
plot(u)




using PlotlyJS, DataFrames


function Sol2DF(pgsol::PowerGridSolution)
    sol = DataFrame(pgsol.dqsol)

    for (ind,val) in enumerate(pgsol.powergrid.nodes)
        ur = ExtractResult(pgsol,Symbol("u_r_"*string(ind)))
        ui = ExtractResult(pgsol,Symbol("u_i_"*string(ind)))
        u = sqrt.(ur.^2 + ui.^2)
        sol[!,Symbol("uabs_"*string(ind))] = u;
    end
    return sol
end

newdf = Sol2DF(pgsol)



function plotallvoltages(pgsol::PowerGridSolution)
    newdf = Sol2DF(pgsol)
    p = Vector{GenericTrace}()
    for (ind,val) in enumerate(pgsol.powergrid.nodes)
        tmp = scatter(x=newdf.timestamp,y=newdf[:,"uabs_"*string(ind)],name=val[1])
        push!(p,tmp)
    end
    display(plot(p))
    return p
end

using MAT
file = matopen("C:\\Users\\liemann\\Desktop\\tmat.mat")
tm = DataFrame(read(file, "tout2")[1:10:end,:],:auto)
close(file)

file = matopen("C:\\Users\\liemann\\Desktop\\vmat.mat")
vm = DataFrame(read(file, "Vsimemt_pe")[1:10:end,:],:auto)
close(file)

for i = 1:size(vm)[2]
    push!(myp,scatter(x=tm[:,1],y=vm[:,i],line=attr(dash="dash")))
end
plot(myp)
