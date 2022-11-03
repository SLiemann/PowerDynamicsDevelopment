begin
    using PlotlyJS, DataFrames
    include("C:/Users/liemann/github/PowerDynamicsDevelopment/src/include_custom_nodes_lines_utilities.jl")
    include("C:/Users/liemann/github/PowerDynamicsDevelopment/src/grids/LTVS_Test_System_N32_GFM.jl")
end


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

pgsol = run