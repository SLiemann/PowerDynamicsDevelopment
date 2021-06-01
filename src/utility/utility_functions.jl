using ModelingToolkit
    
function getPreFaultVoltages(pg::PowerGrid,ic_prefault::Array{Float64,1},ic_endfault::Array{Float64,1})
    ind = getVoltageSymbolPositions(pg)
    ic_endfault[ind] = ic_prefault[ind]
    return ic_endfault
end

function getVoltageSymbolPositions(pg::PowerGrid)
    n = map(x->PowerDynamics.variable_index(pg.nodes,x,:u_r),collect(keys(pg.nodes)))
    append!(n,map(x->PowerDynamics.variable_index(pg.nodes,x,:u_i),collect(keys(pg.nodes))))
end

function getSymbolPosition(pg::PowerGrid,syms::Array{Symbol,1})
    return sort(indexin(syms,rhs(pg).syms))
end

function AddZerosIntoSolution(pg1::PowerGrid,pg2::PowerGrid,sol)
    key1 = collect(keys(pg1.nodes))
    key2 = collect(keys(pg2.nodes))
    n = findfirst(x->x==nothing,indexin(key1,key2))
    len_bus = length(symbolsof(pg1.nodes[key1[n]]))
    first_ind = PowerDynamics.variable_index(pg1.nodes,key1[n],:u_r)
    start_u = findfirst(x-> x== length(sol.u[end]),length.(sol.u))
    for i in start_u:length(sol.u)
        splice!(sol.u[i],first_ind:first_ind-1,zeros(len_bus))
        for j in sol.k[i]
            splice!(j,first_ind:first_ind-1,zeros(len_bus))
        end
    end
    return sol
end

function PiModel(y::Complex{Num}, y_shunt_km, y_shunt_mk, t_km, t_mk)
    Π = Array{Union{Complex{Float64},Complex{Num}}}(undef,2,2)
    Π[1, 1] = - abs2(t_km) * (y + y_shunt_km) # Our sign convention is opposite for the source of the edge
    Π[1, 2] = conj(t_km) * t_mk * y # Our sign convention is opposite for the source of the edge
    Π[2, 1] = - conj(t_mk) * t_km * y
    Π[2, 2] = abs2(t_mk) * (y + y_shunt_mk)
    Π
end

function PiModel(y, y_shunt_km, y_shunt_mk, t_km, t_mk::Num)
    Π = Array{Union{Complex{Float64},Complex{Num}}}(undef,2,2)
    Π[1, 1] = - abs2(t_km) * (y + y_shunt_km) # Our sign convention is opposite for the source of the edge
    Π[1, 2] = conj(t_km) * t_mk * y # Our sign convention is opposite for the source of the edge
    Π[2, 1] = - conj(t_mk) * t_km * y
    Π[2, 2] = abs2(t_mk) * (y + y_shunt_mk)
    Π
end

function PiModel(y, y_shunt_km, y_shunt_mk, t_km::Num, t_mk)
    Π = Array{Union{Complex{Float64},Complex{Num}}}(undef,2,2)
    Π[1, 1] = - abs2(t_km) * (y + y_shunt_km) # Our sign convention is opposite for the source of the edge
    Π[1, 2] = conj(t_km) * t_mk * y # Our sign convention is opposite for the source of the edge
    Π[2, 1] = - conj(t_mk) * t_km * y
    Π[2, 2] = abs2(t_mk) * (y + y_shunt_mk)
    Π
end
