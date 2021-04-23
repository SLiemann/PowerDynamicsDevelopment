function getPreFaultVoltages(pg::PowerGrid,ic_prefault::Array{Float64,1},ic_endfault::Array{Float64,1})
    ind = getVoltageSymbolPositions(pg)
    ic_endfault[ind] = ic_prefault[ind]
    return ic_endfault
end

function getVoltageSymbolPositions(pg::PowerGrid)
    u_syms = [[Symbol("u_r_",i) for i in 1:length(pg.nodes)]; [Symbol("u_i_",i) for i in 1:length(pg.nodes)]]
    return getSymbolPosition(pg,u_syms)
end

function getSymbolPosition(pg::PowerGrid,syms::Array{Symbol,1})
    return sort(indexin(syms,rhs(pg).syms))
end
