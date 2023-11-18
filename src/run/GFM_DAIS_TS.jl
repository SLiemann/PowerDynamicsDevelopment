using PlotlyJS, DataFrames
using ModelingToolkit   

Sbase = 8000e6
Ubase = 400e3
Zbase = Ubase^2/Sbase

begin
    include("C:/Users/liemann/github/PowerDynamicsDevelopment/src/include_custom_nodes_lines_utilities.jl")
    include("C:/Users/liemann/github/PowerDynamicsDevelopment/src/grids/LTVS_Test_System_N32 GFM_TS.jl")
    nothing
end

pg0,ic0 = Initialize_N32_GFM_TS();
@time pgsol0, evr_sol = simulate_LTVS_N32_simulation_TS(pg0,ic0,(0.0,0.21),(1200.0)/Zbase);
plotallvoltages(pgsol0)
plot([myplot(pgsol0,"bus_gfm",:LVRT),plotv(pgsol0,["bus_gfm"])[1]])

plot(myplot(pgsol0,"bus_gfm",:q_imax))
plot(myplot(pgsol0,"bus_gfm",:udc))


###### Trajectory Sensitivity Analysis ######
include("C:/Users/liemann/github/PowerDynamicsDevelopment/src/sensitivity_analyses/Local_Sensitivity.jl")
mtk = modelingtoolkitize(pgsol0.dqsol.prob)
sym_states = states(mtk)
symp = parameters(mtk)
eq = equations(mtk)

#rhs(pg0).syms .=> sym_states
#p_float = GetParamsGFM_TS(pg0)[23]
# iset_abs = 7, idc0 = 19, q_imax = 32, q_idcmax = 33,
# imax_csa = 22, idcmax = 23

hs = Vector{Vector{Num}}(undef,0)
ident = Num.(vcat(sym_states,symp))
ident[32] = 1.0
push!(hs,ident)
ident = Num.(vcat(sym_states,symp))
ident[32] = 0.0
push!(hs,ident)
ident = Num.(vcat(sym_states,symp))
ident[33] = 1.0
push!(hs,ident)
ident = Num.(vcat(sym_states,symp))
ident[33] = 0.0
push!(hs,ident)

s = Vector{Num}(undef,0)
push!(s,(sym_states[7] - symp[22] - 10*sym_states[32]))
push!(s,((abs(sym_states[7]) - symp[22])*sym_states[32]))
push!(s,(sym_states[19] - symp[23] - 10*sym_states[33]))
push!(s,((abs(sym_states[19]) - symp[23])*sym_states[33]))

@time tmp =  InitTrajectorySensitivity([mtk],pgsol0.dqsol);

@time hybrid_sen,Δτ = CalcHybridTrajectorySensitivity([mtk],pgsol0.dqsol,evr_sol,s,hs);