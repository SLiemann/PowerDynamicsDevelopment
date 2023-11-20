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
@time pgsol0, evr_sol = simulate_LTVS_N32_simulation_TS(pg0,ic0,(0.0,0.02),(20.0+1im*20)/Zbase);
plotallvoltages(pgsol0);
plot(myplot(pgsol0,"bus_gfm",:q_imax))
plot(myplot(pgsol0,"bus_gfm",:q_idcmax))
plot(myplot(pgsol0,"bus_gfm",:idc0))
plot(myplot(pgsol0,"bus_gfm",:udc))


###### Trajectory Sensitivity Analysis ######
include("C:/Users/liemann/github/PowerDynamicsDevelopment/src/sensitivity_analyses/Local_Sensitivity.jl")
mtk1 = modelingtoolkitize(pgsol0.dqsol.prob)
mtk2 = modelingtoolkitize(ODEProblem(rhs(GetPostFaultLTVSPG(pg0)), ic0, (0.0, 1.0), pgsol0.dqsol.prob.p))
mtk = [mtk1,mtk2];
sym_states = states(mtk[1])
symp = parameters(mtk[1])
eq = equations(mtk[1]);

# iset_abs = 7, idc0 = 19, q_imax = 32, q_idcmax = 33,
# imax_csa = 22, idcmax = 23

# Bei h(x,q,p) darauf achten, dass keine algebraischen Zustände resettet werden

d,a  = GetFactorisedSymbolicStates(mtk[1]);
hs = Vector{Vector{Num}}(undef,0)
ident = Num.(vcat(d,symp))
ident[32] = 1.0
push!(hs,ident)
ident = Num.(vcat(d,symp))
ident[32] = 0.0
push!(hs,ident)
ident = Num.(vcat(d,symp))
ident[33] = 1.0
push!(hs,ident)
ident = Num.(vcat(d,symp))
ident[33] = 0.0
push!(hs,ident)

s = Vector{Num}(undef,0)
push!(s,(sym_states[7] - symp[22] - 10*sym_states[32]))
push!(s,((abs(sym_states[7]) - symp[22])*sym_states[32]))
push!(s,(sym_states[19] - symp[23] - 10*sym_states[33]))
push!(s,((abs(sym_states[19]) - symp[23])*sym_states[33]))

# 17 = θ, 20=0.005
# 8 = Kp_droop, 9 = Kp_uset
@time tmp =  InitTrajectorySensitivity(mtk,pgsol0.dqsol,[17,20],[8,9]);

@time hybrid_sen,Δτ = CalcHybridTrajectorySensitivity(mtk,pgsol0.dqsol,evr_sol,s,hs,[17,20],[8,9]);

size(hybrid_sen[1])
plot(pgsol0.dqsol.t,hybrid_sen[4][17,:])