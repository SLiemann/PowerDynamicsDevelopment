using PlotlyJS, Plots, DataFrames
using ModelingToolkit   

Sbase = 8000e6
Ubase = 400e3
Zbase = Ubase^2/Sbase

begin
    include("C:/Users/liemann/github/PowerDynamicsDevelopment/src/include_custom_nodes_lines_utilities.jl")
    include("C:/Users/liemann/github/PowerDynamicsDevelopment/src/grids/LTVS_Test_System_N32 PEL_TS.jl")
    nothing
end
# unstable case: share_pe = 0.3; Rf = 10, Xf = 0
pg0,ic0 = Initialize_N32_PEL_TS(share_pe= 0.30);
@time pgsol0, evr_sol = simulate_LTVS_N32_simulation_PEL_TS(pg0,ic0,(0.0,1.0),(20.0+1im*20)/Zbase);
plotallvoltages(pgsol0);

myplot(pgsol0,"bus_load",:q1)
myplot(pgsol0,"bus_load",[:p1,:ps]);
myplot(pgsol0,"bus_load",:tsum);
myplot(pgsol0,"bus_load",:ton);
myplot(pgsol0,"bus_load",:toff);
myplot(pgsol0,"bus_load",:vofft2);
myplot(pgsol0,"bus_load",:Vabstoff);
myplot(pgsol0,"bus_load",:q_on);

###### Trajectory Sensitivity Analysis ######
include("C:/Users/liemann/github/PowerDynamicsDevelopment/src/sensitivity_analyses/Local_Sensitivity.jl")
mtk1 = modelingtoolkitize(pgsol0.dqsol.prob);
mtk2 = modelingtoolkitize(ODEProblem(rhs(GetPostFaultLTVSPG_TS(pg0)), ic0, (0.0, 1.0), pgsol0.dqsol.prob.p))
mtk = [mtk1,mtk2];
sym_states = states(mtk[1]);
symp = parameters(mtk[1]);
eq = equations(mtk[1]);

d,a  = GetFactorisedSymbolicStates(mtk[1]);
hs = Vector{Vector{Num}}(undef,0)
ident = Num.(vcat(d,symp))
push!(hs,ident) #dummy


s = Vector{Num}(undef,0)
push!(s,(sym_states[1]))

# 17 = θ, 20=udc
# 8 = Kp_droop, 9 = Kp_uset
pgsysm = string.(rhs(pg0).syms)
display.(collect(eachindex(pgsysm)) .=> pgsysm)
u_sensis = [18,19];
p_sensis = [6,7];

@time hybrid_sen,Δτ = CalcHybridTrajectorySensitivity(mtk,pgsol0.dqsol,evr_sol,s,hs,u_sensis,p_sensis);

PlotlyJS.plot(pgsol0.dqsol.t,hybrid_sen[3][18,:])


### Saving Sensitivities
using MATLAB

labels_p = [
    "rfault", #1
    "xfault", #2
    "line_1", #3
    "line_2", #4
    "tap_pos", #5
    "Cd", #6
    "Pdc", #7
    ];

sensi_labels = [pg_labels[u_sensis];labels_p[p_sensis]];
state_labels = pg_labels;
odesol = pgsol0.dqsol[:,:];
path = "C:\\Users\\liemann\\github\\PowerDynamicsDevelopment\\src\\results\\"

write_matfile(path*"sensis_PEL_nPFC_BC.mat"; odesol = odesol,sensis = hybrid_sen, delta_tau = Δτ, sensi_labels=sensi_labels,state_labels =state_labels,evr=evr_sol) 