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
@time pgsol0, evr_sol = simulate_LTVS_N32_simulation_TS(pg0,ic0,(0.0,0.6),(20.0+1im*20)/Zbase);
plotallvoltages(pgsol0);
plot(myplot(pgsol0,"bus_gfm",:q_imax))
plot(myplot(pgsol0,"bus_gfm",:q_idcmax))
plot(myplot(pgsol0,"bus_gfm",:idc0))
plot(myplot(pgsol0,"bus_gfm",:udc))


###### Trajectory Sensitivity Analysis ######
include("C:/Users/liemann/github/PowerDynamicsDevelopment/src/sensitivity_analyses/Local_Sensitivity.jl")
mtk1 = modelingtoolkitize(pgsol0.dqsol.prob)
mtk2 = modelingtoolkitize(ODEProblem(rhs(GetPostFaultLTVSPG_TS(pg0)), ic0, (0.0, 1.0), pgsol0.dqsol.prob.p))
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

# 17 = θ, 20=udc
# 8 = Kp_droop, 9 = Kp_uset
pg_labels = string.(rhs(pg0).syms)
u_sensis = [17,20];
p_sensis = [8,9,10,18,19,20,21,22,23];

@time hybrid_sen,Δτ = CalcHybridTrajectorySensitivity(mtk,pgsol0.dqsol,evr_sol,s,hs,u_sensis,p_sensis);

size(hybrid_sen[1])

plot(pgsol0.dqsol.t,hybrid_sen[2][25,:])


### Saving Sensitivities
using MATLAB

labels_p = [
    "rfault", #1
    "xfault", #2
    "line_1", #3
    "line_2", #4
    "tap_pos", #5
    "p0set", #6
    "u0set", #7
    "Kp_droop", #8
    "Kp_uset", #9
    "Ki_uset", #10
    "Kdc", #11
    "gdc", #12
    "cdc", #13
    "xlf", #14
    "rf", #15
    "xcf", #16
    "Tdc", #17
    "Kp_u", #18
    "Ki_u", #19
    "Kp_i", #20
    "Ki_i", #21
    "imax_csa", #22
    "imax_dc", #23
    "LVRT_on", #24
    ];

sensi_labels = [pg_labels[u_sensis];labels_p[p_sensis]];
state_labels = pg_labels;
odesol = pgsol0.dqsol[:,:];
path = "C:\\Users\\liemann\\github\\PowerDynamicsDevelopment\\src\\results\\"

write_matfile(path*"sensis_GFM_SECM_BC.mat"; odesol = odesol,sensis = hybrid_sen, delta_tau = Δτ, sensi_labels=sensi_labels,state_labels =state_labels,evr=evr_sol) 