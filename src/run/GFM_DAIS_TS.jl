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
@time pgsol2, evr_sol = simulate_LTVS_N32_simulation_TS(pg0,ic0,(0.0,200.0),1e6*(20.0+1im*20)/Zbase);
p1 = plotallvoltages(pgsol1);
myplot(pgsol2,"bus_gfm",:i_abs)
plot(myplot(pgsol0,"bus_gfm",:q_imax))
plot(myplot(pgsol0,"bus_gfm",:q_idcmax))
plot(myplot(pgsol0,"bus_gfm",:idc0))
plot([myplot(pgsol0,"bus_gfm",:θ),myplot(pgsol1,"bus_gfm",:θ)])

plot([myplot(pgsol0,"bus_gfm",:P0),myplot(pgsol0_per,"bus_gfm",:P0)])


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
x, dp = extract_local_sensitivities(sensi_ad);
pg_labels = string.(rhs(pg0).syms);
state_labels = pg_labels;
path = "C:\\Users\\liemann\\github\\PowerDynamicsDevelopment\\src\\results\\"

write_matfile(path*"STVS_AD_sensis_droop_5300_SECM_BC.mat"; odesol = x[:,:],sensis = dp, sensi_labels=labels_p,state_labels =pg_labels,evr=evr_sol, time=sensi_ad.t) 

write_matfile(path*"Approximated_Kd_imax_droop_5142_SECM_BC.mat"; odesol_kd = pgsol1.dqsol[:,:], odesol_imax = pgsol2.dqsol[:,:],state_labels =pg_labels, time_kd=pgsol1.dqsol.t,time_imax = pgsol2.dqsol.t) 

#####  Approximated Solution
using MAT
sensi_Data = matread("\\\\fs0\\home\\liemann\\Diss\\Results\\621_Short_Term_GFM\\sensis_GFM_SECM_BC.mat")
sensis = sensi_Data["sensis"]
sol = sensi_Data["odesol"];
Δkp = 0.5*pi
sol_appr = ApproximatedTrajectory(sensi_Data["odesol"],sensis[3],Δkp)
ode_apprx = deepcopy(pgsol0.dqsol);
ode_apprx[:,:] = sol_appr;
pgsol_apprx = PowerGridSolution(ode_apprx,pg0);

####### Solution from Automatic Differentiation
x, dp = extract_local_sensitivities(sensi_ad)
da = dp[8]

PlotlyJS.plot(sensi_ad.t,dp[22]')

sol_appr_ad = ApproximatedTrajectory(x[:,:],da,Δkp);

x5r = x[15,:] .+ da[15,:]*Δkp
x5i = x[16,:] .+ da[16,:]*Δkp
x5 = hypot.(x5r,x5i)
px5 = PlotlyJS.scatter(x=pgsol0_ad.t,y=x5);
PlotlyJS.plot([plotv(pgsol0,"bus_gfm"),px5,plotv(pgsol_apprx,"bus_gfm"),plotv(pgsol0_per,"bus_gfm")])

zx = 17
zs = :P0
xz = x[zx,:] .+ da[zx,:]*Δkp
px5 = PlotlyJS.scatter(x=pgsol0_ad.t,y=xz);
PlotlyJS.plot([myplot(pgsol0,"bus_gfm",zs),px5,myplot(pgsol_apprx,"bus_gfm",zs),myplot(pgsol0_per,"bus_gfm",zs)])

PlotlyJS.plot(PlotlyJS.scatter(x=pgsol0_ad.t,y=da[zx,:]))


x, dp = extract_local_sensitivities(pgsol0);
da = dp[8]