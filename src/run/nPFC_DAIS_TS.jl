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
pg0,ic0 = Initialize_N32_PEL_TS(share_pe= 0.3);
@time sensi_ad, evr_sol = simulate_LTVS_N32_simulation_PEL_TS(pg0,ic0,(0.0,150.0),(20.0+1im*20)/Zbase);
plotallvoltages(sensi_ad);

myplot(pgsol0,"bus_load",:q1)
myplot(sensi_ad,"bus_sm",:ifd)
myplot(sensi_ad,"bus_sm",:timer)

myplot(pgsol0,"bus_load",[:p1,:ps]);
myplot(pgsol0,"bus_load",[:q1,:qs]);

myplot(pgsol0,"bus_load",:ton);
myplot(pgsol0,"bus_load",:toff);
myplot(pgsol0,"bus_load",:vofft2);
myplot(pgsol0,"bus_load",:voff);
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

x, dp = extract_local_sensitivities(sensi_ad)
v3 = sqrt.(x[7,:].^2 + x[8,:].^2);

v3_appr = sqrt.((x[7,:].+dp[7][7,:].*0.2).^2 + (x[8,:].+dp[7][8,:].*0.2).^2);

Plots.plot(sensi_ad.t,v3)
Plots.plot!(sensi_ad.t,v3_appr)


ps = x[19,:]
ps_appr = ps + dp[7][19,:].*0.1
Plots.plot(-ps)
Plots.plot!(-ps_appr)


### Saving Sensitivities
using MATLAB
x, dp = extract_local_sensitivities(sensi_ad)
labels_p = [
    "rfault", #1
    "xfault", #2
    "line_1", #3
    "line_2", #4
    "tap_pos", #5
    "Cd", #6
    "Pdc", #7
    ];
state_labels = string.(rhs(pg0).syms)
path2 = "C:\\Users\\liemann\\github\\PowerDynamicsDevelopment\\src\\results\\"
write_matfile(path2*"AD_sensis_GENTPJ_PEL_share_0.3_BC_1em4_1s_NOqoffset.mat"; odesol = x[:,:],sensis = dp, sensi_labels=labels_p,state_labels =state_labels,evr=evr_sol,time =sensi_ad.t) 

write_matfile(path2*"LTVS_odesol_GENTPJ_PEL_share_0.4_noBC_withOffset_2em4.mat"; odesol = pgsol0.dqsol[:,:],state_labels =state_labels,time =pgsol0.dqsol.t) 

write_matfile(path2*"LTVS_odesol_GENTPJ_PEL_share_0.3_BC_NoOffset_1em3.mat"; odesol = sensi_ad.dqsol[:,:],state_labels =state_labels,time =sensi_ad.dqsol.t) 

### 
using MAT
sensi_Data = matread(datapath*"sensis_PEL_nPFC_BC.mat")

sensis = sensi_Data["sensis"]
sol = sensi_Data["odesol"];
ω0 = 100*pi;
share_pe = 0.3;
xcpu = 0.036; 
Cpu = 1/(xcpu*ω0) * share_pe
ΔCd = 0.1*Cpu
sol_appr = ApproximatedTrajectory(sensi_Data["odesol"],sensis[3],ΔCd)
ode_apprx[:,:] = sol_appr;
pgsol_apprx = PowerGridSolution(ode_apprx,pg0);

PlotlyJS.plot([plotv(pgsol0,"bus_load"),plotv(pgsol_apprx,"bus_load"),plotv(pgsol0_per,"bus_load")])
PlotlyJS.plot([myplot(pgsol0,"bus_load",:p1),myplot(pgsol_apprx,"bus_load",:p1),myplot(pgsol0_per,"bus_load",:p1)])
PlotlyJS.plot([myplot(pgsol0,"bus_load",:ps),myplot(pgsol_apprx,"bus_load",:ps),myplot(pgsol0_per,"bus_load",:ps)])
PlotlyJS.plot([myplot(pgsol0,"bus_load",:qs),myplot(pgsol_apprx,"bus_load",:qs),myplot(pgsol0_per,"bus_load",:qs)])


####### Solution from Automatic Differentiation
x, dp = extract_local_sensitivities(pgsol0_ad);
da = dp[6]
pad = PlotlyJS.scatter(x=pgsol0_ad.t,y=da[15,:]);
xts, dpts = extract_local_sensitivities(pgsol0_adts);
dats = dpts[6]
padts = PlotlyJS.scatter(x=pgsol0_adts.t,y=dats[15,:]);
PlotlyJS.plot([pad,padts])


x, dp = extract_local_sensitivities(pgsol0_ad)
da = dp[6]
sol_appr_ad = ApproximatedTrajectory(x[:,:],da,ΔCd)

u5 = hypot.(sol_appr[9,:],sol_appr[10,:]);
p5 = PlotlyJS.scatter(x=pgsol0_ad.ty=u5);

PlotlyJS.plot([plotv(pgsol0,"bus_load"),plotv(pgsol_apprx,"bus_load"),plotv(pgsol0_per,"bus_load"),p5])

x5r = x[9,:] .+ da[9,:]*ΔCd
x5i = x[10,:] .+ da[10,:]*ΔCd
x5 = hypot.(x5r,x5i)
px5 = PlotlyJS.scatter(x=pgsol0_ad.t,y=x5);
PlotlyJS.plot([plotv(pgsol0,"bus_load"),px5,plotv(pgsol_apprx,"bus_load"),plotv(pgsol0_per,"bus_load"),p5])

zx = 19
zs = :qs
xz = x[zx,:] .+ da[zx,:]*ΔCd
px5 = PlotlyJS.scatter(x=pgsol0_ad.t,y=xz);
PlotlyJS.plot([px5,myplot(pgsol_apprx,"bus_load",zs),myplot(pgsol0_per,"bus_load",zs)])
