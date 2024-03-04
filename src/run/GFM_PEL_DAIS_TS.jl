using PlotlyJS, DataFrames
using ModelingToolkit   

Sbase = 8000e6
Ubase = 400e3
Zbase = Ubase^2/Sbase

begin
    include("C:/Users/liemann/github/PowerDynamicsDevelopment/src/include_custom_nodes_lines_utilities.jl")
    include("C:/Users/liemann/github/PowerDynamicsDevelopment/src/grids/LTVS_Test_System_N32 GFM_PEL_TS.jl")
    nothing
end

pg0,ic0 = Initialize_N32_GFM_PEL_TS();
@time sensi_ad, evr_sol = simulate_LTVS_N32_simulation_N32_GFM_PEL_TS(pg0,ic0,(0.0,80.0),1e6*(20.0+1im*20)/Zbase);
jp = plotallvoltages(pgsol0);
toff = plot(myplot(pgsol0,"bus_load",:toff))
toff = plot(myplot(pgsol0,"bus_load",:ton))
plot(myplot(pgsol0,"bus_load",:q_on))
plot(myplot(pgsol0,"bus_load",:vofft2))
plot(myplot(pgsol0,"bus_load",:voff))
myplot(pgsol0,"bus_load",:ps)
myplot(pgsol0,"bus_load",:q1)

myplot(pgsol0,"bus_gfm",:q_imax)
myplot(pgsol0,"bus_gfm",:i_abs)
plot(myplot(pgsol0,"bus_gfm",:q_imax))
plot(myplot(pgsol0,"bus_gfm",:q_idcmax))
plot(myplot(pgsol0,"bus_gfm",:i_abs))
plot(myplot(pgsol0,"bus_gfm",:idc0))
myplot(pgsol0,"bus_gfm",:P0)
plot([myplot(pgsol0,"bus_gfm",:θ),myplot(pgsol0_per,"bus_gfm",:θ)])
plot([myplot(pgsol0,"bus_gfm",:P0),myplot(pgsol0_per,"bus_gfm",:P0)])

## Compare to EMT-Simulink
using MAT
datapath = "\\\\fs0\\home\\liemann\\Diss\\Results\\623_GFM_PEL\\"
mdata = matread(datapath*"PEL_0_PEL_share_0.3_Zf_20.mat")
vpel = mdata["Vpel"]
mp = [
        PlotlyJS.scatter(x=vpel[:,1],y=vpel[:,2],name="V0 (EMT)"),
        PlotlyJS.scatter(x=vpel[:,1],y=vpel[:,3],name="V1 (EMT)"),
        PlotlyJS.scatter(x=vpel[:,1],y=vpel[:,4],name="V2 (EMT)"),
        PlotlyJS.scatter(x=vpel[:,1],y=vpel[:,5],name="V_GFM (EMT)"),
        PlotlyJS.scatter(x=vpel[:,1],y=vpel[:,6],name="V_load (EMT)")
    ]
append!(mp,jp)
PlotlyJS.plot(mp)


### Saving Sensitivities
using MATLAB
x, dp = extract_local_sensitivities(sensi_ad)
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
    "Cd", #25
    "Pdc", #26
    ];
state_labels = string.(rhs(pg0).syms)
params = GetParamsGFM_GFM_PEL_TS(pg0)
path2 = "C:\\Users\\liemann\\github\\PowerDynamicsDevelopment\\src\\results\\"
write_matfile(path2*"AD_sensis_GFM_PEL_share_0.3_noBC_with_qoffset_80s.mat"; odesol = x[:,:],sensis = dp, sensi_labels=labels_p,state_labels =state_labels,params=params, time = sensi_ad.t) 

write_matfile(path2*"odesol_LTVS_GFM_PEL_share_0.3_noBC_with_qoffset.mat"; odesol = pgsol0.dqsol[:,:],state_labels =state_labels,params=params, time = pgsol0.dqsol.t) 