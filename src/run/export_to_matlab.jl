using MATLAB
using FileIO
using PlotlyJS

tmp = load("droop_I100_ownred_R_13_02_0_X_11_02_0.jld")

Rv= tmp["Rverlauf"];
Xv = tmp["Xverlauf"];
XRv= tmp["XR"];
XRt= tmp["XR_t"];

solv = Array(Sol2DFonlyVoltages(pgsol0))
#px,py = DetermineBoundary(xr,Rverlauf,Xverlauf)
x,y = DetermineBoundary(xr,Rverlauf,Xverlauf)
write_matfile("droop_LTVS_secm5.mat"; sol = solv) 
plot(x,y)