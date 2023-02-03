using MATLAB
using FileIO
using PlotlyJS

tmp = load("VSM_I070_R_20_1_0_X_15_1_0.jld")

Rv= tmp["Rverlauf"];
Xv = tmp["Xverlauf"];
XRv= tmp["XR"];
XRt= tmp["XR_t"];

px,py = DetermineBoundary(xr,Rverlauf,Xverlauf)
px,py = DetermineBoundary(XRv,Rv,Xv)
write_matfile("droop_I100_own_reduction.mat"; xj=px, yj=py) 
plot(px,py)