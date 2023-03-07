using MATLAB
using FileIO
using PlotlyJS

tmp = load("droop_I100_ownred_R_13_02_0_X_11_02_0.jld")

Rv= tmp["Rverlauf"];
Xv = tmp["Xverlauf"];
XRv= tmp["XR"];
XRt= tmp["XR_t"];

#px,py = DetermineBoundary(xr,Rverlauf,Xverlauf)
px,py = DetermineBoundary(XRv,Rv,Xv)
write_matfile("droop_I100_own_reduction_absprio_fein.mat"; xj=px, yj=py) 
plot(px,py)