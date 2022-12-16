using MATLAB
using FileIO
using PlotlyJS

tmp = load("droop_I00_250_R_10_1_0_X_10_1_0.jld")

Rv= tmp["Rverlauf"];
Xv = tmp["Xverlauf"];
XRv= tmp["XR"];
XRt= tmp["XR_t"];

px,py = plotxkrit(reverse(XRv),collect(reverse(Rv)),collect(reverse(Xv)));

write_matfile("droop_I00_250.mat"; xj=px, yj=py) 
plot(px,py)