labels_p = [
    "Kp_droop", #1
    "Kq_droop", #2
    "ωf_P", #3
    "ωf_Q", #4
    "xlf", #5
    "rf", #6
    "xcf", #7
    "Kp_u", #8
    "Ki_u", #9
    "Kp_i", #10
    "Ki_i", #11
    "imax", #12
    "Kvi", #13
    "σXR", #14
    "K_vq", #15
    "imax_csa",#16
]
syms = rhs(pg).syms
look_on = 16
sel_params = [1, 12, 13, 15, 16]
sel_params_low = append!(collect(2:4), append!(collect(8:11), 14))

z1, z2 = GetAbsVoltageSensis(
    pgsol,
    :u_r_4,
    :u_i_4,
    toll_tap,
    sel_params,
    ones(length(sel_params)),
)

z1_low, z2_low = GetAbsVoltageSensis(
    pgsol,
    :u_r_4,
    :u_i_4,
    toll_tap,
    sel_params_low,
    ones(length(sel_params_low)),
)

plot(z2_low[:,8],labels=labels_p[sel_params_low],legend=:bottomright)
plot(z1_low[:,8],labels=labels_p[sel_params_low],legend=:topleft)

plot(pgsol, "bus4", :v, label = "base case")
plot!(pgsol_per, "bus4", :v, label = "real perturbed")
