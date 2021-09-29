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
sel_params = [1,12,13,15,16]
param_val = GFC_LTVS_params_TapParam()[sel_params]

plot(pgsol.dqsol.t[1:end-1],toll_tap[sel_params[1]][look_on,1:end], title = "Sensis of $(String(syms[look_on]))",
    label = labels_p[sel_params[1]],
    legend = :outertopright,
    size = (1000,750),
    #xlims = (1.5,65.0),
    xlims = (0.9,1.5),
    ylims = (-25,100))
    #ylims = (-1.5,1.5))
for i in sel_params[2:end]
    display(plot!(pgsol.dqsol.t[1:end-1],toll_tap[i][look_on,1:end], label = labels_p[i]))
end

z1,z2 = GetAbsVoltageSensis(pgsol,:u_r_4,:u_i_4,toll_tap,sel_params,ones(length(sel_params))) #param_val

plot(pgsol.dqsol.t[1:end-1],z1)
plot(pgsol.dqsol.t[1:end-1],z2)
#Voltage Approximation
begin
    display(plot(pgsol,"bus4",:v,label = "Original"))
    display(plot!(pgsol_per,"bus4",:v, label ="real perturbed"))
    sens_ur, sens_ui = GetVoltageSensis(toll_tap,7,8)
    ur_or = ExtractResult(pgsol,:u_r_4)[1:end-1]
    ui_or = ExtractResult(pgsol,:u_i_4)[1:end-1]
    ur_appr = ur_or .+ sens_ur[:,1]*(0.02*0.02)
    ui_appr = ui_or .+ sens_ui[:,1]*(0.02*0.02)
    u_abs = hypot.(ur_appr,ui_appr)
    plot!(pgsol.dqsol.t[1:end-1],u_abs,label="approximated",linestyle=:dash)
end

#Current Approximation
plot(pgsol,"bus4",:i_abs,label = "Original")
display(plot!(pgsol_per,"bus4",:i_abs, label ="real perturbed"))
iabs_or = ExtractResult(pgsol,:i_abs_4)[1:end-1]
iabs_appr = iabs_or .+ toll_tap[16][16,:]*(1.1*0.01)
plot!(pgsol.dqsol.t[1:end-1],iabs_appr,label="approximated",linestyle=:dash)


############# PLOT Paper ###################
using LaTeXStrings
using Plots.PlotMeasures
using Measures

default(guidefont = (8, "Times"))
begin
    xt0913 = ["0.90", "1.0","1.10", L"t/s \rightarrow", "1.30"]
    yt06_11 = ["0.60", "0.70", "0.80", "0.90", L"$|\underline{u}_f|/p.u. \uparrow$", "1.10"]
    ytm0_80 = []
    for i in 0.0:20.0:80.0
        append!(ytm0_80,[string(i)])
    end
    ytm0_80[end-1] = L"$x_\mathrm{uf}/p.u. \uparrow$"
    ytm0_m280 = []
    for i in 0.0:-70.0:-280
        append!(ytm0_m280,[string(i)])
    end
    ytm0_m280[2] = L"$x_\mathrm{uf}/p.u. \uparrow$"

    color_mat = RGB(0.556863,0.8,0.419608)
    grey = RGB(0.5,0.5,0.5)
    font_size = 10

    plot(pgsol,"bus4",:v,label = "original",color = grey)
    plot!(pgsol_per,"bus4",:v, label ="perturbed", color = color_mat)
    p1 = plot!(
        pgsol.dqsol.t[1:end-1],
        u_abs,
        legend = (0.615,0.36),
        xlims = (0.9, 1.3),
        xticks = (0.9:0.1:1.3, xt0913),
        ylims = (0.60,1.1),
        yticks = (0.6:0.10:1.10, yt06_11),
        color = :black,
        xlabel = "",
        linestyle = :dash,
        label = "approximated"
    )
    annotate!((0.9,0.063),text("(b)",10,"Times"))


    labels_sens_u = [
        L"$k_\mathrm{p}$",
        L"$i_\mathrm{VImax}$",
        L"$k_\mathrm{pVI}$",
        L"$k_\mathrm{vq}$",
        L"$i_\mathrm{CSAmax}$",
    ]


    plot(
        pgsol.dqsol.t[1:end-1],
        z1[:,5],
        label = labels_sens_u[5],
        legend = (0.2,0.6),
        xlims = (0.9, 1.3),
        xticks = (0.9:0.1:1.3, xt0913),
        yticks = (0:20:80,ytm0_80),
        ylims = (-5,80),
        color = :black
        )
    plot!(pgsol.dqsol.t[1:end-1],z1[:,3],label = labels_sens_u[3],color = color_mat)
    p2 = plot!(pgsol.dqsol.t[1:end-1],z1[:,4],label = labels_sens_u[4],color = grey)
    annotate!((0.9,0.17),text("(a)",10,"Times"))

    plot(
        pgsol.dqsol.t[1:end-1],
        z2[:,2],
        label = labels_sens_u[2],
        legend = (0.2,0.6),
        xlims = (0.9, 1.3),
        xticks = (0.9:0.1:1.3, xt0913),
        yticks = (0:-70:-280.0,ytm0_m280),
        ylims = (-290,10),
        color = color_mat
        )
    p3 = plot!(pgsol.dqsol.t[1:end-1],z2[:,1],label = labels_sens_u[1],color = :black)
    annotate!((0.9,0.17),text("(a)",10,"Times"))

    l = @layout [[b; c] a]
    plot(p2,p3,p1, layout = l,
        framestyle = :box,
        linewidth = 1.0,
        size = (700,250), #(600,270)
        legendfont=font(10,"Times"),
        #xlabel = "\$x_1\$",
        labelfontsize = 10,
        #linestyle = :dot,
        #xlims = (-2.5,2.5),
        #xlims = (-1,1),
        xtickfont = font(10, "Times"),
        #ylabel = "\$x_2\$",
        #ylims = (-2.5,2.5),
        #ylims = (-1,1),
        #linecolor = RGB(146/255,208/255,80/255),
        #linecolor = RGB(244/255,177/255,131/255),
        #linecolor = RGB(0/255,0/255,255/255),
        font ="Times",
        ytickfont = font(10, "Times"),
        grid = true,
        gridalpha = 0.25,
        gridstyle = :dash,
        left_margin = 0mm,
        bottom_margin = 0mm,
        right_margin = 2.0mm,
        top_margin = 0mm
    )
end
savefig("C:\\Users\\liemann\\Desktop\\PSCC-EMT-Modell\\STVS.svg")

begin #Unstable case
    color_mat = RGB(0.556863,0.8,0.419608)
    grey = RGB(0.5,0.5,0.5)
    font_size = 10
    xt0913 = ["0.90", "1.0","1.10", L"t/s \rightarrow", "1.30"]
    yt06_11 = ["0.60", "0.70", "0.80", "0.90",L"$|\underline{u}_f|/p.u. \uparrow$", "1.10"]
    plot(pgsol,"bus4",:v,label = "original", color = grey)
    p1 = plot!(pgsol_un,
        "bus4",
        :v,
        label = L"$k_\mathrm{p}=0.022$",
        legend = (0.37,0.9),
        xlims = (0.9, 1.3),
        xticks = (0.9:0.1:1.3, xt0913),
        ylims = (0.60,1.1),
        yticks = (0.6:0.10:1.10, yt06_11),
        color = :black,
        xlabel = ""
        )
    annotate!((0.9,0.1),text("(a)",10,"Times"))


    plot(pgsol_stkvq,"bus4",:v,label = L"$k_\mathrm{vq}=0.11$", color = color_mat)
    plot!(pgsol_stkvi,"bus4",:v,label = L"$k_\mathrm{pVI}=0.0605$", color = :black, linestyle = :dash)
    p2 = plot!(pgsol,
            "bus4",
            :v,
            label = "base case",
            legend = (0.61,0.4),
            xlims = (0.9, 1.3),
            xticks = (0.9:0.1:1.3, xt0913),
            ylims = (0.60,1.1),
            yticks = (0.6:0.10:1.10, yt06_11),
            color = grey,
            xlabel = ""
            )
    annotate!((0.9,0.1),text("(b)",10,"Times"))

    l = @layout [a b]
    plot(p1,p2, layout = l,
        framestyle = :box,
        linewidth = 1.0,
        size = (700,250),
        legendfont=font(10,"Times"),
        labelfontsize = 10,

        xtickfont = font(10, "Times"),

        font ="Times",
        ytickfont = font(10, "Times"),
        grid = true,
        gridalpha = 0.25,
        gridstyle = :dash,
        left_margin = 0mm,
        bottom_margin = 0mm,
        right_margin = 2.0mm,
        top_margin = 0mm
    )
end
savefig("C:\\Users\\liemann\\Desktop\\PSCC-EMT-Modell\\STVS_un_st.svg")
