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
sel_params = collect(1:16)#[1,12,13,15,16]

toll = load("C:/Users/liemann/Desktop/sens_pscc_kq_1em1_dt_1em2.jld")
toll = toll["sens"]

toll_tap = load("C:/Users/liemann/Desktop/sens_pscc_kq_1em1_dt_1em2_tap_param.jld")
toll_tap = toll_tap["sens"]

plot(pgsol.dqsol.t[1:end-1],toll[sel_params[1]][look_on,1:end], title = "Sensis of $(String(syms[look_on]))",
    label = labels_p[sel_params[1]],
    legend = :outertopright,
    #size = (1000,750),
    #xlims = (1.5,65.0),
    xlims = (5,70),
    ylims = (-2,2))
    #ylims = (-1.5,1.5))
for i in sel_params[2:end]
    display(plot!(pgsol.dqsol.t[1:end-1],toll[i][look_on,1:end], label = labels_p[i]))
end


plot(pgsol.dqsol.t[1:end-1],toll_tap[sel_params[1]][look_on,1:end], title = "Sensis of $(String(syms[look_on]))",
    label = labels_p[sel_params[1]],
    legend = :outertopright,
    #size = (1000,750),
    #xlims = (1.5,65.0),
    xlims = (5,72),
    linestyle = :dash,
    ylims = (-2,2))
    #ylims = (-1.5,1.5))
for i in sel_params[2:end]
    display(plot!(pgsol.dqsol.t[1:end-1],toll_tap[i][look_on,1:end], label = labels_p[i],linestyle = :dash))
    sleep(2)
end

z1,z2 = GetAbsVoltageSensis(pgsol,:u_r_4,:u_i_4,toll_tap,sel_params,ones(length(sel_params)))

plot(pgsol.dqsol.t[1:end-1],first(eachcol(z1)),xlims = (5,70),ylims = (-1,10), label = labels_p[sel_params[1]], legend=(0.1,0.9))
for (ind,val) in enumerate(eachcol(z1[:,2:end]))
    display(plot!(pgsol.dqsol.t[1:end-1],val,xlims = (5,70),ylims = (-1,10), label = labels_p[sel_params[ind+1]]))
end

plot(pgsol.dqsol.t[1:end-1],first(eachcol(z2)),xlims = (5,70),ylims = (-10,1), label = labels_p[sel_params[1]], legend=(0.1,0.9))
for (ind,val) in enumerate(eachcol(z2[:,2:end]))
    display(plot(pgsol.dqsol.t[1:end-1],val,xlims = (5,70),ylims = (-12,1), label = labels_p[sel_params[ind+1]]))
end

plot(pgsol.dqsol.t[1:end-1],first(eachcol(z1)),xlims = (5,70),ylims = (-10,10), label = labels_p[sel_params[1]], legend=(0.1,0.9))
plot!(pgsol.dqsol.t[1:end-1],first(eachcol(z2)),xlims = (5,70), legend=(0.1,0.9))
for ind in 2:16
    display(plot(pgsol.dqsol.t[1:end-1],z1[:,ind],xlims = (5,70),ylims = (-10,10), label = labels_p[ind]))
    display(plot!(pgsol.dqsol.t[1:end-1],z2[:,ind]))
end


sens_ur,sens_ui = GetVoltageSensis(toll_tap,7,8)

sens_ur = sens_ur[:,16]
sens_ui = sens_ui[:,16]


plot(pgsol.dqsol.t[1:end-1],z1[:,16],xlims = (5,70),ylims = (-10,10), label = "z1",legend=(0.1,0.9))
plot!(pgsol.dqsol.t[1:end-1],z2[:,16], label="z2")
plot!(pgsol.dqsol.t[1:end-1],sens_ur,label ="sens_ur")
plot!(pgsol.dqsol.t[1:end-1],sens_ui,label ="sens_ui")
plot!(pgsol.dqsol.t[1:end-1],hypot.(sens_ur,sens_ui),label ="hypot")
plot!(pgsol.dqsol.t[1:end-1],sens_ur+sens_ui,label ="sub")


using LaTeXStrings
using Plots.PlotMeasures
using Measures
color_mat = RGB(0.556863,0.8,0.419608)
grey1 = RGB(0.1,0.1,0.1)
grey2 = RGB(0.5,0.5,0.5)
grey3 = RGB(0.7,0.7,0.7)
font_size = 10
xt090 = ["0","15","30","45","60", L"t/s \rightarrow", "90"]
yt07_1 = [ "0.70", "0.80", L"$|\underline{u}|/p.u. \uparrow$", "1.0"]
#Plot 1-1 Verlauf Spannung + wenn droop deaktiviert
begin

    plot(pgsol,"bus4",:v, label = L"$|u_\mathrm{f}|$",color =grey3)
    plot!(pgsol,"bus3",:v, label = L"$|u_\mathrm{load}|$", color =grey2)
    plot!(pgsol,"bus2",:v, label = L"|u_\mathrm{PCC}|",color =grey1)


    plot!(pgsol_droop,"bus2",:v, label = "",color =grey1,linestyle=:dash)
    plot!(pgsol_droop,"bus3",:v, label = "", color =grey2,linestyle=:dash)
    plot!(pgsol_droop,
        "bus4",
        :v,
         label = "",
         color =grey3,
         linestyle=:dash,
         xlims=(0,90),
         ylims=(0.7,1.0),
         legend = (0.4,0.4),
         framestyle = :box,
         linewidth = 1.0,
         #size = (700,250), #(600,270)
         legendfont=font(10,"Times"),
         labelfontsize = 10,
         xtickfont = font(10, "Times"),
         font ="Times",
         ytickfont = font(10, "Times"),
         grid = true,
         gridalpha = 0.25,
         gridstyle = :dash,
         xticks = (0:15:90, xt090),
         yticks = (0.7:0.1:1.0, yt07_1),
         left_margin = 0mm,
         bottom_margin = 0mm,
         right_margin = 2.0mm,
         top_margin = 0mm)
         vline!([62.5,62.5],color=:black,linestyle=:dash,label="")
         quiver!([52.5], [0.725], quiver=([10], [0.025]),color=:black,)
         annotate!([42],[0.725],"droop deactivated",font(10, "Times"))
end



savefig("C:\\Users\\liemann\\Desktop\\test.svg")

save
#Plot 1-2 U-Sensis der interessanten Parameter

#Plot 2-1 approximierte Trajektorie für den sensitivsten Parameter

#Plot 2-2 Original, ein deutlich verbesserter und ein deutlich verschlechter Verlauf
