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
sel_params = collect(1:16)#[2,12,13,15,16]

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
xlims!(tmp-0.1,tmp+0.1)
z1[12083,1]
for (ind,val) in enumerate(eachcol(z1[:,2:end]))
    display(plot!(pgsol.dqsol.t[1:end-1],val,xlims = (5,70),ylims = (-1,10), label = labels_p[sel_params[ind+1]]))
end

plot(pgsol.dqsol.t[1:end-1],first(eachcol(z2)),xlims = (5,70),ylims = (-10,1), label = labels_p[sel_params[1]], legend=(0.1,0.9))
for (ind,val) in enumerate(eachcol(z2[:,2:end]))
    display(plot(pgsol.dqsol.t[1:end-1],val,xlims = (5,70),ylims = (-12,1), label = labels_p[sel_params[ind+1]]))
end

plot(pgsol.dqsol.t[1:end-1],first(eachcol(z1)),xlims = (5,70),ylims = (-10,10), label = labels_p[sel_params[1]], legend=(0.1,0.9))
plot!(pgsol.dqsol.t[1:end-1],first(eachcol(z2)),xlims = (5,70), legend=(0.1,0.9))
for ind in sel_params
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
xt090 = ["0","15","30","45","60", L"\textit{t}/s \rightarrow", "90"]
yt07_1 = [ "0.70", "0.80", L"|\underline{\textit{u}}|/p.u.↑", "1.0"]
#Plot 1-1 Verlauf Spannung + wenn droop deaktiviert
#using UnicodeFun
#using LatexSVG
#using PGFPlots
using Plots; pgfplotsx()

begin
    plot(pgsol,"bus4",:v,  label = L"|\underline{u}_\mathrm{f}|" ,color =grey3)
    plot!(pgsol,"bus3",:v, label = L"|\underline{u}_\mathrm{load}|", color =grey2)
    plot!(pgsol,"bus2",:v, label = L"|\underline{u}_\mathrm{PCC}|",color =grey1)


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
         legend = (0.1,0.55),
         framestyle = :box,
         linewidth = 0.7,
         size = (500,220), #(600,270)
         legendfont=font(10,"Times",halign=:left),
         labelfontsize = 10,
         xtickfont = font(10, "Times"),
         font ="Times",
         ytickfont = font(10, "Times"),
         grid = true,
         gridalpha = 0.15,
         gridstyle = :dash,
         xticks = (0:15:90, xt090),
         yticks = (0.7:0.1:1.0, yt07_1),
         left_margin = 0mm,
         bottom_margin = 0mm,
         right_margin = 2.0mm,
         xlabel = "",
         top_margin = 0mm)

         vline!([62.5,62.5],color=:black,linestyle=:dash,label="")
         quiver!([53.5], [0.725], quiver=([9], [0.025]),color=:black,label="")
         annotate!([42],[0.725],"droop deactivated",font(10, "Times"))
end
#savefig("C:\\Users\\liemann\\Desktop\\test.svg")
savefig("C:\\Users\\liemann\\Desktop\\PSCC-EMT-Modell\\LTVS_droop_new.svg")

save
#Plot 1-2 U-Sensis der interessanten Parameter
begin #bei 12083 endet es bei den Sensis, bei 12084 startet es
    plot(pgsol,"bus4",:v,label="base case",xlims=(5,72),ylims=(0.9,1.0))
    #plot!(pgsol_kp,"bus4",:v,label="Kp")
    #plot!(pgsol_kq,"bus4",:v,label="Kq")
    #plot!(pgsol_imax,"bus4",:v,label="imax")
    #plot!(pgsol_kvi,"bus4",:v,label="kvi")
    #plot!(pgsol_kvq,"bus4",:v,label="kvq")
    plot!(pgsol_imaxcsa,"bus4",:v,label="csa")
end
begin
    xt070 = ["10","20","30","40","50",L"\textit{t}/s \rightarrow", "70"]
    ytm10_0= [ L"-10", L"-7.5", L"-5.0", L"\textit{x}\textsubscript{uf}/p.u.↑", "0.0"] #x_\mathrm{uf}/p.u. ↑
    ytm5_10= [ "-5.0","-2.50","0.0", "2.5", "5.0",L"\textit{x}\textsubscript{uf}/p.u.↑", "10"]
    ind1 = 12083
    uabs_sens = Array{Float64}(undef,length(pgsol.dqsol.t[1:end-1]),6)
    #Kp
    uabs_sens[1:ind1,1] = z1[1:ind1,1]
    uabs_sens[ind1+1:end,1] = z2[ind1+1:end,1]
    #Kq
    uabs_sens[1:ind1,2] = z2[1:ind1,2]
    uabs_sens[ind1+1:end,2] = z1[ind1+1:end,2]
    #imax
    uabs_sens[1:ind1,3] = z1[1:ind1,12]
    uabs_sens[ind1+1:end,3] = z2[ind1+1:end,12]
    #Kvi
    uabs_sens[1:ind1,4] = z1[1:ind1,13]
    uabs_sens[ind1+1:end,4] = z1[ind1+1:end,13]
    #Kvq
    uabs_sens[1:ind1,5] = z1[1:ind1,15]
    uabs_sens[ind1+1:end,5] = z1[ind1+1:end,15]
    #imax_csa
    uabs_sens[1:ind1,6] = z1[1:ind1,16]
    uabs_sens[ind1+1:end,6] = z1[ind1+1:end,16]

    labels_p = [
        L"k_\mathrm{p}", #1
        L"k_\mathrm{q}", #2
        "ωf_P", #3
        "ωf_Q", #4
        "xlf", #5
        "rf", #6
        "xcf", #7
        "Kp_u", #8
        "Ki_u", #9
        "Kp_i", #10
        "Ki_i", #11
        L"i_\mathrm{VImax}", #12
        L"k_\mathrm{pVI}", #13
        "σXR", #14
        L"k_\mathrm{vq}", #15
        L"i_\mathrm{CSAmax}",#16
        ]
    tmp = filter(x->x>-11,uabs_sens[:,1])
    plot(pgsol.dqsol.t[1:length(tmp)],
        tmp,
        xlims = (10,71),
        ylims = (-10,1),
        label = labels_p[1],
        legend=(0.2,0.4),
        color=:black,
        linestyle=:solid
      )
    tmp = filter(x->x>-11,uabs_sens[:,3])
    p1 = plot!(pgsol.dqsol.t[1:length(tmp)],tmp,label = L"i_\mathrm{VImax}" ,color=color_mat,linestyle=:solid,
        xticks = (10:10:70, xt070),
        yticks = (-10:2.5:0,ytm10_0),
        xlabel = "",
        )

    tmp = filter(x->x<11,uabs_sens[:,2])
    plot( pgsol.dqsol.t[1:length(tmp)],tmp,label = labels_p[2],color=:black,linestyle=:solid)
    tmp = filter(x->x<11,uabs_sens[:,6])
    plot!(pgsol.dqsol.t[1:length(tmp)],tmp,label = labels_p[16],color=:black,linestyle=:dash)
    tmp = filter(x->x<11,uabs_sens[:,5])
    plot!(pgsol.dqsol.t[1:length(tmp)],tmp,label = labels_p[15],color=color_mat,linestyle=:dash)
    tmp = filter(x->x<11,uabs_sens[:,4])
    p2 = plot!(pgsol.dqsol.t[1:length(tmp)],tmp,label = labels_p[13],color=color_mat,linestyle=:solid,
        xticks = (10:10:70, xt070),
        yticks = (-5:2.5:10, ytm5_10),
        xlabel = "",
        xlims = (10,71),
        ylims = (-5,10),
        legend=(0.2,0.9),
        )

    l = @layout [a b]
    plot(p1,p2, layout = l,
        framestyle = :box,
        linewidth = 0.7,
        size = (500,220), #(600,270)
        legendfont=font(10,"Times",halign=:left),
        labelfontsize = 10,
        xtickfont = font(10, "Times"),
        font ="Times",
        ytickfont = font(10, "Times"),
        grid = true,
        gridalpha = 0.15,
        gridstyle = :dash,
        left_margin = 0mm,
        bottom_margin = 0mm,
        right_margin = 2.0mm,
        top_margin = 0mm
    )
end
savefig("C:\\Users\\liemann\\Desktop\\PSCC-EMT-Modell\\LTVS_sensis_new.svg")


#Plot 2-1 approximierte Trajektorie für den sensitivsten Parameter
begin #Kvi = 0.055*1.2
    xt070 = ["10","20","30","40","50", "60","70", L"\textit{t}/s \rightarrow", "90"]
    yt09_10= [ "0.85","0.90",  L"|\underline{\textit{u}}\textsubscript{f}|/p.u.↑", "1.0"]

    plot(pgsol,"bus4",:v,label="base case",xlims=(5,90),ylims=(0.85,1.0), color = grey2)
    plot!(pgsol_kvi,"bus4",:v,label="perturbed", color = color_mat)

    sens_ur, sens_ui = GetVoltageSensis(toll_tap,7,8)
    ur_or = ExtractResult(pgsol,:u_r_4)[1:end-1]
    ui_or = ExtractResult(pgsol,:u_i_4)[1:end-1]

    p_old = pgsol.powergrid.nodes["bus4"].Kvi
    p_new = pgsol_kvi.powergrid.nodes["bus4"].Kvi
    #ur_appr = ur_or .+ sens_ur[:,15]*(p_new-p_old)
    #ui_appr = ui_or .+ sens_ui[:,15]*(p_new-p_old)
    #u_abs = hypot.(ur_appr,ui_appr)
    #plot!(pgsol.dqsol.t[1:end-1],u_abs,label="Kq -approximated",linestyle=:dash)#,xlims=(5,72),ylims=(-3,3)

    uabs = hypot.(ur_or,ui_or) .+ uabs_sens[:,4]*(p_new-p_old)
    tmp = filter(x->x<1.1,uabs)
    p2 = plot!(pgsol.dqsol.t[1:length(tmp)],
    tmp,
    label="approximated",
    linewidth = 0.7,
    linestyle=:dash,
    framestyle = :box,
    legend = (0.2,0.55),
    color = :black,
    size = (500,220),
    legendfont=font(10,"Times",halign=:left),
    labelfontsize = 10,
    xtickfont = font(10, "Times"),
    font ="Times",
    ytickfont = font(10, "Times"),
    grid = true,
    gridalpha = 0.15,
    gridstyle = :dash,
    xticks = (10:10:90, xt070),
    yticks = (0.85:0.05:1.0, yt09_10),
    xlabel = "",
    )
end
savefig("C:\\Users\\liemann\\Desktop\\PSCC-EMT-Modell\\LTVS_appr_new.svg")

#Plot 2-2 Original, ein deutlich verbesserter und ein deutlich verschlechter Verlauf
iabs_t = iabs[12082+3:end]
plot(iabs_t)
plot(pgsol.dqsol.t[1:end-1],uabs_sens[:,1],xlims=(5,72),ylims=(-10,5))


sol_sensi_per = deepcopy(pgsol.dqsol)
for (ind,val) in enumerate(collect(eachcol(toll_tap[15])))
    sol_sensi_per.u[ind] .+= val*(0.005)
end
pgsol_appr = PowerGridSolution(sol_sensi_per,pg)
plot!(pgsol_appr,"bus4",:v)
