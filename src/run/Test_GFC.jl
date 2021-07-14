using PowerDynamics
using DifferentialEquations
using Plots
using JLD
#using ModelingToolkit
#using CSV
#using DataFrames

begin
    include("C:/Users/liemann/github/PowerDynamicsDevelopment/src/include_costum_nodes_lines_utilities.jl")
    include("C:/Users/liemann/github/PowerDynamicsDevelopment/src/grids/GFC_Test_Grid.jl")
    include("C:/Users/liemann/github/PowerDynamicsDevelopment/src/sensitivity_analyses/Local_Sensitivity.jl")

    pg = GFC_Test_Grid()
    U,δ,ic0 = PowerFlowClassic(pg, iwamoto = true, max_tol = 1e-7)
    Ykk = NodalAdmittanceMatrice(pg)
    Uc = U.*exp.(1im*δ/180*pi)
    I_c = Ykk*Uc
    S = conj(Ykk*Uc).*Uc

    pg1 ,ic = InitializeInternalDynamics(pg,I_c,ic0)
    params = GFC_params()
    prob = ODEProblem(rhs(pg1),ic,(0.0,5.0),params)
end
@time pgsol,evr = simGFC(prob)

mtk_normal = GetMTKSystem(pg,(0.0,1.0),params)
mtk_fault = GetMTKSystem(GFC_Test_Grid(y_new = zfault()),(0.0,1.0),params)
mtk = [mtk_normal,mtk_fault]
s = GetTriggCondsGFCTest(mtk_normal)
h = GetStateResFunGFCTest(mtk_normal)
@time sens = CalcHybridTrajectorySensitivity(mtk,pgsol.dqsol,params,evr,s,h,[],collect(1:15))
save("C:/Users/liemann/Desktop/Sens_GFC_Test/sens_all_p0set_3_dt_1em2.jld", "sens", sens,"ic0",ic0,"p_pre",params,"evr",evr,"sensis_p",collect(1:15))

plot(pgsol.dqsol.t[1:end-1],sens[1][14,:])

plot!(pgsol,["bus3"],:i_abs, label = "Kvi = " * string(pg.nodes["bus3"].Kvi) * ", K_vq = " *string(pg.nodes["bus3"].K_vq))
ylims!((2.5,2.99))
xlims!((0.99,2.5))
plot(pgsol,collect(keys(pg.nodes))[2:end],:v)
plot(pgsol,["bus3"],:iabs)
plot(pgsol,["bus3"],:p)
plot(pgsol,["bus3"],:q)
plot(pgsol,["bus3"],:θ)
plot(pgsol,["bus3"],:φ)
plot(pgsol,["bus3"],:ω)
plot(pgsol,["bus3"],:Qm)
xlims!((4.9,5.50))
plot(pgsol,["bus3"],:e_uq)
plot(pgsol,["bus3"],:e_ud)

sol = ExtractResult(pgsol,:ω_3)

f, freq = DFT(sol,pgsol.dqsol.t)
plot(freq,abs.(f), xlim=(-10, 10))
#Calculating approximated trajectory
#sol_original = deepcopy(pgsol.dqsol)
#pgsol_or = deepcopy(pgsol)
sol_appr = deepcopy(sol_original)
for (ind,val) in enumerate(collect(eachcol(sens[1])))
    sol_appr.u[ind] .+= val*-(0.003)
end
pgsol_tmp = PowerGridSolution(sol_appr,pg)
plot!(pgsol_tmp,"bus3",:i_abs, label = "approximated")
plot!(pgsol,"bus3",:i_abs, label = "real perturbed")
plot(pgsol_or,"bus3",:i_abs, label = "Original")
plot!(pgsol_tmp,collect(keys(pg.nodes)),:v,legend = fals


#Plot all parameter sensitivities
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
    ]
syms = rhs(pg).syms
look_on = 1
plot(pgsol.dqsol.t[1:end-1],sens[1][look_on,1:end], title = "Sensis of $(String(syms[look_on]))",
    label = labels_p[1],
    legend = :outertopright,
    size = (1000,750))
    #xlims!((1.9,2.3))
for i in collect(3:15)
    display(plot!(pgsol.dqsol.t[1:end-1],sens[i][look_on,1:end], label = labels_p[i]))
    #sleep(3.0)
end


#Calculating
Imax = 1.0
I_tresh = 0.95
ΔXR = 10.0
V0 = 1.0
a = (Imax - I_tresh)^2 * (1 + ΔXR^2)
b = 2.0*(Imax - I_tresh) * (ΔXR)
c = - V0^2/Imax^2

kp = (-b + sqrt(b^2 - 4.0*a*c))/(2*a)


id = 0.95
iq = 0.4
iabs = hypot(id,iq)
ϕ = atan(iq,id)

iabs = 1
ids = iabs*cos(ϕ)
iqs = iabs*sin(ϕ)
ϕ2 = atan(iqs,ids)
