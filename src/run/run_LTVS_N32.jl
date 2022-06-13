using PowerDynamics
#using Plots

begin
    include("C:/Users/liemann/github/PowerDynamicsDevelopment/src/include_custom_nodes_lines_utilities.jl")
    include("C:/Users/liemann/github/PowerDynamicsDevelopment/src/grids/LTVS_Test_System_N32.jl")
end
begin
    pg = LTVS_Test_System_N32()
    Qmax   = [Inf, Inf, Inf,Inf,Inf*53.0*sqrt(1-0.837735^2),Inf]
    Qmin   = -Qmax
    U,δ,ic0 = PowerFlowClassic(pg,iwamoto = true, Qmax = Qmax, Qmin = Qmin, Qlimit_iter_check =3,max_tol = 1e-7,iter_max = 30)
    #Uc = U.*exp.(1im*δ/180*pi
    #Ykk = NodalAdmittanceMatrice(pg)
    #Ic = abs.(Ykk*Uc./5.5)
    #S  = round.(Uc.*(conj.(Ykk)*conj.(Uc)),digits=3)
    #pg, ic0 = InitializeInternalDynamics(pg,ic0)

    #pg, ic0 = GetInitializedLTVSSystem()
    #plot(pgsol,"bus4",:i_abs,label = "I-Original",xlims=(5,70),ylims=(1.02,1.11), legend = (0.5,0.1))
    #display(plot!(pgsol_per,"bus4",:i_abs, label ="real perturbed"))
    #display(plot(pgsol_stkvi,"bus4",:v,label = "Enhanced"))
    #display(plot!(pgsol_stkvq,"bus4",:v,label = "Enh kvq"))
    #display(plot(pgsol,"bus4",:v,label = "U-Original",xlims=(5,70),ylims=(0.9,1.01), legend = (0.5,0.1)))
    #display(plot!(pgsol_per,"bus4",:v,label = "real perturbed"))
    #display(plot!(pgsol_per,"bus4",:v, label ="real perturbed")) #linestyle = :dash
end

Uc = U.*exp.(1im*δ/180*pi)
Ykk = NodalAdmittanceMatrice(pg)
S  = round.(Uc.*(conj.(Ykk)*conj.(Uc)),digits=5)
abs.(S)
