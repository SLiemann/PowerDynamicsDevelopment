using PowerDynamics
using Plots
using DifferentialEquations
using JLD

function test()
    #include("C:/Users/liemann/github/PowerDynamicsDevelopment/src/include_costum_nodes_lines_utilities.jl")
    include("C:/Users/liemann/github/PowerDynamicsDevelopment/src/operationpoint/PowerFlow.jl")
    include("C:/Users/liemann/github/PowerDynamicsDevelopment/src/operationpoint/InitializeInternals.jl")
    include("C:/Users/liemann/github/PowerDynamicsDevelopment/src/grids/LTVS_Test_System.jl")
    include("C:/Users/liemann/github/PowerDynamicsDevelopment/src/utility/utility_functions.jl")
    include("C:/Users/liemann/github/PowerDynamicsDevelopment/src/sensitivity_analyses/Module_Local_Sensitivity.jl")
    return nothing
end

function CalcLTVSSensis(;time = 2e-6,n=2)
    pg = GFC_LTVS_Test_System(nTap = 5)
    Qmax   = [Inf, Inf, Inf,Inf, Inf,Inf*sqrt(1-0.95^2)]
    Qmin   = -Qmax
    U,δ,ic0 = PowerFlowClassic(pg,iwamoto = true, Qmax = Qmax, Qmin = Qmin, Qlimit_iter_check = 2,max_tol = 1e-6)
    pg, ic0 = InitializeInternalDynamics(pg,ic0)

    pgsol,evr  = run_LTVS_simulation(pg,ic0,(0.0,time))

    mtk_normal = GetMTKLTVSSystem(pg_state = "gfc_normal")
    mtk_fault = GetMTKLTVSSystem(pg_state = "gfc_fault")
    mtk_postfault = GetMTKLTVSSystem(pg_state = "gfc_postfault")
    mtk = [mtk_normal; mtk_fault; mtk_postfault]

    s = GetTriggCondsLTVS(mtk_normal)
    h = GetStateResFunLTVS(mtk_normal)
    p_pre = GFC_LTVS_params()
    sensis_p = collect(1:n)
    #@time toll = Main.MyLocalSensi.CalcHybridTrajectorySensitivity(mtk,pgsol.dqsol,p_pre,evr,s,h,[],sensis_p)
    @time InitCalcHybrid(mtk,pgsol.dqsol,p_pre,evr,s,h,[],sensis_p)
end

function InitCalcHybrid(
    mtk::Vector{ODESystem},
    sol::ODESolution,
    p_pre::Vector{Float64},
    evr::Matrix{Float64},
    s::Vector{Equation},
    h::Vector{Matrix{Equation}},
    u0_sensi::Vector{Union{Int64,Any}},
    p_sensi::Vector{Int64},
  )
  mtk0 = mtk[1] # it is assumed that the first element is the initial system
  ic = sol.prob.u0
  @time xx0_k, yx0_k, sym_states,sym_params, A_states, D_states, M, N, O, symp, Δt,len_sens, f, g, J =
      Main.MyLocalSensi.InitTrajectorySensitivity(mtk0, ic, p_pre, u0_sensi, p_sensi)
  xx0 = [i[1] for i in xx0_k]
  yx0 = [i[1] for i in yx0_k]
  @time f_all,g_all,J_all,M_all,N_all,O_all = Main.MyLocalSensi.GetEqsJacobianSensMatrices(mtk,xx0,yx0,u0_sensi,p_sensi)
  Fx_all, Fy_all, Gx_all, Gy_all = J_all
  hx,hy,sx,sy = Main.MyLocalSensi.CalcTriggerAndStateResetJacobians(mtk0,s,h)
  sensis = Vector{Array{Float64}}(undef, len_sens)
  for i = 1:length(sensis)
    sensis[i] = Array{Float64}(
      undef,
      size(D_states)[1] + size(A_states)[1],
      size(sol)[2] - 1,
    )
  end
  ind_sol = vcat(1,setdiff(indexin(evr[:,1],sol.t),[nothing]),length(sol.t))
  #ind_sol = [1]
  #for i in evr[:,1] # DifferentialEquations.jl has multiple time points
  #    ind_sol = vcat(ind_sol,findall(x->x==i,sol.t)[end])
  #end
  #ind_sol = vcat(ind_sol,length(sol.t))

  Fx_pre = Fx_all[1]
  Fy_pre = Fy_all[1]
  Gx_pre = Gx_all[1]
  Gy_pre = Gy_all[1]
  J_pre = Fx_pre, Fy_pre, Gx_pre, Gy_pre

  f_pre = f_all[1]
  g_pre = g_all[1]
  Mc = copy(M)
  return nothing
end


mtk_fault = GetMTKLTVSSystem(pg_state = "gfc_fault")
eqs, aeqs, D_states, A_states = Main.MyLocalSensi.GetSymbolicEquationsAndStates(mtk_fault)
@time Main.MyLocalSensi.GetParameterJacobian(mtk_fault,parameters(mtk_fault),states(mtk_fault)[[]],16,eqs,aeqs)

@parameters xx0[1:7, 1:16]
@parameters yx0[1:12, 1:16]
xx0 = Symbolics.scalarize(xx0)
yx0 = Symbolics.scalarize(yx0)

@time Main.MyLocalSensi.GetEqsJacobianSensMatrices(mtk,xx0,yx0,[],collect(1:16))
