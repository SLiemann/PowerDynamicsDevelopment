using PowerDynamics
using ModelingToolkit
using LinearAlgebra

function plot_sensi(time,sensis)
   p = plot(layout = (size(sensis)[1],1))
   for (index, value) in enumerate(sensis)
      plot!(time,value',subplot=index)
   end
   return p
end

function my_rhs(x::Equation)
   x.rhs
end

function my_lhs(x::Equation)
   x.lhs
end

function TimeDomainSensitivies(
    pg::PowerGrid,
    time_interval ::Tuple{Float64,Float64},
    ic::Array{Float64,1},
    p::Array{Float64,1},
    sensis_u0_pd::Array{Int64,1},
    sensis_p_pd::Array{Int64,1},
  )
    prob = ODEProblem(rhs(pg),ic,time_interval,p)
    sol  = solve(prob,Rodas4())
    TimeDomainSensitivies(pg,time_interval,ic,p,sensis_u0_pd,sensis_p_pd,sol)
end

function TimeDomainSensitivies(
  pg::PowerGrid,
  time_interval::Tuple{Float64,Float64},
  ic::Array{Float64,1},
  p::Array{Float64,1},
  sensis_u0_pd::Array{Int64,1},
  sensis_p_pd::Array{Int64,1},
  sol::ODESolution,
)
  mtsys   = GetMTKSystem(pg,time_interval,ic,p)
  TimeDomainSensitivies(mtsys,time_interval,ic,p,sensis_u0_pd,sensis_p_pd,sol)
end

function TimeDomainSensitivies(
  mtsys::ODESystem,
  time_interval::Tuple{Float64,Float64},
  ic::Array{Float64,1},
  p::Array{Float64,1},
  sensis_u0_pd::Array{Int64,1},
  sensis_p_pd::Array{Int64,1},
  sol::ODESolution,
)
  fulleqs = equations(mtsys)
  state = states(mtsys)
  params = parameters(mtsys)
  eqs, aeqs, D_states, A_states = GetSymbolicEquationsAndStates(fulleqs, state)

  #it is assumed that state and rhs(powergrid).syms have the same order
  #sensis_u0 = state[indexin(sensis_u0_pd, rhs(pg).syms)]
  sensis_u0 = state[sensis_u0_pd]
  #sensis_p_pd is here a list with indices of the parameters p
  sensis_p = params[sensis_p_pd]

  #dict from states and parameters with their starting values
  symu0 = state .=> ic
  symp = params .=> p

  Fx, Fy, Gx, Gy = GetSymbolicFactorizedJacobian(eqs, aeqs, D_states, A_states)

  Diff_u0 = Differential.(sensis_u0)
  Diff_p = Differential.(sensis_p)
  len_sens = size(sensis_u0)[1] + size(sensis_p)[1]
  Fp = Array{Num}(undef, size(eqs)[1], len_sens)
  Gp = Array{Num}(undef, size(aeqs)[1], len_sens)

  Fp[:, 1:size(Diff_u0)[1]] .= Num(0)
  Gp[:, 1:size(Diff_u0)[1]] .= Num(0)
  for (ind, val) in enumerate(Diff_p)
    Fp[:, ind+size(Diff_u0)[1]] =
      Num.(expand_derivatives.(map(val, my_rhs.(eqs))))
    Gp[:, ind+size(Diff_u0)[1]] =
      Num.(expand_derivatives.(map(val, my_rhs.(aeqs))))
  end

  @parameters Δt
  @parameters xx0[1:size(D_states)[1], 1:len_sens] #xx0 are the sensitivities regargind differential states
  @parameters yx0[1:size(A_states)[1], 1:len_sens] #yx0 are the sensitivities regargind algebraic states
  M = [
    Δt/2*Fx-I Δt/2*Fy
    Gx Gy
  ]
  N =
    isempty(aeqs) ? -xx0 - Δt / 2 * (Fx * xx0 + Fy * yx0 + Fp) :
    [
      -xx0 - Δt / 2 * (Fx * xx0 + Fy * yx0 + Fp)
      zeros(size(A_states)[1], len_sens)
    ]
  O = [
    -Δt / 2 * Fp
    -Gp
  ]

  #Initialisierung: xx0 enthält die Sensis für x0 und p für x
  xx0_k = xx0 .=> 0.0
  xx0_f = zeros(size(xx0)[1], len_sens)
  ind = setdiff(indexin(sensis_u0, D_states),[nothing])
  for i = 1:length(ind)
    xx0_k[i, ind[i]] = xx0_k[i, ind[i]][1] => 1.0
    xx0_f[i, ind[i]] = 1.0
  end
  # Bei den Sensis für y werden zuerst die dy/x0 Sensi initialisiert
  Gy_float = Substitute(Gy, [symu0; symp])
  # for increasing calculation of inv(Gy)
  yx0_t0 = -inv(Gy_float) * (Gx * xx0_f[:, 1:size(sensis_u0)[1]])
  yp_t0 =
    -inv(Gy_float) *
    (Gp * vcat(zeros(size(sensis_u0)[1], size(sensis_p)[1]), I))
  yx0_k = yx0 .=> Substitute([yx0_t0 yp_t0], [symu0; symp])

  sensi = Vector{Array{Float64}}(undef, len_sens)
  for i = 1:length(sensi)
    sensi[i] = Array{Float64}(
      undef,
      size(D_states)[1] + size(A_states)[1],
      size(sol)[2] - 1,
    )
  end

  for i = 1:size(sol)[2]-1
    dt = sol.t[i+1] - sol.t[i]

    uk = state .=> sol.u[i]
    uk1 = state .=> sol.u[i+1]

    Mfloat = Substitute(M, [uk1; symp; Δt => dt])
    Nfloat = Substitute(N, [uk; symp; vec(xx0_k); vec(yx0_k); Δt => dt])
    Ofloat = Substitute(O, [uk1; symp; Δt => dt])
    res = inv(Mfloat) * (Nfloat + Ofloat)

    for j = 1:length(sensi)
      ind_y = setdiff(indexin(A_states, state), [nothing])
      ind_x = setdiff(indexin(D_states, state), [nothing])
      sensi[j][ind_x, i] = res[1:size(D_states)[1], j]
      sensi[j][ind_y, i] = res[size(D_states)[1]+1:end, j]
    end

    xx0_k = xx0 .=> res[1:size(D_states)[1], :]
    yx0_k = yx0 .=> res[size(D_states)[1]+1:end, :]
  end
  return sensi
end

GetSymbolicEquationsAndStates(mtsys::ODESystem) =
  GetSymbolicEquationsAndStates(equations(mtsys), states(mtsys))
function GetSymbolicEquationsAndStates(
  fulleqs::Array{Equation,1},
  state::Vector{Term{Real, Base.ImmutableDict{DataType, Any}}},
)
  aeqs = Vector{Equation}()
  eqs = Vector{Equation}()
  A_states = similar(state,0)
  D_states = similar(state,0)
  for (index, value) in enumerate(fulleqs)
    if my_lhs(value) !== 0
      push!(eqs, value)
      push!(D_states, state[index])
    elseif my_lhs(value) === 0
      push!(aeqs, value)
      push!(A_states, state[index])
    else
      error("Can not interprete LHS of equation; $value")
    end
  end
  return eqs, aeqs, D_states, A_states
end

GetSymbolicStates(mtsys::ODESystem) = GetSymbolicStates(equations(mtsys), states(mtsys))
function GetSymbolicStates(
  fulleqs::Array{Equation,1},
  state::Vector{Term{Real, Base.ImmutableDict{DataType, Any}}},
  )
  A_states = similar(state,0)
  D_states = similar(state,0)
  for (index, value) in enumerate(fulleqs)
    if my_lhs(value) !== 0
      push!(D_states, state[index])
    elseif my_lhs(value) === 0
      push!(A_states, state[index])
    else
      error("Can not interprete LHS of equation; $value")
    end
  end
  return D_states, A_states
end

GetSymbolicEquations(mtsys::ODESystem) = GetSymbolicEquations(equations(mtsys))
function GetSymbolicEquations(fulleqs::Array{Equation,1})
  aeqs = Vector{Equation}()
  eqs = Vector{Equation}()

  for (index, value) in enumerate(fulleqs)
    if my_lhs(value) !== 0
      push!(eqs, value)
    elseif my_lhs(value) === 0
      push!(aeqs, value)
    else
      error("Can not interprete LHS of equation; $value")
    end
  end
  return eqs, aeqs
end

function GetJacobian(
  eqs::Array{Equation},
  states::Vector{Term{Real, Base.ImmutableDict{DataType, Any}}},
)
    Fx = Array{Num}(undef, size(eqs)[1], size(states)[1])
    Diff_states = Differential.(states)
    for (ind, val) in enumerate(Diff_states)
      Fx[:, ind] = Num.(expand_derivatives.(map(val, my_rhs.(eqs))))
    end
    return Fx
end

function GetSymbolicFactorizedJacobian(
  eqs::Array{Equation,1},
  aeqs::Array{Equation,1},
  D_states::Vector{Term{Real, Base.ImmutableDict{DataType, Any}}},
  A_states::Vector{Term{Real, Base.ImmutableDict{DataType, Any}}},
)
  Fx = GetJacobian(eqs,D_states)
  Fy = GetJacobian(eqs,A_states)

  Gx = GetJacobian(aeqs,D_states)
  Gy = GetJacobian(aeqs,A_states)

  return Fx, Fy, Gx, Gy
end

function GetSymbolicFactorizedJacobian(mtsys::ODESystem)
  fulleqs = equations(mtsys)
  symstates = states(mtsys)
  x_eqs, y_eqs, x, y = GetSymbolicEquationsAndStates(fulleqs, symstates)
  GetSymbolicFactorizedJacobian(x_eqs, y_eqs, x, y)
end

function Substitute(syms::Union{Vector{Num},Array{Num}}, subs_args::Matrix{Pair{Num,Float64}})
  return Symbolics.value.(substitute.(syms, (subs_args,)))
end
function Substitute(syms::Union{Vector{Num},Array{Num}}, subs_args::Vector{Pair{Num,Float64}})
  return Symbolics.value.(substitute.(syms, (subs_args,)))
end
function Substitute(syms::Union{Vector{Num},Array{Num}}, subs_args::Dict{SymbolicUtils.Symbolic{Real}, Float64})
  return Symbolics.value.(substitute.(syms, (subs_args,)))
end
function Substitute(syms::Union{Vector{Num},Array{Num}}, subs_args::Vector{Pair{SymbolicUtils.Symbolic{Real}, Float64}}) #SymbolicUtils.Symbolic{Real} ::Array{Pair{Num,Float64},1}
  return Symbolics.value.(substitute.(syms, (subs_args,)))
end
function Substitute(syms::Matrix{Num}, subs_args::Dict{Any, Any}) #SymbolicUtils.Symbolic{Real} ::Array{Pair{Num,Float64},1}
  return Symbolics.value.(substitute.(syms, (subs_args,)))
end

function GetMTKSystem(pg::PowerGrid, time_interval::Tuple{Float64,Float64}, p)
  U,δ,ic0,cu = PowerFlowClassic(pg,iwamoto = false)
  pg_new, ic =  InitializeInternalDynamics(pg,ic0)
  prob = ODEProblem(rhs(pg_new), ic, time_interval, p)
  new_f = ODEFunction(
    prob.f.f,
    syms = prob.f.syms,
    mass_matrix = Int.(prob.f.mass_matrix),
  )
  ODEProb = ODEProblem(new_f, ic, time_interval, p)
  return modelingtoolkitize(ODEProb)
end

function InitTrajectorySensitivity(
  mtsys::ODESystem,
  ic::Array{Float64,1},
  p::Array{Float64,1},
  sensis_u0_pd::Array{Any,1}, #array of indices
  sensis_p_pd::Array{Int64,1}, #array of indices
)
  fulleqs = equations(mtsys)
  sym_states = states(mtsys)
  sym_params = parameters(mtsys)
  eqs, aeqs, D_states, A_states = GetSymbolicEquationsAndStates(fulleqs, sym_states)
  sensis_u0 = sym_states[sensis_u0_pd]
  #sensis_p_pd is here a list with indices of the parameters p
  sensis_p = sym_params[sensis_p_pd]

  #dict from states and parameters with their starting values
  symu0 = sym_states .=> ic
  symp = sym_params .=> p

  Fx, Fy, Gx, Gy = GetSymbolicFactorizedJacobian(eqs, aeqs, D_states, A_states)

  len_sens = size(sensis_u0)[1] + size(sensis_p)[1]
  Diff_u0 = Differential.(sensis_u0)
  Diff_p = Differential.(sensis_p)
  Fp = Array{Num}(undef, size(eqs)[1], len_sens)
  Gp = Array{Num}(undef, size(aeqs)[1], len_sens)

  Fp[:, 1:size(Diff_u0)[1]] .= Num(0)
  Gp[:, 1:size(Diff_u0)[1]] .= Num(0)

  for (ind, val) in enumerate(Diff_p)
    Fp[:, ind+size(Diff_u0)[1]] =
      Num.(expand_derivatives.(map(val, my_rhs.(eqs))))
    Gp[:, ind+size(Diff_u0)[1]] =
      Num.(expand_derivatives.(map(val, my_rhs.(aeqs))))
  end

  @parameters Δt
  @parameters xx0[1:size(D_states)[1], 1:len_sens] #xx0 are the sensitivities regarding differential states
  @parameters yx0[1:size(A_states)[1], 1:len_sens] #yx0 are the sensitivities regarding algebraic states
  xx0 = Symbolics.scalarize(xx0)
  yx0 = Symbolics.scalarize(yx0)
  M,N,O = TrajectorySensitivityMatrices([Fx, Fy, Gx, Gy],Fp,Gp,xx0,yx0,aeqs,A_states,len_sens)

  #Initialisierung: xx0 enthält die Sensis für x0 und p bezüglich Differentialzustände
  xx0_k = xx0 .=> 0.0
  xx0_f = zeros(size(xx0)[1], len_sens)
  ind = setdiff(indexin(sensis_u0, D_states),[nothing])
  for i = 1:length(ind)
    xx0_k[i, ind[i]] = xx0_k[i, ind[i]][1] => 1.0
    xx0_f[i, ind[i]] = 1.0
  end
  # Bei den Sensis für y werden zuerst die dy/x0 Sensi initialisiert
  Gy_float = Substitute(Gy, [symu0; symp])
  # for increasing calculation of inv(Gy)
  yx0_t0 = -inv(Gy_float) * (Gx * xx0_f[:, 1:size(sensis_u0)[1]])
  yp_t0 =
    -inv(Gy_float) *
    (Gp * vcat(zeros(size(sensis_u0)[1], size(sensis_p)[1]), I))
  yx0_k = yx0 .=> Substitute([yx0_t0 yp_t0], [symu0; symp])

  eqs = Num.(my_rhs.(eqs))
  aeqs = Num.(my_rhs.(aeqs))
  return xx0_k,yx0_k,sym_states,sym_params,A_states,D_states,M,N,O,symp,Δt,len_sens, eqs,aeqs,(Fx,Fy,Gx,Gy)
end

function TrajectorySensitivityMatrices(
  J::Vector{Matrix{Num}},
  Fp::Matrix{Num},
  Gp::Matrix{Num},
  xx0::Matrix{Num},
  yx0::Matrix{Num},
  aeqs::Vector{Equation},
  A_states::Vector{Term{Real,Base.ImmutableDict{DataType,Any}}},
  len_sens::Int64,
)
    @parameters Δt
    Fx,Fy,Gx,Gy = J
    M = [
      Δt/2*Fx-I Δt/2*Fy
      Gx Gy
    ]
    N =
      isempty(aeqs) ? -xx0 - Δt / 2 * (Fx * xx0 + Fy * yx0 + Fp) :
      [
        -xx0 - Δt / 2 * (Fx * xx0 + Fy * yx0 + Fp)
        zeros(size(A_states)[1], len_sens)
      ]
    O = [
      -Δt / 2 * Fp
      -Gp
    ]
   return M,N,O
end

function ContinuousSensitivity(
  sol::ODESolution,
  xx0_k::Matrix{Pair{Num,Float64}},
  yx0_k::Matrix{Pair{Num,Float64}},
  sym_states::Vector{Term{Real,Base.ImmutableDict{DataType,Any}}},
  A_states::Vector{Term{Real,Base.ImmutableDict{DataType,Any}}},
  D_states::Vector{Term{Real,Base.ImmutableDict{DataType,Any}}},
  M::Matrix{Num},
  N::Matrix{Num},
  O::Matrix{Num},
  symp::Vector{Pair{Sym{Real, Base.ImmutableDict{DataType, Any}}, Float64}},
  Δt::Num,
  len_sens::Int64,
)
  sensi = Vector{Array{Float64}}(undef, len_sens)
  for i = 1:length(sensi)
    sensi[i] = Array{Float64}(
      undef,
      size(D_states)[1] + size(A_states)[1],
      size(sol)[2] - 1,
    )
  end
  xx0 = [i[1] for i in xx0_k]
  yx0 = [i[1] for i in yx0_k]
  ind_y = setdiff(indexin(A_states, sym_states), [nothing])
  ind_x = setdiff(indexin(D_states, sym_states), [nothing])
  for i = 1:size(sol)[2]-1
    dt = sol.t[i+1] - sol.t[i]

    uk = sym_states .=> sol.u[i]
    uk1 = sym_states .=> sol.u[i+1]

    Mfloat = Float64.(Substitute(M, [uk1; symp; Δt => dt]))
    Nfloat = Float64.(Substitute(N, [uk; symp; vec(xx0_k); vec(yx0_k); Δt => dt]))
    Ofloat = Float64.(Substitute(O, [uk1; symp; Δt => dt]))
    res = inv(Mfloat) * (Nfloat + Ofloat)

    for j = 1:length(sensi)
      sensi[j][ind_x, i] = res[1:size(D_states)[1], j]
      sensi[j][ind_y, i] = res[size(D_states)[1]+1:end, j]
    end

    xx0_k = xx0 .=> res[1:size(D_states)[1], :]
    yx0_k = yx0 .=> res[size(D_states)[1]+1:end, :]
  end
  return sensi, xx0_k, yx0_k
end

function CalcTriggerAndStateResetJacobians(mtk::ODESystem,s::Vector{Equation},h::Vector{Matrix{Equation}})
    eqs, aeqs, x, y = GetSymbolicEquationsAndStates(mtk)
    hx = Array{Array{Num}}(undef,length(h),1)
    hy = similar(hx)
    sx = Array{Array{Num}}(undef,length(s),1)
    sy = similar(sx)
    for i=1:length(h)
        hx[i] = GetJacobian(h[i],x)
        hy[i] = GetJacobian(h[i],y)
    end
    for i=1:length(s)
        sx[i] = GetJacobian([s[i]],x)
        sy[i] = GetJacobian([s[i]],y)
    end
    return hx,hy,sx,sy
end

function CalcSensitivityAfterJump(
    sym_states::Vector{Term{Real,Base.ImmutableDict{DataType,Any}}},
    sym_params::Vector{Sym{Real,Base.ImmutableDict{DataType,Any}}},
    xx0_pre::VecOrMat{Float64},
    yx0_pre::VecOrMat{Float64},
    x0_pre::VecOrMat{Float64},
    x0_post::VecOrMat{Float64},
    p_pre::VecOrMat{Float64},
    p_post::VecOrMat{Float64},
    f_pre::VecOrMat{Num},
    f_post::VecOrMat{Num},
    g_pre::VecOrMat{Num},
    g_post::VecOrMat{Num},
    J_pre::NTuple{4,Matrix{Num}},
    J_post::NTuple{4,Matrix{Num}},
    hx::Matrix{Num},
    hy::Matrix{Num},
    sx::Matrix{Num},
    sy::Matrix{Num},
)
    fx_pre, fy_pre, gx_pre, gy_pre = J_pre
    fx_post, fy_post, gx_post, gy_post = J_post

    subs_pre = [sym_states .=> x0_pre; sym_params .=> p_pre]
    subs_post = [sym_states .=> x0_post; sym_params .=> p_post]

    f_pre_float = Float64.(Substitute(f_pre, subs_pre))
    f_post_float = Float64.(Substitute(f_post, subs_post))

    fx_pre_float = Float64.(Substitute(fx_pre, subs_pre))
    fy_pre_float = Float64.(Substitute(fy_pre, subs_pre))
    gx_pre_float = Float64.(Substitute(gx_pre, subs_pre))
    gy_pre_float = Float64.(Substitute(gy_pre, subs_pre))

    gx_post_float = Float64.(Substitute(gx_post, subs_post))
    gy_post_float = Float64.(Substitute(gy_post, subs_post))

    hx_pre = Float64.(Substitute(hx, subs_pre))
    hy_pre = Float64.(Substitute(hy, subs_pre))
    sx_pre = Float64.(Substitute(sx, subs_pre))
    sy_pre = Float64.(Substitute(sy, subs_pre))

    gygx = inv(gy_pre_float) * gx_pre_float

    hx_star = hx_pre- hy_pre * gygx

    s_star = sx_pre - sy_pre * gygx

    τx0 = s_star * xx0_pre
    tmp = s_star * f_pre_float
    if sum(tmp) != 0.0
        τx0 = s_star * xx0_pre ./ (tmp)
    else
        τx0 = 0.0.*τx0
    end

    xx0_post = hx_star  * xx0_pre - (f_post_float - hx_star * f_pre_float) * τx0
    yx0_post = -inv(gy_post_float) * gx_post_float * xx0_post

    #return xx0_pre, yx0_pre
    return xx0_post, yx0_post
end

function CalcHybridTrajectorySensitivity(
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
    xx0_k, yx0_k, sym_states,sym_params, A_states, D_states, M, N, O, symp, Δt,len_sens, f, g, J =
        InitTrajectorySensitivity(mtk0, ic, p_pre, u0_sensi, p_sensi)
    xx0 = [i[1] for i in xx0_k]
    yx0 = [i[1] for i in yx0_k]
    f_all,g_all,J_all,M_all,N_all,O_all = GetEqsJacobianSensMatrices(mtk,xx0,yx0,u0_sensi,p_sensi)
    Fx_all, Fy_all, Gx_all, Gy_all = J_all
    hx,hy,sx,sy = CalcTriggerAndStateResetJacobians(mtk0,s,h)
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

    for i = 1:length(ind_sol)-1
        sol_part = sol[ind_sol[i]:ind_sol[i+1]]
        sensi_part,xx0_k,yx0_k = ContinuousSensitivity(
                                    sol_part,
                                    xx0_k,
                                    yx0_k,
                                    sym_states, #at first, assumed to be constant
                                    A_states, #at first, assumed to be constant
                                    D_states, #at first, assumed to be constant
                                    M, #can change after event
                                    N, #can change after event
                                    O, #can change after event
                                    symp,
                                    Δt,
                                    len_sens,
                                )
        for j = 1:length(sensi_part)
            sensis[j][:,ind_sol[i]:ind_sol[i+1]-1] = sensi_part[j]
        end
        if ind_sol[i+1] ≠ ind_sol[end]
            display("Event at t = $(evr[i,1])")

            xx0_pre = [i[2] for i in xx0_k]
            yx0_pre = [i[2] for i in yx0_k]
            x0_pre = sol.u[ind_sol[i+1]]    # ind_sol[i+1]-1 is before jump
            x0_post = sol.u[ind_sol[i+1]+1]    # ind_sol[i+1] is after jump
            p_post = evr[i,3:end-2]
            hx_tmp = hx[Int(evr[i,end])]
            hy_tmp = hy[Int(evr[i,end])]
            sx_tmp = sx[Int(evr[i,end-1])]
            sy_tmp = sy[Int(evr[i,end-1])]

            Fx_post = Fx_all[Int(evr[i,2])]
            Fy_post = Fy_all[Int(evr[i,2])]
            Gx_post = Gx_all[Int(evr[i,2])]
            Gy_post = Gy_all[Int(evr[i,2])]
            J_post = Fx_post, Fy_post, Gx_post, Gy_post

            f_post = f_all[Int(evr[i,2])]
            g_post = g_all[Int(evr[i,2])]

            xx0_post, yx0_post = CalcSensitivityAfterJump(
                                    sym_states,
                                    sym_params,
                                    xx0_pre,
                                    yx0_pre,
                                    x0_pre,
                                    x0_post,
                                    p_pre,
                                    p_post,
                                    f_pre,
                                    f_post,
                                    g_pre,
                                    g_post,
                                    J_pre,
                                    J_post,
                                    hx_tmp,
                                    hy_tmp,
                                    sx_tmp,
                                    sy_tmp,
                                )
        else
            return sensis
        end
        xx0_k = xx0 .=> xx0_post
        yx0_k = yx0 .=> yx0_post
        symp = sym_params .=> p_post
        p_pre = p_post

        Mc = copy(M)
        M,N,O = M_all[Int(evr[i,2])], N_all[Int(evr[i,2])], O_all[Int(evr[i,2])]

        Fx_pre = deepcopy(Fx_post)
        Fy_pre = deepcopy(Fy_post)
        Gx_pre = deepcopy(Gx_post)
        Gy_pre = deepcopy(Gy_post)
        J_pre = Fx_pre, Fy_pre, Gx_pre, Gy_pre

        f_pre = deepcopy(f_post)
        g_pre = deepcopy(g_post)
    end
end

function GetEqsJacobianSensMatrices(mtk::Vector{ODESystem},xx0::Matrix{Num},yx0::Matrix{Num},u0_sens::Vector{Union{Int64,Any}},p_sens::Vector{Int64})
    f = Vector{Array{Num,1}}(undef,length(mtk))
    g = similar(f)

    Fx = Vector{Array{Num,2}}(undef,length(mtk))
    Fy = similar(Fx)
    Gx = similar(Fx)
    Gy = similar(Fx)
    Fp = similar(Fx)
    Gp = similar(Fx)

    M = similar(Fx)
    N = similar(Fx)
    O = similar(Fx)

    len_sens = size(u0_sens)[1] + size(p_sens)[1]
    sensis_u0 = states(mtk[1])[u0_sens]
    sensis_p = parameters(mtk[1])[p_sens]
    Diff_u0 = Differential.(sensis_u0)
    Diff_p = Differential.(sensis_p)

    for (ind,val) in enumerate(mtk)
        fulleqs = equations(val)
        symstates = states(val)
        f_tmp, g_tmp, x, y = GetSymbolicEquationsAndStates(fulleqs, symstates)
        Fx[ind],Fy[ind],Gx[ind],Gy[ind] = GetSymbolicFactorizedJacobian(f_tmp, g_tmp, x, y)

        Fp[ind] = Num(0).*zeros(length(x),len_sens) #init with zeros
        Gp[ind] = Num(0).*zeros(length(y),len_sens)
        for (ind2, val2) in enumerate(Diff_p)
          Fp[ind][:, ind2+size(Diff_u0)[1]] =
            Num.(expand_derivatives.(map(val2, my_rhs.(f_tmp))))
          Gp[ind][:, ind2+size(Diff_u0)[1]] =
            Num.(expand_derivatives.(map(val2, my_rhs.(g_tmp))))
        end

        M[ind],N[ind],O[ind] = TrajectorySensitivityMatrices([Fx[ind],Fy[ind],Gx[ind],Gy[ind]],Fp[ind],Gp[ind],xx0,yx0,g_tmp,y,len_sens)
        f[ind] = my_rhs.(f_tmp)
        g[ind] = my_rhs.(g_tmp)
    end

    return f,g,[Fx,Fy,Gx,Gy],M,N,O
end
