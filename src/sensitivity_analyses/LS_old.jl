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

function TimeDomainSensitivies(pg::PowerGrid,time_interval,ic,p,sensis_u0_pd,sensis_p_pd)
    prob = ODEProblem(rhs(pg),ic,time_interval,p)
    sol  = solve(prob,Rodas4())
    TimeDomainSensitivies(pg,time_interval,ic,p,sensis_u0_pd,sensis_p_pd,sol)
end

function TimeDomainSensitivies(pg::PowerGrid,time_interval,ic,p,sensis_u0_pd,sensis_p_pd,sol)
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
  state::Vector{Term{Real,Nothing}},
)
  aeqs = Vector{Equation}()
  eqs = Vector{Equation}()
  A_states = Vector{Term{Real,Nothing}}()
  D_states = Vector{Term{Real,Nothing}}()
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

function GetJacobian(
  eqs::Array{Equation},
  states::Array{Term{Real,Nothing},1},
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
  D_states::Array{Term{Real,Nothing},1},
  A_states::Array{Term{Real,Nothing},1},
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

function Substitute(syms::Array{Num}, subs_args) #SymbolicUtils.Symbolic{Real} ::Array{Pair{Num,Float64},1}
  return Symbolics.value.(substitute.(syms, (subs_args,)))
end


function GetMTKSystem(pg::PowerGrid, time_interval::Tuple{Float64,Float64}, p::Array{Float64,1})
  U,δ,ic0 = PowerFlowClassic(pg,iwamoto = false)
  Ykk = NodalAdmittanceMatrice(pg)
  Uc = U.*exp.(1im*δ/180*pi)
  I_c = Ykk*Uc
  pg_new, ic =  InitializeInternalDynamics(pg,I_c,ic0)
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
  sensis_u0_pd::Array{Int64,1},
  sensis_p_pd::Array{Int64,1},
)
  fulleqs = equations(mtsys)
  sym_states = states(mtsys)
  sym_params = parameters(mtsys)
  eqs, aeqs, D_states, A_states = GetSymbolicEquationsAndStates(fulleqs, sym_states)

  #it is assumed that state and rhs(powergrid).syms have the same order
  #sensis_u0 = state[indexin(sensis_u0_pd, rhs(pg).syms)]
  sensis_u0 = sym_states[sensis_u0_pd]
  #sensis_p_pd is here a list with indices of the parameters p
  sensis_p = sym_params[sensis_p_pd]

  #dict from states and parameters with their starting values
  symu0 = sym_states .=> ic
  symp = sym_params .=> p

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
  @parameters xx0[1:size(D_states)[1], 1:len_sens] #xx0 are the sensitivities regarding differential states
  @parameters yx0[1:size(A_states)[1], 1:len_sens] #yx0 are the sensitivities regarding algebraic states
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

  eqs = Num.(my_rhs.(eqs))
  aeqs = Num.(my_rhs.(aeqs))
  return xx0_k,yx0_k,sym_states,sym_params,A_states,D_states,M,N,O,symp,Δt,len_sens, eqs,aeqs,(Fx,Fy,Gx,Gy)
end

function ContinuousSensitivity(sol,xx0_k,yx0_k,sym_states,A_states,D_states,M,N,O,symp,Δt,len_sens)
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

    Mfloat = Substitute(M, [uk1; symp; Δt => dt])
    Nfloat = Substitute(N, [uk; symp; vec(xx0_k); vec(yx0_k); Δt => dt])
    Ofloat = Substitute(O, [uk1; symp; Δt => dt])
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

function CalcTriggerAndStateResetJacobians(mtk::ODESystem,s,h)
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
    sym_states,
    sym_params,
    xx0_pre,
    yx0_pre,
    x0_pre,
    x0_post,
    p_pre,
    p_post,
    f,
    g,
    J,
    hx,
    hy,
    sx,
    sy,
)
    fx, fy, gx, gy = J

    subs_pre = [sym_states .=> x0_pre; sym_params .=> p_pre]
    subs_post = [sym_states .=> x0_post; sym_params .=> p_post]

    f_pre = Substitute(f, subs_pre)
    f_post = Substitute(f, subs_post)

    fx_pre = Substitute(fx, subs_pre)
    fy_pre = Substitute(fy, subs_pre)
    gx_pre = Substitute(gx, subs_pre)
    gy_pre = Substitute(gy, subs_pre)

    gx_post = Substitute(gx, subs_post)
    gy_post = Substitute(gy, subs_post)

    hx_pre = Substitute(hx, subs_pre)
    hy_pre = Substitute(hy, subs_pre)
    sx_pre = Substitute(sx, subs_pre)
    sy_pre = Substitute(sy, subs_pre)

    gygx = inv(gy_pre) * gx_pre

    hx_star = hx_pre - hy_pre * gygx

    s_star = sx_pre - sy_pre * gygx

    τx0 = s_star * xx0_pre
    tmp = s_star * f_pre
    if sum(tmp) != 0.0
        τx0 = s_star * xx0_pre ./ (tmp)
    else
        τx0 = 0.0.*τx0
    end

    xx0_post = hx_star  * xx0_pre - (f_post - hx_star * f_pre) * τx0
    yx0_post = -inv(gy_post) * gx_post * xx0_post

    return xx0_post, yx0_post
end

function CalcHybridTrajectorySensitivity(mtk,sol,p_pre,evr,s,h,u0_sensi,p_sensi)
    ic = sol.prob.u0
    xx0_k, yx0_k, sym_states,sym_params, A_states, D_states, M, N, O, symp, Δt,len_sens, f, g, J =
        InitTrajectorySensitivity(mtk, ic, p_pre, u0_sensi, p_sensi)
    display(xx0_k)
    display(yx0_k)
    display(sym_states)
    display(sym_params)
    display(A_states)
    display(D_states)
    display(symp)
    display(len_sens)
    display(Δt)
    xx0 = [i[1] for i in xx0_k]
    yx0 = [i[1] for i in yx0_k]
    #fx,fy,gx,gy = J
    hx,hy,sx,sy = CalcTriggerAndStateResetJacobians(mtk,s,h)
    sensis = Vector{Array{Float64}}(undef, len_sens)
    for i = 1:length(sensis)
      sensis[i] = Array{Float64}(
        undef,
        size(D_states)[1] + size(A_states)[1],
        size(sol)[2] - 1,
      )
    end
    ind_sol = vcat(1,setdiff(indexin(evr[:,1],sol.t).+1,[nothing]),length(sol.t))

    #ind_sol = [1]
    #for i in evr[:,1] # DifferentialEquations.jl has multiple time points
    #    ind_sol = vcat(ind_sol,findall(x->x==i,sol.t)[end])
    #end
    #ind_sol = vcat(ind_sol,length(sol.t))

    @progress for i = 1:length(ind_sol)-1
        sol_part = sol[ind_sol[i]:ind_sol[i+1]]
        sensi_part,xx0_k,yx0_k = ContinuousSensitivity(
                                    sol_part,
                                    xx0_k,
                                    yx0_k,
                                    sym_states,
                                    A_states,
                                    D_states,
                                    M,
                                    N,
                                    O,
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

            xx0_post, yx0_post = CalcSensitivityAfterJump(
                                    sym_states,
                                    sym_params,
                                    xx0_pre,
                                    yx0_pre,
                                    x0_pre,
                                    x0_post,
                                    p_pre,
                                    p_post,
                                    f,
                                    g,
                                    J,
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
    end
end