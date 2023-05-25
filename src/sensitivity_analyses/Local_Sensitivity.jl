using PowerDynamics
using ModelingToolkit
using LinearAlgebra


# Annahmen für die Trajektorien-Sensitivität
# 1. das Netzmodell muss als DAIS aufgebaut sein:
#   1.1 diskrete Zustände als Dummy-DGL
#   1.2 Differentialgleichungen bleiben gleich (können aber diskrete Zustände haben)
#   1.3 durch die Reset-Funktionen ändern sich nur diskrete Zustände
#   1.4 algebraische Gleichungen können sich 'beliebig' ändern
#   1.5 die Anzahl der Systemzuständen bleibt jedoch immer gleich (zur Not null/NaN )
# 2. Parameter ändern sich NICHT zur Laufzeit
# 3. Aufbau des Event-Recorders (evr)
#   3.1 [t_i mtk_i s_i h_i]
#   3.1.1 t_i = Zeit, mtk_i = 'aktives Netz' bei Änd. der alg. Gl., s_i = switch cond., h_i = reset function
#   3.2 Initialisierung: evr = Array{Float64}(undef,0,4)
#   3.3 Eintrag machen (Bsp.) evr = vcat(evr, [integrator.t 1 1 1])
#   3.2 der erste Eintrag ist immer evr[1] = [t_0 1 0 0]
#   3.2.1 also bei t0 gehts los
# 4. die Simulation endet nicht mit einem relevanten Event
# 5. Bei den Events von DifferentialEquations.jl ist es so, dass es im Lösungsobjekt später
# zwei gleiche Zeitpunkte gibt. Also Δt ist dazwischen null. Der erste Zeitpunkt davon wird für τ- und 
# der zweite dann für τ+ genommen. Es sollten nicht mehr als zwei gleiche Zeitpunkte auftreten.
# Das wird bisher jedoch nicht überprüft.

function plot_sensi(time,sensis)
   p = plot(layout = (size(sensis)[1],1))
   for (index, value) in enumerate(sensis)
      plot!(time,value',subplot=index)
   end
   return p
end

my_rhs(x::Equation) = x.rhs
my_rhs(x::Num) = x
my_lhs(x::Equation) = x.lhs
my_lhs(x::Num) = x

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

function GetFactorisedSymbolicStates(mtsys::ODESystem)
  fulleqs = equations(mtsys)
  state = states(mtsys)
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
  eqs::Array{Num},
  states::Vector{Num},
)
    Fx = Array{Num}(undef, size(eqs)[1], size(states)[1])
    Diff_states = Differential.(states)
    for (ind, val) in enumerate(Diff_states)
      Fx[:, ind] = Num.(expand_derivatives.(map(val, eqs)))
    end
    return Fx
end

function GetJacobian(
  eqs::Array{Num},
  states::Vector{Term{Real, Base.ImmutableDict{DataType, Any}}},
)
    Fx = Array{Num}(undef, size(eqs)[1], size(states)[1])
    Diff_states = Differential.(states)
    for (ind, val) in enumerate(Diff_states)
      Fx[:, ind] = Num.(expand_derivatives.(map(val, eqs)))
    end
    return Fx
end

function GetJacobian(
  eqs::Array{Equation},
  states::Vector{SymbolicUtils.Symbolic{Real}},
)
    Fx = Array{Num}(undef, size(eqs)[1], size(states)[1])
    Diff_states = Differential.(states)
    for (ind, val) in enumerate(Diff_states)
      Fx[:, ind] = Num.(expand_derivatives.(map(val, my_rhs.(eqs))))
    end
    return Fx
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
  D_states::Vector{SymbolicUtils.Symbolic{Real}},
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
function Substitute(syms::Matrix{Num}, subs_args::Vector{Pair{_A, Float64} where _A})
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

function InitTrajectorySensitivity(mtk::Vector{ODESystem},sol::ODESolution)
  sym_states = states(mtk[1]);
  sym_params = parameters(mtk[1]);

  D0_states, A_states = GetFactorisedSymbolicStates(mtk[1]);
  D_states = Num.(vcat(D0_states,sym_params));
  A_states = Num.(A_states);

  D0_indices= Int64.(setdiff(indexin(D0_states, sym_states), [nothing]));
  A_indices = Int64.(setdiff(indexin(A_states, sym_states), [nothing]));

  #Hier anpassen, welche Zustände Parameter berücksichtigt werden soll
  sensis_x = D_states; #anpassen bei variabler Größe
  sensis_y = A_states; #anpassen bei variabler Größe?

  #dict from states and parameters with their starting values
  sym_x = sensis_x .=> vcat(sol.prob.u0[D0_indices],sol.prob.p);
  sym_y = sensis_y .=> sol.prob.u0[A_indices]; 
  len_sens = length(sym_x) #+ length(sym_y) #den Einfluss der algebraischen Variablen steckt indirekt in den Gleichungen!!

  @parameters Δt t
  @parameters xx0[1:length(D_states), 1:len_sens] #xx0 are the symbolic sensitivities regarding differential states
  @parameters yx0[1:length(A_states), 1:len_sens] #yx0 are the symbolic sensitivities regarding algebraic states
  xx0 = Symbolics.scalarize(xx0)
  yx0 = Symbolics.scalarize(yx0)

  # die Sensitivität eines Parameters/Zustand muss immer auf alle Zustände/Parameter berechnet werden: x_x1, p_x1, y_x1 
  # somit müssen die Ableitungen für ALLE Zustände berechnet werden!
  matrices =  GetEqsJacobianSensMatrices(mtk,xx0,yx0) #anpassen wenn nicht alle Variablen gerechnet werden sollen
  F_all,G_all,J_all,M,N = matrices;
  Fx_all,Fy_all,Gx_all,Gy_all = J_all;

  #Initialisierung der differntialen symbolischen Sensitivitäten: xx0 enthält die Sensis für x, q und p 
  xx0_k = xx0 .=> 0.0
  for i in diagind(xx0)
      xx0_k[i] = xx0[i] => 1.0
  end

  #Initialisierung der algebraischen symbolischen Sensitivitäten
  Gy_float = Float64.(Substitute(Gy_all[1], [sym_x; sym_y; t => sol.t[1]]))
  Gx_float = Float64.(Substitute(Gx_all[1], [sym_x; sym_y; t => sol.t[1]]))
  yx0f = -inv(Gy_float) * Gx_float
  yx0_k = yx0 .=> yx0f

  #Objekterstelleung der Trajektorien Sensitivität
  sensis = Vector{Array{Float64}}(undef, len_sens)
  for i = 1:length(sensis)
    sensis[i] = zeros(Float64,size(D_states)[1] + size(A_states)[1],size(sol)[2])
  end
  # Initialisierung der Trajektorien Sensitivität
  for i= eachindex(D_states)
    sensis[i][i,1] = 1.0
  end
  for i= eachindex(A_states) .+ length(D_states)
    sensis[i-length(D_states)][i,1] = yx0f[i-length(D_states)]
  end
  return xx0_k,yx0_k,A_states,D_states, A_indices, D0_indices, matrices,Δt,len_sens, sensis
end

function TrajectorySensitivityMatrices(J::Vector{Matrix{Num}},xx0::Matrix{Num},yx0::Matrix{Num})
    @parameters Δt
    Fx,Fy,Gx,Gy = J
    M = [
      Δt/2*Fx-I Δt/2*Fy
      Gx Gy
    ]
    N =
      isempty(yx0) ? -xx0 - Δt / 2 * (Fx * xx0 + Fy * yx0) :
      [
        -xx0 - Δt / 2 * (Fx * xx0 + Fy * yx0)
        zeros(size(yx0)[1], size(yx0)[2])
      ]
   return M,N
end

function ContinuousSensitivity(
  sol::ODESolution,
  xx0_k::Matrix{Pair{Num,Float64}},
  yx0_k::Matrix{Pair{Num,Float64}},
  A_states::Vector{Num},
  D_states::Vector{Num},
  A_indices::Vector{Int64},
  D0_indices::Vector{Int64},
  M::Matrix{Num},
  N::Matrix{Num},
  Δt::Num,
  len_sens::Int64,
)
  #lokaler Speicher für die berechneten Trajektorien Senstivitäten
  sensi = Vector{Array{Float64}}(undef, len_sens)
  # die sensis sind hier eins kürzer, da ab k+1 bis zum Ende gerechnet wird!
  # t_k ist der Zeitpunkt zuvor (also entweder Zeitpuntk des Sprungs oder Startzeit)
  for i = 1:length(sensi)
    sensi[i] = 0.0.*Array{Float64}(
      undef,
      size(D_states)[1] + size(A_states)[1],
      size(sol)[2]-1,
    )
  end
  xx0 = [i[1] for i in xx0_k]
  yx0 = [i[1] for i in yx0_k]
  #ind_y = Int64.(setdiff(indexin(A_states, sym_states), [nothing]))
  #ind_x = Int64.(setdiff(indexin(D_states, sym_states), [nothing]))
  @parameters t

  for i = 2:size(sol)[2]
    dt = sol.t[i] - sol.t[i-1]

    xk  = D_states .=> vcat(sol[D0_indices,i-1],sol.prob.p)
    xk1 = D_states .=> vcat(sol[D0_indices,i],sol.prob.p)

    yk  = A_states .=> sol[A_indices,i-1]
    yk1 = A_states .=> sol[A_indices,i]

    Mfloat = Float64.(Substitute(M, [xk1; yk1; Δt => dt; t => sol.t[i]]))
    Nfloat = Float64.(Substitute(N, [xk; yk; vec(xx0_k); vec(yx0_k); Δt => dt; t => sol.t[i-1]]))
    res = inv(Mfloat) * Nfloat 

    for j = 1:length(sensi)
      sensi[j][1:length(D_states), i-1] = res[1:length(D_states), j]
      sensi[j][length(D_states)+1:end, i-1] = res[length(D_states)+1:end, j]
    end
    # Lesezeichen
    xx0_k = xx0 .=> res[1:size(D_states)[1],:]
    yx0_k = yx0 .=> res[size(D_states)[1]+1:end,:]
  end
  return sensi, xx0_k, yx0_k
end

function CalcTriggerAndStateResetJacobians(s::Vector{Num},h::Vector{Vector{Num}},D_states,A_states)
  hx = Array{Array{Num}}(undef,length(h),1)
  hy = similar(hx)
  sx = Array{Array{Num}}(undef,length(s),1)
  sy = similar(sx)

  for i=eachindex(hx)
    hx[i] = GetJacobian(h[i],D_states)
    hy[i] = GetJacobian(h[i],A_states)
  end 
  for i=eachindex(s)
    sx[i] = GetJacobian([s[i]],D_states)
    sy[i] = GetJacobian([s[i]],A_states)
  end  
  return hx,hy,sx,sy
end

function CalcSensitivityAfterJump(
    xy_k::Vector{Pair{Num, Float64}},
    xy_k1::Vector{Pair{Num, Float64}},
    xx0_k_float::VecOrMat{Float64},
    f_pre::VecOrMat{Num},
    f_post::VecOrMat{Num},
    J_pre::NTuple{4,Matrix{Num}} ,
    J_post::NTuple{4,Matrix{Num}},
    hx::Matrix{Num},
    hy::Matrix{Num},
    sx::Matrix{Num},
    sy::Matrix{Num},
)
    fx_pre, fy_pre, gx_pre, gy_pre = J_pre
    fx_post, fy_post, gx_post, gy_post = J_post

    f_pre_float = Float64.(Substitute(f_pre, xy_k))
    f_post_float = Float64.(Substitute(f_post, xy_k1))

    gx_pre_float = Float64.(Substitute(gx_pre, xy_k))
    gy_pre_float = Float64.(Substitute(gy_pre, xy_k))

    gx_post_float = Float64.(Substitute(gx_post, xy_k1))
    gy_post_float = Float64.(Substitute(gy_post, xy_k1))

    hx_pre = Float64.(Substitute(hx, xy_k))
    hy_pre = Float64.(Substitute(hy, xy_k))
    sx_pre = Float64.(Substitute(sx, xy_k))
    sy_pre = Float64.(Substitute(sy, xy_k))

    gygx = inv(gy_pre_float) * gx_pre_float

    hx_star = hx_pre- hy_pre * gygx
    s_star = sx_pre - sy_pre * gygx

    τx0_nom = s_star * xx0_k_float
    τx0_denom = s_star * f_pre_float

    if τx0_denom[1] == 0.0
      error("Trajectory is not hitting the switching surface transversally.")
    end
    τx0 = -τx0_nom ./ τx0_denom

    xx0_post = hx_star  * xx0_k_float - (f_post_float - hx_star * f_pre_float) * τx0
    yx0_post = -inv(gy_post_float) * gx_post_float * xx0_post

    return xx0_post, yx0_post
end

function GetEqsJacobianSensMatrices(mtk::Vector{ODESystem},xx0::Matrix{Num},yx0::Matrix{Num})
  f = Vector{Array{Num,1}}(undef,length(mtk))
  g = similar(f)

  Fx = Vector{Array{Num,2}}(undef,length(mtk))
  Fy = similar(Fx)
  Gx = similar(Fx)
  Gy = similar(Fx)

  M = similar(Fx)
  N = similar(Fx)

  for (ind,val) in enumerate(mtk)
      fulleqs = equations(val)
      symstates = states(val)
      sym_params = parameters(val)
      eqs, aeqs, x, y = GetSymbolicEquationsAndStates(fulleqs, symstates)

      x = vcat(x,sym_params);
      @variables t
      D= Differential(t);
      for i=eachindex(sym_params) #hinzufügen der Dummy-Gleichungen für die Parameter
          eqs = vcat(eqs,D(sym_params[i])~0.0)
      end
      Fx[ind],Fy[ind],Gx[ind],Gy[ind] = GetSymbolicFactorizedJacobian(eqs, aeqs, x, y)

      M[ind],N[ind] = TrajectorySensitivityMatrices([Fx[ind],Fy[ind],Gx[ind],Gy[ind]],xx0,yx0)
      f[ind] = my_rhs.(eqs)
      g[ind] = my_rhs.(aeqs)
  end

  return f,g,[Fx,Fy,Gx,Gy],M,N
end

function CalcContinuousSensitivity(mtk::Vector{ODESystem},sol::ODESolution)
    xx0_k,yx0_k,A_states,D_states, A_indices, D0_indices,matrices,Δt,len_sens, sensis = InitTrajectorySensitivity(mtk,sol);
    f_all, g_all, J_all, M_all, N_all = matrices

    sensi_part,xx0_k,yx0_k = ContinuousSensitivity(
      sol,
      xx0_k,
      yx0_k,
      A_states, 
      D_states, 
      A_indices,
      D0_indices,
      M_all[1], 
      N_all[1], 
      Δt,
      len_sens)

    return sensi_part
end

function CalcHybridTrajectorySensitivity(
  mtk::Vector{ODESystem},
  sol::ODESolution,
  evr::Matrix{Float64},
  s::Vector{Num},
  h::Vector{Vector{Num}},
)
    #u0_sensi::Vector{Union{Int64,Any}},
    #p_sensi::Vector{Int64},

    xx0_k,yx0_k,A_states,D_states, A_indices, D0_indices,matrices,Δt,len_sens, sensis = InitTrajectorySensitivity(mtk,sol);
    f_all, g_all, J_all, M_all, N_all = matrices
    Fx_all, Fy_all, Gx_all, Gy_all = J_all
    xx0 = [i[1] for i in xx0_k]
    yx0 = [i[1] for i in yx0_k]

    Fx_pre = Fx_all[1]
    Fy_pre = Fy_all[1]
    Gx_pre = Gx_all[1]
    Gy_pre = Gy_all[1]
    J_pre = Fx_pre, Fy_pre, Gx_pre, Gy_pre

    f_pre = f_all[1]
    g_pre = g_all[1]

    Mcont = copy(M_all[1])
    Ncont = copy(N_all[1])

    #Bestimmung der partiellen Ableitungen für switchinig conditions and reset functions
    hx,hy,sx,sy = CalcTriggerAndStateResetJacobians(s,h,D_states,A_states)

    # Bestimmung der Zeitpunkte der Events im Loesungsobjekt + Anfang & Ende 
    #ind_sol = vcat(2,setdiff(indexin(evr[:,1],sol.t),[nothing]),length(sol.t))

    # ind1 = Anfang der Indizes für die kontinuierlichen Sensis
    ind1 = vcat(1,setdiff(indexin(evr[:,1],sol.t),[nothing]).+1)
    # ind2 = Ende der Indizes für die kontinuierlichen Sensis
    ind2 = vcat(setdiff(indexin(evr[:,1],sol.t),[nothing]),length(sol.t))
    # zwischen den Anfang und End Indizes liegt das Event mit Zustandsänderung
    ind_sol = Int.(hcat(ind1,ind2))

    @parameters t

    for i = 1:size(ind_sol)[1]
        sol_part = sol[ind_sol[i,1]:ind_sol[i,2]]
        sensi_part,xx0_k,yx0_k = ContinuousSensitivity(
                                                      sol_part,
                                                      xx0_k,
                                                      yx0_k,
                                                      A_states, 
                                                      D_states, 
                                                      A_indices,
                                                      D0_indices,
                                                      Mcont, 
                                                      Ncont, 
                                                      Δt,
                                                      len_sens)
        for j = eachindex(sensi_part)
            sensis[j][:,ind_sol[i,1]+1:ind_sol[i,2]] = sensi_part[j]
        end
        if ind_sol[i,2] ≠ ind_sol[end,2]
            display("Event at t = $(evr[i,1])")

            xx0_k_float = [i[2] for i in xx0_k]

            xk = D_states .=> vcat(sol[D0_indices,ind_sol[i,2]],sol.prob.p)    #hier τ-
            xk1 = D_states .=> vcat(sol[D0_indices,ind_sol[i,2]+1],sol.prob.p) #hier τ+
        
            yk  = A_states .=> sol[A_indices,ind_sol[i,2]]   #hier τ-
            yk1 = A_states .=> sol[A_indices,ind_sol[i,2]+1] #hier τ+

            #p_post = evr[i,3:end-2]
            hx_tmp = hx[Int(evr[i,4])]
            hy_tmp = hy[Int(evr[i,4])]
            sx_tmp = sx[Int(evr[i,3])]
            sy_tmp = sy[Int(evr[i,3])]

            Fx_post = Fx_all[Int(evr[i,2])]
            Fy_post = Fy_all[Int(evr[i,2])]
            Gx_post = Gx_all[Int(evr[i,2])]
            Gy_post = Gy_all[Int(evr[i,2])]
            J_post = Fx_post, Fy_post, Gx_post, Gy_post

            f_post = f_all[Int(evr[i,2])]
            g_post = g_all[Int(evr[i,2])]

            xx0_k1_float, yx0_k1_float = CalcSensitivityAfterJump(
                                    [xk; yk;  t => sol.t[ind_sol[i,2]]],
                                    [xk1;yk1; t => sol.t[ind_sol[i,2]+1]],
                                    xx0_k_float,
                                    f_pre,
                                    f_post,
                                    J_pre,
                                    J_post,
                                    hx_tmp,
                                    hy_tmp,
                                    sx_tmp,
                                    sy_tmp,
                                )
            # speichern der Sensis nach dem Sprung
            for k in size(xx0_k1_float)[2]
              sensis[k][:,ind_sol[i,2]+1] = vcat(xx0_k1_float[:,k],yx0_k1_float[:,k])
            end
        else
            return sensis
        end

        xx0_k = xx0 .=> xx0_k1_float
        yx0_k = yx0 .=> yx0_k1_float

        Mcont = deepcopy(M_all[Int(evr[i,2])])
        Ncont = deepcopy(N_all[Int(evr[i,2])])

        Fx_pre = deepcopy(Fx_post)
        Fy_pre = deepcopy(Fy_post)
        Gx_pre = deepcopy(Gx_post)
        Gy_pre = deepcopy(Gy_post)
        J_pre = Fx_pre, Fy_pre, Gx_pre, Gy_pre

        f_pre = deepcopy(f_post)
        g_pre = deepcopy(g_post)
    end
end

function SortFDResults(res::Array{Float64},num_states)
  res = res'
  res_sort = Vector{Array{Float64}}(undef,0)
  len_t = Int(size(res)[2]/num_states)
  tmp_array = Array{Float64}(undef,0,len_t-2)
 
  for i=1:size(res)[1]
    for j=0:num_states-1
        tmp_array = [tmp_array;res[i,j*len_t+1+1:(j+1)*len_t-1]']
    end
    push!(res_sort,tmp_array)
    tmp_array = Array{Float64}(undef,0,len_t-2)
  end
  return res_sort
end
