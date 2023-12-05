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
#   3.1 [t_i mtk_i s_i h_i parameters]
#   3.1.1 t_i = Zeit, mtk_i = 'aktives Netz' bei Änd. der alg. Gl., s_i = switch cond., h_i = reset function
#   3.2 Initialisierung: evr = Array{Float64}(undef,0,4)
#   3.3 Eintrag machen (Bsp.) evr = vcat(evr, [integrator.t 1 1 1 integrator.t])
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
    if my_lhs(value) === 0 || my_lhs(value) === 0.0
      push!(aeqs, value)
      push!(A_states, state[index])
    elseif my_lhs(value) !== 0 || my_lhs(value) !== 0.0
        push!(eqs, value)
        push!(D_states, state[index])
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

GetDiscreteStatePositions(mtsys::ODESystem) = GetDiscreteStatePositions(equations(mtsys))
function GetDiscreteStatePositions(fulleqs::Array{Equation,1})
  ind = Vector{Int64}()
  for (index, value) in enumerate(fulleqs)
    if my_lhs(value) !== 0 && my_rhs(value) === 0.0
      push!(ind, index)
    end
  end
  return ind
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
#function Substitute(syms::Matrix{Num}, subs_args::Dict{Any, Any}) #SymbolicUtils.Symbolic{Real} ::Array{Pair{Num,Float64},1}
#  return Symbolics.value.(substitute.(syms, (subs_args,)))
#end
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

function InitTrajectorySensitivity(mtk::Vector{ODESystem},sol::ODESolution,ind_d_sens::Vector{Int64},ind_p_sens::Vector{Int64})
  sym_states = states(mtk[1]);
  sym_params = parameters(mtk[1]);

  D0_states, A_states = GetFactorisedSymbolicStates(mtk[1]);
  D_states = Num.(vcat(D0_states,sym_params));
  A_states = Num.(A_states);

  D0_indices= Int64.(setdiff(indexin(D0_states, sym_states), [nothing]));
  A_indices = Int64.(setdiff(indexin(A_states, sym_states), [nothing]));

  #Hier anpassen, welche Zustände Parameter berücksichtigt werden soll
  sensis_x = Num.(vcat(sym_states[ind_d_sens],sym_params[ind_p_sens])); #anpassen bei variabler Größe

  #dict from states and parameters with their starting values
  sym_x = sensis_x .=> vcat(sol.prob.u0[ind_d_sens],sol.prob.p[ind_p_sens])
  len_sens = length(sym_x) # #den Einfluss der algebraischen Variablen steckt indirekt in den Gleichungen!!

  @parameters Δt t
  @parameters xx0[1:length(D_states), 1:len_sens] #xx0 are the symbolic sensitivities regarding differential states
  @parameters yx0[1:length(A_states), 1:len_sens] #yx0 are the symbolic sensitivities regarding algebraic states
  xx0 = Symbolics.scalarize(xx0)
  yx0 = Symbolics.scalarize(yx0)

  # die Sensitivität eines Parameters/Zustand muss immer auf alle Zustände/Parameter berechnet werden: x_x1, p_x1, y_x1 
  # somit müssen die Ableitungen für ALLE Zustände berechnet werden!
  matrices =  GetEqsJacobianSensMatrices(mtk,xx0,yx0); #anpassen wenn nicht alle Variablen gerechnet werden sollen
  F_all,G_all,J_all,M,N = matrices;
  Fx_all,Fy_all,Gx_all,Gy_all = J_all;

  #Initialisierung der differntialen symbolischen Sensitivitäten: xx0 enthält die Sensis für x, q und p 
  xx0_k = xx0 .=> 0.0
  xx0_float = zeros(size(xx0))
  ind_dp_sens = zeros(Int64,len_sens,1)
  for (ind,val) in enumerate(ind_d_sens)
      d_state = sym_states[val]
      ind_tmp = Int64.(setdiff(indexin([d_state], D_states), [nothing]))[1]
      ind_dp_sens[ind] = deepcopy(ind_tmp)
      xx0_k[ind_tmp,ind] = xx0[ind_tmp,ind] => 1.0
      xx0_float[ind_tmp,ind] = 1.0
  end

  for (ind,val) in enumerate(ind_p_sens)
    p_state = sym_params[val]
    ind_tmp = Int64.(setdiff(indexin([p_state], D_states), [nothing]))[1]
    ind_dp_sens[ind+length(ind_d_sens)] = deepcopy(ind_tmp)
    xx0_k[ind_tmp,ind+length(ind_d_sens)] = xx0[ind_tmp,ind+length(ind_d_sens)] => 1.0
    xx0_float[ind_tmp,ind+length(ind_d_sens)] = 1.0
  end

  sym_all_y = Num.(A_states) .=> sol.prob.u0[A_indices]; 
  sym_all_x = Num.(D_states) .=> [sol.prob.u0[D0_indices]; sol.prob.p];

  #Initialisierung der algebraischen symbolischen Sensitivitäten
  Gy_float = Float64.(Substitute(Gy_all[1], [sym_all_x; sym_all_y; t => sol.t[1]]))
  Gx_float = Float64.(Substitute(Gx_all[1], [sym_all_x; sym_all_y; t => sol.t[1]]))
  yx0f = -inv(Gy_float) * Gx_float * xx0_float; #wichtig hier mit xx0_float zu multiplizieren (bei Hisken wird dies als Einheitsmatrix angenommne)
  yx0_k = yx0 .=> yx0f

  #Objekterstelleung der Trajektorien Sensitivität
  sensis = Vector{Array{Float64}}(undef, len_sens)
  for i = 1:length(sensis)
    sensis[i] = zeros(Float64,size(D_states)[1] + size(A_states)[1],size(sol)[2])
  end
  # Initialisierung der Trajektorien Sensitivität
  for (ind,val) in enumerate(ind_dp_sens)
    sensis[ind][val,1] = 1.0
  end
  for i= 1:length(sensis)
    sensis[i][length(D_states)+1:end,1] = yx0f[:,i]
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
  discrete_state_pos::Vector{Pair{Num, Float64}}
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

    return xx0_post, yx0_post, τx0
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

# function CalcContinuousSensitivity(mtk::Vector{ODESystem},sol::ODESolution)
#     xx0_k,yx0_k,A_states,D_states, A_indices, D0_indices,matrices,Δt,len_sens, sensis = InitTrajectorySensitivity(mtk,sol);
#     f_all, g_all, J_all, M_all, N_all = matrices

#     sensi_part,xx0_k,yx0_k = ContinuousSensitivity(
#       sol,
#       xx0_k,
#       yx0_k,
#       A_states, 
#       D_states, 
#       A_indices,
#       D0_indices,
#       M_all[1], 
#       N_all[1], 
#       Δt,
#       len_sens)

#     return sensi_part
# end

function CalcHybridTrajectorySensitivity(
  mtk::Vector{ODESystem},
  sol::ODESolution,
  evr::Matrix{Float64},
  s::Vector{Num},
  h::Vector{Vector{Num}},
  ind_d_sensi::Vector{Int64},
  ind_p_sensi::Vector{Int64}
  )

    xx0_k,yx0_k,A_states,D_states, A_indices, D0_indices,matrices,Δt,len_sens, sensis = InitTrajectorySensitivity(mtk,sol,ind_d_sensi,ind_p_sensi);

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

    # ALTE HERANGEHENSWEISE
    # # ind1 = Anfang der Indizes für die kontinuierlichen Sensis
    # ind1 = vcat(1,setdiff(indexin(evr[:,1],sol.t),[nothing]).+1)
    # # ind2 = Ende der Indizes für die kontinuierlichen Sensis
    # ind2 = vcat(setdiff(indexin(evr[:,1],sol.t),[nothing]),length(sol.t))
    # # zwischen den Anfang und End Indizes liegt das Event mit Zustandsänderung
    # ind_sol = Int.(hcat(ind1,ind2))

    # NEUE HERANGEHENSWEIS hier wird bei jedem zu beachtenden Callback die Option save_positions=(false,true)
    # gesetzt damit der letzte gefundende Index genau dem Sprung nach dem Callback darstellt. 
    ind_sol = zeros(Int64,size(evr)[1]+1,2)
    ind_sol[1,1] = 1;
    for (ind,val) in enumerate(evr[:,1])
        tmp_ind = findall(==(val),pgsol0.dqsol.t) 
        ind_sol[ind,2] = tmp_ind[end-1]
        ind_sol[ind+1,1] = tmp_ind[end]
    end
    ind_sol[end,2] = length(pgsol0.dqsol.t);

    @parameters t

    τx0 = zeros(size(evr)[1],len_sens)

    ind_q = GetDiscreteStatePositions(mtk[1]); # sollte überall gleich sein.
    sym_params = parameters(mtk[1]); # sollte überall gleich sein.
    all_states = states(mtk[1]);# sollte überall gleich sein.

    for i = 1:size(ind_sol)[1]
        sol_part = sol[ind_sol[i,1]:ind_sol[i,2]]
        # Hier ist viel Potential für Fehler!!
        p_selection = sol.prob.p; #nur Initialisierung
        if i !=1
          p_selection = evr[i-1,5:end]
        end

        presubs = vcat(Num.(all_states[ind_q]),sym_params) .=> [sol_part[1][ind_q]; p_selection]
        display("Zeitpunkt der kontinuierlichen Sensis: $(sol_part.t[1])");
        display("Werte der diskreten Zustände: $(sol_part[2][ind_q]) von $(all_states[ind_q])");
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
                                                      len_sens,
                                                      presubs)
        for j = eachindex(sensi_part)
            sensis[j][:,ind_sol[i,1]+1:ind_sol[i,2]] = sensi_part[j]
        end

        if ind_sol[i,2] ≠ ind_sol[end,2]
            display("Event at t = $(evr[i,1])")

            xx0_k_float = [i[2] for i in xx0_k]

            xk = D_states .=> vcat(sol[D0_indices,ind_sol[i,2]],p_selection)    #hier τ-
            xk1 = D_states .=> vcat(sol[D0_indices,ind_sol[i,2]+1],evr[i,5:end]) #hier τ+
        
            yk  = A_states .=> sol[A_indices,ind_sol[i,2]]   #hier τ-
            yk1 = A_states .=> sol[A_indices,ind_sol[i,2]+1] #hier τ+
            
            Fx_post = Fx_all[Int(evr[i,2])]
            Fy_post = Fy_all[Int(evr[i,2])]
            Gx_post = Gx_all[Int(evr[i,2])]
            Gy_post = Gy_all[Int(evr[i,2])]
            J_post = Fx_post, Fy_post, Gx_post, Gy_post

            f_post = f_all[Int(evr[i,2])]
            g_post = g_all[Int(evr[i,2])]

            if !(iszero(evr[i,4])) && !(iszero(evr[i,3])) #should only be calculated by a real event (not a fake one, e.g. post-fault system)
              hx_tmp = hx[Int(evr[i,4])]
              hy_tmp = hy[Int(evr[i,4])]
              sx_tmp = sx[Int(evr[i,3])]
              sy_tmp = sy[Int(evr[i,3])]

              xx0_k1_float, yx0_k1_float, tmp_τx0 = CalcSensitivityAfterJump(
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
            else
              # if not a real event occurs
              xx0_k1_float = xx0_k_float
              yx0_k1_float = [i[2] for i in yx0_k]
              tmp_τx0 = zeros(1,len_sens) #no jump, no shitft
            end

            # speichern der Sensis nach dem Sprung
            for k in 1:size(xx0_k1_float)[2]
              sensis[k][:,ind_sol[i,2]+1] = vcat(xx0_k1_float[:,k],yx0_k1_float[:,k])
            end
            τx0[i,:] = tmp_τx0
        else
            #deleting parameter sensitivities
            len_params = length(parameters(mtk[1]))
            len_d0 = length(D0_indices)
            len_states =  length(states(mtk[1]))
            tmp = zeros(Float64,len_states,length(sol))
            for k in 1:length(sensis)
              tmp[D0_indices,:] = sensis[k][1:len_d0,:] #x_x0
              tmp[A_indices,:] = sensis[k][len_d0+len_params+1:end,:] #y_x0
              sensis[k] = deepcopy(tmp)
            end
            return sensis, τx0
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


function ApproximatedTrajectory(sol::ODESolution,sens::VecOrMat{Float64},Δx0::Float64)
  ApproximatedTrajectory(sol[:,:],sens,Δx0)
end

function ApproximatedTrajectory(sol::VecOrMat{Float64},sens::VecOrMat{Float64},Δx0::Float64)
  sol_appr = similar(sol[:,:])
  for i in 1:size(sol)[2]
    sol_appr[:,i] = sol[:,i] + sens[:,i].*Δx0
  end
  return sol_appr
end

function TrajectoryRefinement(mtk::Vector{ODESystem},sol::ODESolution,evr::Matrix{Float64},sensis::Vector{Array{Float64}},Δτ::Matrix{Float64},x0_ind::Int64,Δx0::Float64)

  # ind1 = Anfang der Indizes für die kontinuierlichen Sensis
  ind1 = vcat(1,setdiff(indexin(evr[:,1],sol.t),[nothing]).+1)
  # ind2 = Ende der Indizes für die kontinuierlichen Sensis
  ind2 = vcat(setdiff(indexin(evr[:,1],sol.t),[nothing]),length(sol.t))
  # zwischen den Anfang und End Indizes liegt das Event mit Zustandsänderung
  ind_sol = Int.(hcat(ind1,ind2,zeros(size(ind1)[1],2)))
  display(ind_sol)
  for i=1:size(ind_sol)[1]-1
      if Δτ[i,x0_ind]*Δx0 >= 0.0
          # WENN EIN EVENT SHIFT ÜBER EIN ANDERES EVENT DRÜBER GEHT  
          # MÜSSTE MAN EIGENTLICH DIE INDIZES WIEDER ANPASSEN UND NICHT SO WEIT GEHEN.
          # ANSONSTEN WERDEN WIEDER ALTE WERTE VERÄNDERT
          ind = findall(>(0),sol.t[ind_sol[i,2]+1] .< sol.t .< sol.t[ind_sol[i,2]+1] + Δτ[i,x0_ind]*Δx0)
          if isempty(ind) # if the shift in time is not that big to reach the next solution point
              ind_sol[i,3] = ind_sol[i,2]+1
              ind_sol[i,4] = ind_sol[i,2]+1
          else
              ind_sol[i,3] = ind[1]-1
              ind_sol[i,4] = ind[end]
          end
      else
          ind = findall(>(0),sol.t[ind_sol[i,2]+1] + Δτ[i,x0_ind]*Δx0 .< sol.t .<= sol.t[ind_sol[i,2]+1])
          if isempty(ind)
              ind_sol[i,3] = ind_sol[i,2]+1
              ind_sol[i,4] = ind_sol[i,2]+1
          else
              ind_sol[i,2] = ind[1]-1
              ind_sol[i,3] = ind[1]
              ind_sol[i,4] = ind[end]-1
          end
      end
  end
  display(ind_sol)

  f_all = Vector{Array{Num,1}}(undef,length(mtk))
  Fx_all = Vector{Array{Num,2}}(undef,length(mtk))
  Fy_all = similar(Fx_all)
  Gx_all = similar(Fx_all)
  Gy_all = similar(Fx_all)
  for (ind,val) in enumerate(mtk)
      fulleqs = equations(val)
      symstates = states(val)
      sym_params = parameters(val)
      eqs, aeqs, x, y = GetSymbolicEquationsAndStates(fulleqs, symstates)

      Fx_all[ind],Fy_all[ind],Gx_all[ind],Gy_all[ind] = GetSymbolicFactorizedJacobian(eqs, aeqs, x, y)
      f_all[ind] = my_rhs.(eqs)
  end
  sol_refin = zeros(size(sol))

  #get indices
  sym_states = states(mtk[1]);
  sym_params = parameters(mtk[1]);

  D0_states, A_states = GetFactorisedSymbolicStates(mtk[1]);

  D0_indices= Int64.(setdiff(indexin(D0_states, sym_states), [nothing]))
  A_indices = Int64.(setdiff(indexin(A_states, sym_states), [nothing]))

  pk = sym_params .=> sol.prob.p

  for i=1:size(ind_sol)[1]
      if i !=size(ind_sol)[1]
          if Δτ[i,x0_ind]*Δx0 >= 0.0
              
              f = f_all[Int(evr[i,2])]
              Gx = Gx_all[Int(evr[i,2])]
              Gy = Gy_all[Int(evr[i,2])]
              
              #Approximation until jump
              if i == 1
                  tmp_sol = sol[:,ind_sol[i,1]:ind_sol[i,2]]
                  tmp_sens = sensis[x0_ind][:,ind_sol[i,1]:ind_sol[i,2]]
                  #display(length(tmp_sol))
                  sol_refin[:,ind_sol[i,1]:ind_sol[i,2]] = ApproximatedTrajectory(tmp_sol,tmp_sens,Δx0)
              else # the solution and sensis after the new jumps have to be used
                  tmp_sol = sol[:,ind_sol[i-1,4]+1:ind_sol[i,2]]
                  tmp_sens = sensis[x0_ind][:,ind_sol[i-1,4]+1:ind_sol[i,2]]
                  sol_refin[:,ind_sol[i-1,4]+1:ind_sol[i,2]] = ApproximatedTrajectory(tmp_sol,tmp_sens,Δx0)
              end

              xk = sym_states .=> sol[:,ind_sol[i,2]]  # the last values BEFORE the original jump -> are used
              f_f  = Float64.(Substitute(f,[xk;pk]))   # for all approximations within t+Δτ
              Gx_f = Float64.(Substitute(Gx,[xk;pk]))
              Gy_f = Float64.(Substitute(Gy,[xk;pk]))

              gygx = inv(Gy_f)*Gx_f*f_f

              x_part = sol_refin[D0_indices,ind_sol[i,2]]
              y_part = sol_refin[A_indices, ind_sol[i,2]]
  
              for (ind,j) in enumerate(sol.t[ind_sol[i,3]:ind_sol[i,4]])
                  dt = j .- sol.t[ind_sol[i,2]]
                  sol_refin[D0_indices, ind_sol[i,3]+ind-1] = x_part + f_f  * dt
                  sol_refin[A_indices,ind_sol[i,3]+ind-1] = y_part + gygx * dt
              end
          else
              #for negative values
              # maybe thats not right here, because the evr saves the NEW active 
              # set of equations, or? So evr[i,] would contain already f+ etc.
              # For Δτ  than a evr[i-1] have to be used and evr[1] as first entry.
              f = f_all[Int(evr[i,2])] # here f+ etc. have to be taken
              Gx = Gx_all[Int(evr[i,2])]
              Gy = Gy_all[Int(evr[i,2])]
              
              #Approximation until jump
               if i == 1
                  tmp_sol = sol[:,ind_sol[i,1]:ind_sol[i,2]]
                  tmp_sens = sensis[x0_ind][:,ind_sol[i,1]:ind_sol[i,2]]
                  sol_refin[:,ind_sol[i,1]:ind_sol[i,2]] = ApproximatedTrajectory(tmp_sol,tmp_sens,Δx0)
               else # the solution and sensis after the new jumps have to be used
                  tmp_sol = sol[:,ind_sol[i-1,4]+1:ind_sol[i,2]]
                  tmp_sens = sensis[x0_ind][:,ind_sol[i-1,4]+1:ind_sol[i,2]]
                  sol_refin[:,ind_sol[i-1,4]+1:ind_sol[i,2]] = ApproximatedTrajectory(tmp_sol,tmp_sens,Δx0)
               end
              xk = sym_states .=> sol[:,ind_sol[i,4]+1]  # the first values after the original jump -> are used
              f_f  = Float64.(Substitute(f,[xk;pk]))   # for all approximations within t+Δτ
              Gx_f = Float64.(Substitute(Gx,[xk;pk]))
              Gy_f = Float64.(Substitute(Gy,[xk;pk]))

              gygx = inv(Gy_f)*Gx_f*f_f
              x_part = sol[D0_indices,ind_sol[i,4]+1] + sensis[x0_ind][D0_indices,ind_sol[i,4]+1] * Δx0
              y_part = sol[A_indices,ind_sol[i,4]+1]  + sensis[x0_ind][A_indices,ind_sol[i,4]+1] * Δx0
  
              for (ind,j) in enumerate(sol.t[ind_sol[i,3]:ind_sol[i,4]])
                  dt = j .- sol.t[ind_sol[i,4]+1]
                  sol_refin[D0_indices, ind_sol[i,3]+ind-1] = x_part + f_f  * dt
                  sol_refin[A_indices,ind_sol[i,3]+ind-1] = y_part + gygx * dt
              end
          end
      else
              tmp_sol = sol[:,ind_sol[i-1,4]+1:ind_sol[i,2]]
              tmp_sens = sensis[x0_ind][:,ind_sol[i-1,4]+1:ind_sol[i,2]]
              sol_refin[:,ind_sol[i-1,4]+1:ind_sol[i,2]] = ApproximatedTrajectory(tmp_sol,tmp_sens,Δx0)
              # p_sol = sol[:,ind_sol[i,1]:ind_sol[i,2]]
              # tmp_sens = sensis[x0_ind][:,ind_sol[i,1]:ind_sol[i,2]]
              # sol_refin[:,ind_sol[i,1]:ind_sol[i,2]] = ApproximatedTrajectory(tmp_sol,tmp_sens,Δx0)
      end
  end
  return sol_refin
end