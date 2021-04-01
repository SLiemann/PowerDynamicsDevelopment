using PowerDynamics: rhs
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

function TimeDomainSensitivies(pg::PowerGrid,time_interval::Tuple{Float64,Float64},ic::Array{Float64,1},p::Array{Float64,1},sensis_u0_pd::Array{Symbol,1},sensis_p_pd::Array{Int64,1},sol::ODESolution)
    mtsys   = GetMTKSystem(pg,time_interval,ic,p)
    fulleqs = equations(mtsys)
    state   = states(mtsys)
    params  = parameters(mtsys)
    eqs,aeqs,D_states,A_states = GetSymbolicEquationsAndStates(fulleqs,state)

    #it is assumed that state and rhs(powergrid).syms have the same order
    sensis_u0 = state[indexin(sensis_u0_pd,rhs(pg).syms)]
    #sensis_p_pd is here a list with indices of the parameters p
    sensis_p = params[sensis_p_pd]

    #dict from states and parameters with their starting values
    symu0 = state .=> ic
    symp  = params .=> p

    Fx,Fy,Gx,Gy = GetSymbolicFactorizedJacobian(eqs,aeqs,D_states,A_states)

    Diff_u0  = Differential.(sensis_u0)
    Diff_p   = Differential.(sensis_p)
    len_sens = size(sensis_u0)[1]+size(sensis_p)[1];
    Fp = Array{Num}(undef,size(eqs)[1],len_sens)
    Gp = Array{Num}(undef,size(aeqs)[1],len_sens)

    Fp[:,1:size(Diff_u0)[1]] .= Num(0)
    Gp[:,1:size(Diff_u0)[1]] .= Num(0)
    for (ind, val) in enumerate(Diff_p)
      Fp[:,ind+size(Diff_u0)[1]] = Num.(expand_derivatives.(map(val,my_rhs.(eqs))))
      Gp[:,ind+size(Diff_u0)[1]] = Num.(expand_derivatives.(map(val,my_rhs.(aeqs))))
    end

    @parameters Δt
    @parameters xx0[1:size(D_states)[1],1:len_sens] #xx0 are the sensitivities regargind differential states
    @parameters yx0[1:size(A_states)[1],1:len_sens] #yx0 are the sensitivities regargind algebraic states
    M = [Δt/2*Fx-I Δt/2*Fy;
           Gx Gy]
    N = isempty(aeqs) ? -xx0-Δt/2*(Fx*xx0+Fy*yx0+Fp) : [-xx0-Δt/2*(Fx*xx0+Fy*yx0+Fp); zeros(size(A_states)[1],len_sens)]
    O = [-Δt/2*Fp;
         -Gp]

    #Initialisierung: xx0 enthält die Sensis für x0 und p für x
    xx0_k = xx0 .=> 0.0
    xx0_f = zeros(size(xx0)[1],len_sens)
    ind = indexin(sensis_u0,D_states)
    for i in 1:length(ind)
      xx0_k[i,ind[i]] = xx0_k[i,ind[i]][1] => 1.0
      xx0_f[i,ind[i]] = 1.0
    end
    # Bei den Sensis für y werden zuerst die dy/x0 Sensi initialisiert
    Gy_float = Substitute(Gy,[symu0; symp])
    # for increasing calculation of inv(Gy)
    yx0_t0 = -inv(Gy_float)*(Gx*xx0_f[:,1:size(sensis_u0)[1]])
    yp_t0  = -inv(Gy_float)*(Gp*vcat(zeros(size(sensis_u0)[1],size(sensis_p)[1]),I))
    yx0_k  = yx0 .=> Substitute([yx0_t0 yp_t0],[symu0; symp])

    sensi = Vector{Array{Float64}}(undef,len_sens)
    for i in 1:length(sensi) sensi[i] = Array{Float64}(undef,size(D_states)[1]+size(A_states)[1],size(sol)[2]-1) end

    for i in 1:size(sol)[2]-1
          dt = sol.t[i+1]-sol.t[i]

          uk  = state .=> sol.u[i]
          uk1 = state .=> sol.u[i+1]

          Mfloat = Substitute(M,[uk1; symp;Δt => dt])
          Nfloat = Substitute(N,[uk;  symp;vec(xx0_k);vec(yx0_k);Δt => dt])
          Ofloat = Substitute(O,[uk1; symp;Δt => dt])
          res  = inv(Mfloat)*(Nfloat+Ofloat)

          for j in 1:length(sensi)
             ind_y = setdiff(indexin(A_states,state), [nothing])
             ind_x = setdiff(indexin(D_states,state), [nothing])
             sensi[j][ind_x,i] = res[1:size(D_states)[1],j]
             sensi[j][ind_y,i] = res[size(D_states)[1]+1:end,j]
          end

          xx0_k = xx0 .=> res[1:size(D_states)[1],:]
          yx0_k = yx0 .=> res[size(D_states)[1]+1:end,:]
    end
   return sensi
end

function GetSymbolicEquationsAndStates(fulleqs::Array{Equation,1},state::Vector{Term{Real,Nothing}})
  aeqs = Vector{Equation}()
  eqs  = Vector{Equation}()
  A_states = Vector{Term{Real,Nothing}}()
  D_states = Vector{Term{Real,Nothing}}()
  for (index, value) in enumerate(fulleqs)
    if my_lhs(value) !==0
       push!(eqs,value)
       push!(D_states,state[index])
    elseif my_lhs(value) ===0
       push!(aeqs,value)
       push!(A_states,state[index])
    else
       error("Can not interprete LHS of equation; $value")
    end
  end
  return eqs,aeqs,D_states,A_states
end

function GetSymbolicFactorizedJacobian(eqs::Array{Equation,1},aeqs::Array{Equation,1},D_states::Array{Term{Real,Nothing},1},A_states::Array{Term{Real,Nothing},1})
  Fx = Array{Num}(undef,size(eqs)[1],size(D_states)[1])
  Fy = Array{Num}(undef,size(eqs)[1],size(A_states)[1])

  Gx = Array{Num}(undef,size(aeqs)[1],size(D_states)[1])
  Gy = Array{Num}(undef,size(aeqs)[1],size(A_states)[1])

  Diff_D_states = Differential.(D_states)
  Diff_A_states = Differential.(A_states)

  for (ind, val) in enumerate(Diff_D_states)
    Fx[:,ind] = Num.(expand_derivatives.(map(val,my_rhs.(eqs))))
    Gx[:,ind] = Num.(expand_derivatives.(map(val,my_rhs.(aeqs))))
  end
  for (ind, val) in enumerate(Diff_A_states)
    Fy[:,ind] = Num.(expand_derivatives.(map(val,my_rhs.(eqs))))
    Gy[:,ind] = Num.(expand_derivatives.(map(val,my_rhs.(aeqs))))
  end
  return Fx,Fy,Gx,Gy
end

function Substitute(syms::Array{Num},subs_args::Array{Pair{Num,Float64},1}) #SymbolicUtils.Symbolic{Real}
   return Symbolics.value.(substitute.(syms,(subs_args,)))
end

function CalcEigenValues(pg::PowerGrid,ic::Array{Float64,1},p::Array{Float64,1};output::Bool = false)
    mtsys     = GetMTKSystem(pg,(0.0,1.0),ic,p)
    fulleqs   = equations(mtsys)
    symstates = states(mtsys)
    x_eqs,y_eqs,x,y = GetSymbolicEquationsAndStates(fulleqs,symstates)
    Fx,Fy,Gx,Gy = GetSymbolicFactorizedJacobian(x_eqs,y_eqs,x,y)
    Fxf,Fyf,Gxf,Gyf = [Substitute(f,[symstates .=> ic0; parameters(mtsys) .=> 1]) for f in [Fx,Fy,Gx,Gy]]
    Af = Fxf -Fyf*inv(Gyf)*Gxf
    EW = eigvals(Af)
    if output
      println("|ID | Real-part | Imag-part | Frequency | Damping Time Constant |")
      index = indexin(x,symstates)
      syms  = rhs(pg).syms[index]
      for (ind,ew) in enumerate(EW)
        println("| $(syms[ind])) | $(round(real(ew),digits =3)) | $(round(imag(ew),digits = 3)) | $(round(abs(imag(ew))/2/pi,digits =3)) | $(round(1.0/abs(real(ew)),digits =3)) |")
      end
    end
    return EW
end

function GetMTKSystem(pg::PowerGrid,time_interval::Tuple{Float64,Float64},ic::Array{Float64,1},p::Array{Float64,1})
  #workaround for modelingtoolkitize (only Int64 entries for mass matrix)
  prob = ODEProblem(rhs(pg),ic,time_interval,p)
  new_f = ODEFunction(prob.f.f, syms = prob.f.syms, mass_matrix = Int.(prob.f.mass_matrix))
  ODEProb = ODEProblem(new_f,ic,time_interval,p)
  return modelingtoolkitize(ODEProb)
end
