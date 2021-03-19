using ModelingToolkit

function plot_sensi(time,sensis)
   p = plot(layout = (size(sensis)[1],1))
   for (index, value) in enumerate(sensis)
      plot!(time,value',subplot=index)
   end
   return p
end

function TimeDomainSensitivies(prob::ODEProblem,ic,p,sensis_u0,sensis_p)
    sol = solve(prob,Rodas4())
    TimeDomainSensitivies(prob,ic,p,sensis_u0,sensis_p,sol)
end

function TimeDomainSensitivies(prob::ODEProblem,ic,p,sensis_u0,sensis_p,sol)
   #sol,eqs,aeqs,D_states,A_states,symu0,symp,sensis_u0,sensis_p
    mtsys   = modelingtoolkitize(prob)
    fulleqs = equations(mtsys)
    state   = states(mtsys)
    params  = parameters(mtsys)

    aeqs = Vector{Equation}()
    eqs  = Vector{Equation}()
    A_states = Vector{Term{Real,Nothing}}()
    D_states = Vector{Term{Real,Nothing}}()
    for (index, value) in enumerate(fulleqs)
      if typeof(value.lhs) == Term{Real,Nothing} #works, but nasty
         push!(aeqs,value)
         push!(A_states,state[index])
      elseif value.lhs == 0
         push!(eqs,value)
         push!(D_states,state[index])
      else
         error("Can not interprete LHS of equation; $value")
      end
    end

    symu0 = Vector{Pair{Num,Float64}}(undef, size(ic))
    for (index,value) in enumerate(ic) symu0[index] = state[index] => ic[index]   end
    symp  = Vector{Pair{Num,Float64}}(undef, size(p))
    for (index,value) in enumerate(p) symu0[index] = params[index] => p[index]   end

    @derivatives D'~t

    len_sens = size(sensis_u0)[1]+size(sensis_p)[1];
    @parameters xx0[1:size(D_states)[1],1:len_sens]
    @parameters yx0[1:size(A_states)[1],1:len_sens]

    Diff_D_states = []
    Diff_A_states = []
    Diff_u0  = []
    Diff_p   = []

    for i in D_states push!(Diff_D_states, Differential(i)) end
    for i in A_states push!(Diff_A_states, Differential(i)) end
    for i in sensis_u0  push!(Diff_u0, Differential(i)) end
    for i in sensis_p   push!(Diff_p, Differential(i)) end

    Fx = Array{Expression}(undef,size(eqs)[1],size(D_states)[1])
    Fy = Array{Expression}(undef,size(eqs)[1],size(A_states)[1])
    Fp = Array{Expression}(undef,size(eqs)[1],len_sens)

    for (index_i, i) in enumerate(eqs)
            for (index_j, j) in enumerate(Diff_D_states)
               Fx[index_i,index_j] = expand_derivatives(j(i.rhs))
            end
            for (index_j, j) in enumerate(Diff_A_states)
               Fy[index_i,index_j] = expand_derivatives(j(i.rhs))
            end
            for (index_j, j) in enumerate(Diff_u0)
               Fp[index_i,index_j] = ModelingToolkit.Constant(0)
            end
            for (index_j, j) in enumerate(Diff_p)
               Fp[index_i,index_j+size(Diff_u0)[1]] = expand_derivatives(j(i.rhs))
            end
    end

    Gx = Array{Expression}(undef,size(aeqs)[1],size(D_states)[1])
    Gy = Array{Expression}(undef,size(aeqs)[1],size(A_states)[1])
    Gp = Array{Expression}(undef,size(aeqs)[1],len_sens)

    for (index_i, i) in enumerate(aeqs)
            for (index_j, j) in enumerate(Diff_D_states)
               Gx[index_i,index_j] = expand_derivatives(j(i.rhs))
            end
            for (index_j, j) in enumerate(Diff_A_states)
               Gy[index_i,index_j] = expand_derivatives(j(i.rhs))
            end
            for (index_j, j) in enumerate(Diff_u0)
               Gp[index_i,index_j] = ModelingToolkit.Constant(0)
            end
            for (index_j, j) in enumerate(Diff_p)
               Gp[index_i,index_j+size(Diff_u0)[1]] = expand_derivatives(j(i.rhs))
            end
    end

    M = [Δt/2*Fx-I Δt/2*Fy;
           Gx Gy]
    N = []
    if isempty(aeqs)
       N = -xx0-Δt/2*(Fx*xx0+Fy*yx0+Fp)
       else
       N = [-xx0-Δt/2*(Fx*xx0+Fy*yx0+Fp);
           zeros(size(A_states)[1],len_sens)]
    end
    O = [-Δt/2*Fp;
         -Gp]

    #Initialisierung: xx0 enthält die Sensis für x0 und p für x
    xx0_k     = Array{Pair{Operation,Float64}}(undef,size(xx0)[1],len_sens)
    xx0_f     = Array{Float64}(undef,size(xx0)[1],len_sens)
    for (index1, value1) in enumerate(sensis_u0) #die Startwertsensis wo dx1/dx1 = 1, sonst 0 z.B. dx1/dx2 =0
       for (index2, value2) in enumerate(D_states)
           if string(value1) == string(value2)
              xx0_f[index2,index1] = 1
              xx0_k[index2,index1] = xx0[index2,index1] => 1
           else
              xx0_f[index2,index1] = 0
              xx0_k[index2,index1] = xx0[index2,index1] => 0
           end
       end
    end
    for (index1, value1) in enumerate(sensis_p) #Alle Startwertsensis dx/dp = 0
       for (index2, value2) in enumerate(D_states)
              xx0_f[index2,index1+size(sensis_u0)[1]] = 0
              xx0_k[index2,index1+size(sensis_u0)[1]] = xx0[index2,index1+size(sensis_u0)[1]] => 0
       end
    end
    # Bei den Sensis für y werden zuerst die dy/x0 Sensi initialisiert
    yx0_t0 = -inv(Gy)*(Gx*xx0_f[:,1:size(sensis_u0)[1]])
    yp_t0 =  -inv(Gy)*(Gp*vcat(zeros(size(sensis_u0)[1],size(sensis_p)[1]),I))
    sym_yx0_k = [yx0_t0 yp_t0]
    yx0_k     = Array{Pair{Operation,Float64}}(undef,size(yx0)[1],len_sens)
    for (index, value) in enumerate(sym_yx0_k) yx0_k[index] = yx0[index] => substitute.(sym_yx0_k[index],([
      ; symp],))[1].value end

    Mfloat = Array{Float64}(undef,size(M)[1],size(M)[2])
    Nfloat = Array{Float64}(undef,size(N)[1],size(N)[2])
    Ofloat = Array{Float64}(undef,size(O)[1],size(O)[2])

    uk  = Vector{Pair{Operation,Float64}}(undef,size(symu0)[1])
    uk1 = Vector{Pair{Operation,Float64}}(undef,size(symu0)[1])

    for (index,value) in enumerate(symu0)
       uk[index] = value[1] => value[2]
       uk1[index] = value[1] => value[2]
    end

    sensi = Vector{Array{Float64}}(undef,len_sens)
    for i in 1:length(sensi) sensi[i] = Array{Float64}(undef,size(D_states)[1]+size(A_states)[1],size(sol)[2]-1) end

    for i in 1:size(sol)[2]-1
          dt = sol.t[i+1]-sol.t[i]

          for (index, value) in enumerate(sol.u[i])
             uk[index]  = uk[index][1] => value
          end
          for (index, value) in enumerate(sol.u[i+1])
             uk1[index] = uk1[index][1] => value
          end

          Mval = substitute.(M,([uk1; symp;Δt => dt],))
          Nval = substitute.(N,([uk;  symp;vec(xx0_k);vec(yx0_k);Δt => dt],))
          Oval = substitute.(O,([uk1; symp;Δt => dt],))

          for j in eachindex(Mval) Mfloat[j] = Mval[j].value end
          for j in eachindex(Nval) Nfloat[j] = Nval[j].value end
          for j in eachindex(Oval) Ofloat[j] = Oval[j].value end
          res  = inv(Mfloat)*(Nfloat+Ofloat)
          for j in 1:length(sensi) sensi[j][:,i] = res[:,j] end
          for j in 1:size(xx0_k)[2] #für jede Spalte
             for (index, value) in enumerate(xx0_k[:,j])
                xx0_k[index,j] = value[1] => res[index,j]
             end
          end
          for j in 1:size(yx0_k)[2] #für jede Spalte
             for (index, value) in enumerate(yx0_k[:,j])
                yx0_k[index,j] = value[1] => res[index+size(xx0_k)[2],j]
             end
          end
    end
   return sensi
end
