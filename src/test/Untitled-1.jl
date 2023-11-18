

mtk = [mtk]
sol = pgsol0.dqsol
evr = evr_sol
s = s
h = hs

sym_states = states(mtk[1]);
sym_params = parameters(mtk[1]);

D0_states, A_states = GetFactorisedSymbolicStates(mtk[1]);
D_states = Num.(vcat(D0_states));
A_states = Num.(A_states);

D0_indices= Int64.(setdiff(indexin(D0_states, sym_states), [nothing]));
A_indices = Int64.(setdiff(indexin(A_states, sym_states), [nothing]));

#Hier anpassen, welche Zustände Parameter berücksichtigt werden soll
sensis_x = D_states; #anpassen bei variabler Größe
sensis_y = A_states; #anpassen bei variabler Größe?

#dict from states and parameters with their starting values
sym_x = sensis_x .=> sol.prob.u0[D0_indices];
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
for i= 1:length(sensis)
  sensis[i][length(D_states)+1:end,1] = yx0f[:,i]
end

return xx0_k,yx0_k,A_states,D_states, A_indices, D0_indices, matrices,Δt,len_sens, sensis







###############################################################################
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
Ncont = copy(N_all[1]);

#Bestimmung der partiellen Ableitungen für switchinig conditions and reset functions
@time hx,hy,sx,sy = CalcTriggerAndStateResetJacobians(s,h,D_states,A_states);

# Bestimmung der Zeitpunkte der Events im Loesungsobjekt + Anfang & Ende 
#ind_sol = vcat(2,setdiff(indexin(evr[:,1],sol.t),[nothing]),length(sol.t))

# ind1 = Anfang der Indizes für die kontinuierlichen Sensis
ind1 = vcat(1,setdiff(indexin(evr[:,1],sol.t),[nothing]).+1)
# ind2 = Ende der Indizes für die kontinuierlichen Sensis
ind2 = vcat(setdiff(indexin(evr[:,1],sol.t),[nothing]),length(sol.t))
# zwischen den Anfang und End Indizes liegt das Event mit Zustandsänderung
ind_sol = Int.(hcat(ind1,ind2))

@parameters t


τx0 = zeros(size(evr)[1],len_sens)

sol_part = sol[1:2];
sol0 = sol_part;
M0 = Mcont;
N0 = Ncont; 

#lokaler Speicher für die berechneten Trajektorien Senstivitäten
sensi = Vector{Array{Float64}}(undef, len_sens)
# die sensis sind hier eins kürzer, da ab k+1 bis zum Ende gerechnet wird!
# t_k ist der Zeitpunkt zuvor (also entweder Zeitpuntk des Sprungs oder Startzeit)
for i = 1:length(sensi)
  sensi[i] = 0.0.*Array{Float64}(
    undef,
    size(D_states)[1] + size(A_states)[1],
    size(sol0)[2]-1,
  )
end
xx0 = [i[1] for i in xx0_k]
yx0 = [i[1] for i in yx0_k]
#ind_y = Int64.(setdiff(indexin(A_states, sym_states), [nothing]))
#ind_x = Int64.(setdiff(indexin(D_states, sym_states), [nothing]))
@parameters t


i=2
  dt = sol0.t[i] - sol0.t[i-1]

  xk  = D_states .=> vcat(sol0[D0_indices,i-1],sol0.prob.p)
  xk1 = D_states .=> vcat(sol0[D0_indices,i],sol0.prob.p)

  yk  = A_states .=> sol0[A_indices,i-1]
  yk1 = A_states .=> sol0[A_indices,i]

  @time Mfloat = Float64.(Substitute(M0, [xk1; yk1; Δt => dt; t => sol0.t[i]]))
  @time Nfloat = Float64.(Substitute(N0, [xk; yk; vec(xx0_k); vec(yx0_k); Δt => dt; t => sol0.t[i-1]]))
  @time res = inv(Mfloat) * Nfloat 

  for j = 1:length(sensi)
    sensi[j][1:length(D_states), i-1] = res[1:length(D_states), j]
    sensi[j][length(D_states)+1:end, i-1] = res[length(D_states)+1:end, j]
  end
  # Lesezeichen
  xx0_k = xx0 .=> res[1:size(D_states)[1],:]
  yx0_k = yx0 .=> res[size(D_states)[1]+1:end,:]


Substitute(M0, [xk1; yk1; Δt => dt; t => sol0.t[i]])
typeof(M0)
typeof([xk1; yk1; Δt => dt; t => sol0.t[i]])

syms = M0
subs_args = [xk1; yk1; Δt => dt; t => sol0.t[i]]
@time substitute.(syms, (subs_args,))
@time Symbolics.value.(substitute.(syms, (subs_args,)))

tmp = [xk1[1]]
@time tmp2 = substitute.(syms, (tmp,));
