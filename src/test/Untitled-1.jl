
mtk = [mtk]
sol = pgsol0.dqsol
evr = evr_sol
s = s
h = hs

sym_states = states(mtk[1]);
sym_params = parameters(mtk[1]);

D0_states, A_states = GetFactorisedSymbolicStates(mtk[1]);
D_states = Num.(vcat(D0_states,sym_params));
A_states = Num.(A_states);

D0_indices= Int64.(setdiff(indexin(D0_states, sym_states), [nothing]));
A_indices = Int64.(setdiff(indexin(A_states, sym_states), [nothing]));

# INPUT: globale  Indices der D States und der Parameters 
# TESTEN: Was passiert wenn eins der beiden gleich null ist???
ind_d_sens = [17,18];
ind_p_sens = 6:7;

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
for (ind,val) in enumerate(ind_d_sens)
    d_state = sym_states[val]
    ind_tmp = Int64.(setdiff(indexin([d_state], D_states), [nothing]))[1]
    xx0_k[ind_tmp,ind] = xx0[ind_tmp] => 1.0
end

for (ind,val) in enumerate(ind_p_sens)
  p_state = sym_params[val]
  ind_tmp = Int64.(setdiff(indexin([p_state], D_states), [nothing]))[1]
  display(ind_tmp)
  xx0_k[ind_tmp,ind+length(ind_d_sens)] = xx0[ind_tmp] => 1.0
end


sym_all_y = Num.(A_states) .=> sol.prob.u0[A_indices]; 
sym_all_x = Num.(D_states) .=> [sol.prob.u0[D0_indices]; sol.prob.p];

#Initialisierung der algebraischen symbolischen Sensitivitäten
Gy_float = Float64.(Substitute(Gy_all[1], [sym_all_x; sym_all_y; t => sol.t[1]]))
Gx_float = Float64.(Substitute(Gx_all[1], [sym_all_x; sym_all_y; t => sol.t[1]]))
yx0f = -inv(Gy_float) * Gx_float * I[1:length(D_states),1:len_sens];
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

#return xx0_k,yx0_k,A_states,D_states, A_indices, D0_indices, matrices,Δt,len_sens, sensis


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

i=1
sol_part = sol[ind_sol[i,1]:ind_sol[i,2]]
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

for i = 2:size(sol0)[2]
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
end
#return sensi, xx0_k, yx0_k

for j = eachindex(sensi)
  sensis[j][:,ind_sol[i,1]+1:ind_sol[i,2]] = sensi[j]
end


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


f_pre_float = Float64.(Substitute(f_pre, [xk; yk;  t => sol.t[ind_sol[i,2]]]))
f_post_float = Float64.(Substitute(f_post, [xk1;yk1; t => sol.t[ind_sol[i,2]+1]]))

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
);


xy_k = [xk; yk;  t => sol.t[ind_sol[i,2]]]
xy_k1= [xk1;yk1; t => sol.t[ind_sol[i,2]+1]]
xx0_k_float=
f_pre=
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



syms = N[1];
subs_args = [xk1; yk1; Δt => dt; t => sol0.t[i]];
for i=1:5
  @time substitute.(syms, (subs_args,));
end
tmp = xk1[12:13]
@time tmp2 = substitute.(syms, (tmp,));
for i=1:5
  @time tmp1 = substitute.(tmp2, (subs_args,));
end
  

@time  substitute.(N[1], (subs_args,));
@time  substitute.(N0, (subs_args,));
