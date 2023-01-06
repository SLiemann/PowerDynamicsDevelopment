using PowerDynamics: _find_operationpoint_rootfind
using PowerDynamics: rhs#, symbolsof
using NLsolve: nlsolve, converged
using IfElse
#include("PowerFlow.jl") # for NodalAdmittanceMatrice

function InitializeInternalDynamics(pg::PowerGrid,ic_lf::Array{Float64,1}) # ,I_c::Matrix{Complex{Float64}})
   Ykk = NodalAdmittanceMatrice(pg)
   Uc  = getComplexBusVoltage(pg,ic_lf)
   I_c = Ykk*Uc

   ind_offset = 1
   pg_new = deepcopy(pg)
   ic0 = deepcopy(ic_lf)
   for (ind,val) in enumerate(pg.nodes)
      len_node_dynamics = length(symbolsof(val[2]))
       if len_node_dynamics != 2
          ic0[ind_offset:ind_offset+len_node_dynamics-1],node = InitNode(val[2],ind,I_c,ic_lf,ind_offset)
          pg_new.nodes[val[1]] = node #Update powergrid
       end
       ind_offset += len_node_dynamics
   end
   return pg_new,ic0
end

function InitNode(SM::FourthOrderEq,ind::Int64,I_c::Vector{Complex{Float64}},ic_lf::Array{Float64,1},ind_offset::Int64)
   v_d_temp = ic_lf[ind_offset]
   v_q_temp = ic_lf[ind_offset+1]
   #Rotor angle
   δ = angle(v_d_temp+1im*v_q_temp+(+1im*(SM.X_q - SM.X_q_dash))*I_c[ind])

   v   = v_d_temp +1im*v_q_temp
   v   = 1im*v*exp(-1im*δ)
   v_d = real(v)
   v_q = imag(v)
   i   = 1im*I_c[ind]*exp(-1im*δ)
   i_d = real(i)
   i_q = imag(i)

   V_f = v_q + (SM.X_d - SM.X_d_dash) * i_d
   pe = SM.P + (SM.X_q_dash - SM.X_d_dash)*i_d*i_q
   node_temp = FourthOrderEq(H=SM.H, P=pe, D=SM.D, Ω=SM.Ω, E_f=V_f, T_d_dash=SM.T_d_dash ,T_q_dash=SM.T_q_dash ,X_q_dash =SM.X_q_dash ,X_d_dash=SM.X_d_dash,X_d=SM.X_d, X_q=SM.X_q)
   return [v_d_temp, v_q_temp, δ, 0.], node_temp
end

function InitNode(DG::GridSideConverter,ind::Int64,I_c::Vector{ComplexF64},ic_lf::Array{Float64,1},ind_offset::Int64)
   v_d_temp = ic_lf[ind_offset]
   v_q_temp = ic_lf[ind_offset+1]
   u0 = v_d_temp+1im*v_q_temp
   i0 = I_c[ind]

   idq = i0*exp(-1im*angle(u0))
   id_temp = real(idq)
   iq_temp = imag(idq)

   # x,y and z are only initialized the right way,
   # if the dynamic limiter is not limiting the output at the stationary point
   # Other cases has to be added later on

   id0 = id_temp/DG.Kgsc
   x_st_temp = id0*DG.Tp
   iq0 = iq_temp/DG.Kgsc
   y_st_temp = iq0*DG.Tq
   z_st_temp = iq0*DG.Tv

   node_temp = GridSideConverter(mode=DG.mode, p_ref=DG.p_ref, q_ref=DG.q_ref, v_ref=DG.v_ref,
        idmax=DG.idmax, iqmax=DG.iqmax, imax=DG.imax,
        Kp=DG.Kp, Tp=DG.Tp, Kq=DG.Kq, Tq=DG.Tq,
        Kv=DG.Kv, Tv=DG.Tv, Kgsc=DG.Kgsc, Tgsc=DG.Tgsc,
        δqv=DG.δqv, v1_max=DG.v1_max, v1_min=DG.v1_min, q_max=DG.q_max);

   # Structure from node: x, y, z, id, iq]
   return [v_d_temp, v_q_temp, x_st_temp, y_st_temp, z_st_temp, id_temp, iq_temp], node_temp
end
function InitNode(GS::SynchronousMachineGENSAL,ind::Int64,I_c::Vector{Complex{Float64}},ic_lf::Array{Float64,1},ind_offset::Int64)
   a = 0.2; b = 2*GS.S12-1.2*2*GS.S10; c = 1.2*GS.S10^2-GS.S12^2
   A = (-b + sqrt(b^2 -4*a*c))/(2*a)
   B = 1/(A^2 - 2*A*GS.S10 + GS.S10^2)
   #sf(x) = B*(x-A)^2

   v_r_temp = ic_lf[ind_offset]
   v_i_temp = ic_lf[ind_offset+1]

   δ_temp = angle(v_r_temp+1im*v_i_temp+(GS.R_a+1im*GS.x_q)*I_c[ind]/(GS.Srated/GS.Sbase))
   print("δ_temp: ", δ_temp)

   v = v_r_temp +1im*v_i_temp
   v = 1im*v*exp(-1im*δ_temp)
   v_d = real(v)
   v_q = imag(v)
   i   = 1im*I_c[ind]*exp(-1im*δ_temp)/(GS.Srated/GS.Sbase)
   i_d = real(i)
   i_q = imag(i)

   ω_temp = 0.
   ψ_qss_temp = -v_d
   print(" ψ_qss_temp = -v_d = ", -v_d, "\n")
   print(" ψ_qss_temp = -i_q*(x_q - x_qss) = ", -i_q*(GS.x_q - GS.x_qss), "\n")
   ψ_dss_temp = v_q

   #E_qs_temp = -(1/B - 2*A/B)/2 + sqrt(abs(((1/B - 2*A/B)/2)^2 + 1/B*(GS.E_fd-B*A^2-(GS.x_d-GS.x_ds)*i_d)))
   E_qs_temp = GS.E_fd +(GS.x_ds-GS.x_d)*i_d
   ψ_ds_temp = (ψ_dss_temp - (GS.x_dss-GS.x_l)/(GS.x_ds-GS.x_l)*E_qs_temp)*(GS.x_ds-GS.x_l)/(GS.x_ds-GS.x_dss)

   node_temp = SynchronousMachineGENSAL(Sbase=GS.Sbase, Srated=GS.Srated, D=GS.D, H=GS.H, P=GS.P, E_fd=GS.E_fd, R_a=GS.R_a, x_d=GS.x_d,
                                       x_q=GS.x_q, x_ds=GS.x_ds, x_dss=GS.x_dss, x_qss=GS.x_qss, x_l=GS.x_l,
                                       T_d0s=GS.T_d0s, T_d0ss=GS.T_d0ss, T_q0ss=GS.T_q0ss, S10=GS.S10, S12=GS.S12);

   # [ω, dω], [δ, dδ], [E_qs, dE_qs], [ψ_ds, dψ_ds], [ψ_qss, dψ_qss]
   return [v_r_temp, v_i_temp, ω_temp, δ_temp, E_qs_temp, ψ_ds_temp, ψ_qss_temp], node_temp
end

function InitNode(SM::SixOrderMarconatoMachine,ind::Int64,I_c::Vector{Complex{Float64}},ic_lf::Array{Float64,1},ind_offset::Int64)
   v_d_temp = ic_lf[ind_offset]
   v_q_temp = ic_lf[ind_offset+1]
   #Rotor angle
   δ = angle(v_d_temp+1im*v_q_temp+(SM.R_a+1im*SM.X_q)*I_c[ind]/(SM.Srated/SM.Sbase)) #+angle(v_d_temp + 1im*v_q_temp)

   v = v_d_temp +1im*v_q_temp
   v = 1im*v*exp(-1im*δ)
   v_d = real(v)
   v_q = imag(v)
   i   = 1im*I_c[ind]*exp(-1im*δ)/(SM.Srated/SM.Sbase)
   i_d = real(i)
   i_q = imag(i)

   #Converstion of short-circuit time constants to open-loop time constants
   T_d0s = SM.T_ds*(SM.X_d/SM.X_ds)
   T_q0s = SM.T_qs*(SM.X_q/SM.X_qs)
   T_d0ss = SM.T_dss*(SM.X_ds/SM.X_dss)
   T_q0ss = SM.T_qss*(SM.X_qs/SM.X_qss)

   #Auxillary variables
   γ_d = T_d0ss * SM.X_dss * (SM.X_d - SM.X_ds) / (T_d0s * SM.X_ds)
   γ_q = T_q0ss * SM.X_qss * (SM.X_q - SM.X_qs) / (T_q0s * SM.X_qs)

   e_qss = v_q + SM.R_a * i_q + SM.X_dss * i_d
   e_dss = v_d + SM.R_a * i_d - SM.X_qss * i_q

   e_ds =  (SM.X_q - SM.X_qs - γ_q) * i_q
   #-(SM.X_qs - SM.X_qss + γ_q) * i_q + e_dss #also valid

   e_qs = 0.
   v_f =  1.
   if SM.T_AA != 0.
      e_qs = (SM.T_AA/T_d0s) * (-(SM.X_d - SM.X_ds - γ_d) *i_d + (T_d0s/SM.T_AA - 1)*(e_qss + (SM.X_ds - SM.X_dss + γ_d) *i_d ))
      v_f  = (T_d0s / SM.T_AA) * (e_qss - e_qs + (SM.X_ds - SM.X_dss + γ_d) *i_d)
   else
      e_qs = e_qss + (SM.X_ds - SM.X_dss + γ_d) * i_d
      v_f  = e_qs  + (SM.X_d  - SM.X_ds  - γ_d) * i_d
   end

   #Pm also needs to be initialized
   Pm = (v_q + SM.R_a * i_q) * i_q + (v_d + SM.R_a * i_d) * i_d

   #Create new bus
   node_temp = SixOrderMarconatoMachine(Sbase=SM.Sbase,Srated=SM.Srated,H=SM.H, P=Pm, D=SM.D, Ω=SM.Ω,
                                           E_f=v_f, R_a=SM.R_a, T_ds=SM.T_ds, T_qs=SM.T_qs, T_dss=SM.T_dss,
                                           T_qss=SM.T_qss, X_d=SM.X_d, X_q=SM.X_q, X_ds=SM.X_ds, X_qs=SM.X_qs,
                                           X_dss=SM.X_dss, X_qss=SM.X_qss, T_AA=SM.T_AA);

   # Structure from node: u_r, u_i, θ, ω, e_ds, e_qs, e_dss,e_qss
   return [v_d_temp, v_q_temp, δ, 0., e_ds, e_qs, e_dss, e_qss], node_temp
end

function InitNode(SM::SixOrderMarconatoMachineAVROEL,ind::Int64,I_c::Vector{Complex{Float64}},ic_lf::Array{Float64,1},ind_offset::Int64)
   v_d_temp = ic_lf[ind_offset]
   v_q_temp = ic_lf[ind_offset+1]#*15/380
   #Rotor angle
   δ = angle(v_d_temp+1im*v_q_temp+(SM.R_a+1im*SM.X_q)*I_c[ind]/(SM.Srated/SM.Sbase))

   v = v_d_temp +1im*v_q_temp
   v = 1im*v*exp(-1im*δ)
   v_d = real(v)
   v_q = imag(v)
   i   = 1im*I_c[ind]*exp(-1im*δ)/(SM.Srated/SM.Sbase)
   i_d = real(i)
   i_q = imag(i)

   #Converstion of short-circuit time constants to open-loop time constants
   T_d0s = SM.T_ds*(SM.X_d/SM.X_ds)
   T_q0s = SM.T_qs*(SM.X_q/SM.X_qs)
   T_d0ss = SM.T_dss*(SM.X_ds/SM.X_dss)
   T_q0ss = SM.T_qss*(SM.X_qs/SM.X_qss)

   #Auxillary variables
   γ_d = T_d0ss * SM.X_dss * (SM.X_d - SM.X_ds) / (T_d0s * SM.X_ds)
   γ_q = T_q0ss * SM.X_qss * (SM.X_q - SM.X_qs) / (T_q0s * SM.X_qs)

   e_qss = v_q + SM.R_a * i_q + SM.X_dss * i_d
   e_dss = v_d + SM.R_a * i_d - SM.X_qss * i_q

   e_ds =  (SM.X_q - SM.X_qs - γ_q) * i_q
   #-(SM.X_qs - SM.X_qss + γ_q) * i_q + e_dss #also valid

   e_qs = 0.
   v_f =  1.
   if SM.T_AA != 0.
      e_qs = (SM.T_AA/T_d0s) * (-(SM.X_d - SM.X_ds - γ_d) *i_d + (T_d0s/SM.T_AA - 1)*(e_qss + (SM.X_ds - SM.X_dss + γ_d) *i_d ))
      v_f  = (T_d0s / SM.T_AA) * (e_qss - e_qs + (SM.X_ds - SM.X_dss + γ_d) *i_d)
   else
      e_qs = e_qss + (SM.X_ds - SM.X_dss + γ_d) * i_d
      v_f  = e_qs  + (SM.X_d  - SM.X_ds  - γ_d) * i_d
   end
   #AVR & OEL
   #Field current in a non-reciprocal system, otherwise would be: ifd = (E_f - T_d0s * de_qs) / (X_d - X-l)
   ifd = v_f  #due to: ifd = (E_f - T_d0s * de_qs)
   Vref = abs(v) + v_f/SM.G1

   #Pm also needs to be initialized
   Pm = (v_q + SM.R_a * i_q) * i_q + (v_d + SM.R_a * i_d) * i_d

   #Create new bus
   node_temp = SixOrderMarconatoMachineAVROEL(Sbase=SM.Sbase,Srated=SM.Srated,H=SM.H, P=Pm, D=SM.D, Ω=SM.Ω, R_a=SM.R_a, T_ds=SM.T_ds, T_qs=SM.T_qs,
                                        T_dss=SM.T_dss, T_qss=SM.T_qss, X_d=SM.X_d, X_q=SM.X_q, X_ds=SM.X_ds,
                                        X_qs=SM.X_qs, X_dss=SM.X_dss, X_qss=SM.X_qss, T_AA=SM.T_AA, V0 = Vref,
                                        Ifdlim = SM.Ifdlim, L1 = SM.L1, G1 = SM.G1, Ta = SM.Ta, Tb = SM.Tb,
                                        G2 = SM.G2, L2 = SM.L2);

   # Structure from node: u_r, u_i, θ, ω, e_ds, e_qs, e_dss,e_qss,ifd,timer,x1,E_f
   return [v_d_temp, v_q_temp, δ, 0., e_ds, e_qs, e_dss, e_qss,ifd,SM.L1,v_f,v_f], node_temp
end

function InitNode(load::Union{SimpleRecoveryLoad,SimpleRecoveryLoadParam},ind::Int64,I_c::Vector{Complex{Float64}},ic_lf::Array{Float64,1},ind_offset::Int64)
   v_d_temp = ic_lf[ind_offset]
   v_q_temp = ic_lf[ind_offset+1]
   v = sqrt(v_d_temp^2 + v_q_temp^2)
   xd = load.P0 - load.Pt * (v^2)
   xq = load.Q0 - load.Qt * (v^2)
   return [v_d_temp,v_q_temp,xd, xq], load
end

function InitNode(GFC::Union{GridFormingConverter,GridFormingConverterParam,GridFormingConverterCSA,GridFormingConverterCSAAntiWindup,GFMCurrentPrio},ind::Int64,I_c::Vector{Complex{Float64}},ic_lf::Array{Float64,1},ind_offset::Int64)
   v_d_temp = ic_lf[ind_offset]
   v_q_temp = ic_lf[ind_offset+1]
   U0 = v_d_temp+1im*v_q_temp

   #The current of the capacitor has to be related, since rf,xlf and xcf are related to Sbase!!!
   i1 = I_c[ind] / (GFC.Srated/GFC.Sbase) + U0/(-1im*GFC.xcf) /(GFC.Srated/GFC.Sbase)
   E = U0 + (GFC.rf + 1im*GFC.xlf) * i1
   θ = angle(U0)
   ω = 0.0

   s = U0 * conj(I_c[ind]) #/ (GFC.Srated/GFC.Sbase)
   p = real(s)
   q = imag(s)
   Q = q
   q0set = q

   idqmeas = I_c[ind]*(cos(-θ)+1im*sin(-θ)) / (GFC.Srated/GFC.Sbase) #1im*
   idmeas = real(idqmeas)
   iqmeas = imag(idqmeas)

   idq = i1*(cos(-θ)+1im*sin(-θ)) #1im*
   id = real(idq)
   iq = imag(idq)

   U0 = U0*(cos(-θ)+1im*sin(-θ))
   udmeas = real(U0) #should be equal to abs(U0)
   uqmeas = imag(U0) #should be zero

   E0 = E*(cos(-θ)+1im*sin(-θ))
   umd = real(E0)
   umq = imag(E0)

   e_id = (umd - udmeas + iq * GFC.xlf) #/ GFC.Ki_i
   e_iq = (umq - uqmeas - id * GFC.xlf) #/ GFC.Ki_i

   e_ud = (id - idmeas + uqmeas / GFC.xcf) #/ GFC.Ki_u #hier müsste es ohne idmeas und iqmeas sein
   e_uq = (iq - iqmeas - udmeas / GFC.xcf) #/ GFC.Ki_u #passt das überhaupt mit dem Srated/Sbase???

   if typeof(GFC) == GridFormingConverterParam
      GFC_new = GridFormingConverterParam(
         Sbase = GFC.Sbase,
         Srated = GFC.Srated,
         p0set = GFC.p0set,
         q0set = q0set, #new
         u0set = GFC.u0set,
         Kp_droop = GFC.Kp_droop,
         Kq_droop = GFC.Kq_droop,
         ωf_P = GFC.ωf_P,
         ωf_Q = GFC.ωf_Q,
         xlf = GFC.xlf,
         rf = GFC.rf,
         xcf = GFC.xcf,
         Kp_u = GFC.Kp_u,
         Ki_u = GFC.Ki_u,
         Kp_i = GFC.Kp_i,
         Ki_i = GFC.Ki_i,
         imax = GFC.imax,
         Kvi = GFC.Kvi,
         σXR = GFC.σXR,
         K_vq = GFC.K_vq,
         p_ind = GFC.p_ind
      )
      return [v_d_temp, v_q_temp,θ,ω,Q,e_ud,e_uq,e_id,e_iq,abs(idq),p], GFC_new
   elseif typeof(GFC) == GridFormingConverterCSA
      GFC_new = GridFormingConverterCSA(
         Sbase = GFC.Sbase,
         Srated = GFC.Srated,
         p0set = GFC.p0set,
         q0set = q0set, #new
         u0set = GFC.u0set,
         Kp_droop = GFC.Kp_droop,
         Kq_droop = GFC.Kq_droop,
         ωf_P = GFC.ωf_P,
         ωf_Q = GFC.ωf_Q,
         xlf = GFC.xlf,
         rf = GFC.rf,
         xcf = GFC.xcf,
         Kp_u = GFC.Kp_u,
         Ki_u = GFC.Ki_u,
         Kp_i = GFC.Kp_i,
         Ki_i = GFC.Ki_i,
         imax = GFC.imax,
         Kvi = GFC.Kvi,
         σXR = GFC.σXR,
         K_vq = GFC.K_vq,
         imax_csa = GFC.imax_csa,
         p_ind = GFC.p_ind
      )
      #,abs(E0),abs(U0/(-1im*GFC.xcf))/(GFC.Srated*GFC.Sbase),p,q
      return [v_d_temp, v_q_temp,θ,ω,Q,e_ud,e_uq,e_id,e_iq,abs(idq)], GFC_new
   elseif typeof(GFC) == GridFormingConverterCSAAntiWindup
      GFC_new = GridFormingConverterCSAAntiWindup(
         Sbase = GFC.Sbase,
         Srated = GFC.Srated,
         p0set = GFC.p0set,
         q0set = q0set, #new
         u0set = GFC.u0set,
         Kp_droop = GFC.Kp_droop,
         Kq_droop = GFC.Kq_droop,
         ωf_P = GFC.ωf_P,
         ωf_Q = GFC.ωf_Q,
         xlf = GFC.xlf,
         rf = GFC.rf,
         xcf = GFC.xcf,
         Kp_u = GFC.Kp_u,
         Ki_u = GFC.Ki_u,
         Kp_i = GFC.Kp_i,
         Ki_i = GFC.Ki_i,
         imax = GFC.imax,
         Kvi = GFC.Kvi,
         σXR = GFC.σXR,
         K_vq = GFC.K_vq,
         imax_csa = GFC.imax_csa,
         p_ind = GFC.p_ind
      )
      #,abs(E0),abs(U0/(-1im*GFC.xcf))/(GFC.Srated*GFC.Sbase),p,q
      return [v_d_temp, v_q_temp,θ,ω,Q,e_ud,e_uq,e_id,e_iq,abs(idq),p], GFC_new
   elseif typeof(GFC) == GFMCurrentPrio
      GFC_new = GFMCurrentPrio(
         Sbase = GFC.Sbase,
         Srated = GFC.Srated,
         p0set = GFC.p0set,
         q0set = q0set, #new
         u0set = GFC.u0set,
         Kp_droop = GFC.Kp_droop,
         Kq_droop = GFC.Kq_droop,
         ωf_P = GFC.ωf_P,
         ωf_Q = GFC.ωf_Q,
         xlf = GFC.xlf,
         rf = GFC.rf,
         xcf = GFC.xcf,
         Kp_u = GFC.Kp_u,
         Ki_u = GFC.Ki_u,
         Kp_i = GFC.Kp_i,
         Ki_i = GFC.Ki_i,
         imax = GFC.imax,
         Kvi = GFC.Kvi,
         σXR = GFC.σXR,
         K_vq = GFC.K_vq,
         imax_csa = GFC.imax_csa,
         iprio = GFC.iprio,
         p_ind = GFC.p_ind
      )
      #,abs(E0),abs(U0/(-1im*GFC.xcf))/(GFC.Srated*GFC.Sbase),p,q
      return [v_d_temp, v_q_temp,θ,ω,Q,e_ud,e_uq,e_id,e_iq,abs(idq),p,q], GFC_new
   else
      GFC_new = GridFormingConverter(
         Sbase = GFC.Sbase,
         Srated = GFC.Srated,
         p0set = GFC.p0set,
         q0set = q0set, #new
         u0set = GFC.u0set,
         Kp_droop = GFC.Kp_droop,
         Kq_droop = GFC.Kq_droop,
         ωf = GFC.ωf,
         xlf = GFC.xlf,
         rf = GFC.rf,
         xcf = GFC.xcf,
         Kp_u = GFC.Kp_u,
         Ki_u = GFC.Ki_u,
         Kp_i = GFC.Kp_i,
         Ki_i = GFC.Ki_i,
         imax = GFC.imax,
         Kvi = GFC.Kvi,
         σXR = GFC.σXR,
      )
      return [v_d_temp, v_q_temp,θ,ω,Q,e_ud,e_uq,e_id,e_iq,abs(E),θ/pi*180.0], GFC_new
   end
end

function InitNode(PE::oPFC,ind::Int64,I_c::Vector{Complex{Float64}},ic_lf::Array{Float64,1},ind_offset::Int64)
   v_d_temp = ic_lf[ind_offset]
   v_q_temp = ic_lf[ind_offset+1]
   U = v_d_temp+1im*v_q_temp
   s = U * conj(I_c[ind]) #/ (GFC.Srated/GFC.Sbase)
   q = imag(s)

   oPFC_new = oPFC(
      Cd = PE.Cd,
      Pdc = PE.Pdc,
      Ulow = PE.Ulow,
      Qn = q*(abs(U)^0.9), #new
      t0 = PE.t0,
      ϵ = PE.ϵ,
      p_ind = PE.p_ind,
   )

   return [v_d_temp, v_q_temp,abs(U),1.0,PE.Pdc,q,0.0,0.0], oPFC_new
end

function InitNode(MC::Union{MatchingControl,MatchingControlRed},ind::Int64,I_c::Vector{Complex{Float64}},ic_lf::Array{Float64,1},ind_offset::Int64)
   v_d_temp = ic_lf[ind_offset]
   v_q_temp = ic_lf[ind_offset+1]
   U0 = v_d_temp+1im*v_q_temp
   θ = angle(U0)

   s = U0 * conj(I_c[ind])/ (MC.Srated/MC.Sbase)
   p = real(s)
   q = imag(s)

   #The current of the capacitor has to be related, since rf,xlf and xcf are related to Sbase!!!
   i1 = I_c[ind] / (MC.Srated/MC.Sbase) + U0/(-1im*MC.xcf) #/(MC.Srated/MC.Sbase)
   E = U0 + (MC.rf + 1im*MC.xlf) * i1

   idqmeas = I_c[ind]*(cos(-θ)+1im*sin(-θ)) / (MC.Srated/MC.Sbase) #1im*
   idmeas = real(idqmeas)
   iqmeas = imag(idqmeas)

   idq = i1*(cos(-θ)+1im*sin(-θ)) #1im*
   id = real(idq)
   iq = imag(idq)

   P_before = real(conj(i1) * E)
   Q_before = imag(conj(i1) * E)
   dP =  P_before - p
   idc0 = MC.gdc + MC.p0set/ (MC.Srated/MC.Sbase) + dP
   p0_new = idc0 - MC.gdc - dP
   udc = 0.0 #ist hier nur das delta

   U0 = U0*(cos(-θ)+1im*sin(-θ))
   udmeas = real(U0) #should be equal to abs(U0)
   uqmeas = imag(U0) #should be zero

   E0 = E*(cos(-θ)+1im*sin(-θ))
   umd = real(E0)
   umq = imag(E0)

   e_id = (umd - udmeas + iq * MC.xlf ) #- id*MC.rf
   e_iq = (umq - uqmeas - id * MC.xlf )#- iq*MC.rf

   e_ud = (id - idmeas + uqmeas / MC.xcf)#/ MC.Ki_u #hier müsste es ohne idmeas und iqmeas sein
   e_uq = (iq - iqmeas - udmeas / MC.xcf) #/ MC.Ki_u #passt das überhaupt mit dem Srated/Sbase???

   if typeof(MC) == MatchingControlRed
      MC_new = MatchingControlRed(
      Sbase = MC.Sbase,
      Srated = MC.Srated,
      p0set = p0_new, #new
      u0set = MC.u0set,
      Kp_uset = MC.Kp_uset,
      Ki_uset = MC.Ki_uset,
      Kdc = MC.Kdc,
      gdc = MC.gdc,
      cdc = MC.cdc,
      xlf = MC.xlf,
      rf = MC.rf,
      xcf =  MC.xcf,
      Tdc = MC.Tdc,
      Kp_u = MC.Kp_u,
      Ki_u = MC.Ki_u,
      Kp_i = MC.Kp_i,
      Ki_i = MC.Ki_i,
      imax_csa = MC.imax_csa,
      imax_dc = MC.imax_dc,
      p_red = MC.p_red,
      LVRT_on = MC.LVRT_on,
      p_ind = MC.p_ind,
      )
  else
      MC_new = MatchingControl(
      Sbase = MC.Sbase,
      Srated = MC.Srated,
      p0set = p0_new, #new
      u0set = MC.u0set,
      Kp_uset = MC.Kp_uset,
      Ki_uset = MC.Ki_uset,
      Kdc = MC.Kdc,
      gdc = MC.gdc,
      cdc = MC.cdc,
      xlf = MC.xlf,
      rf = MC.rf,
      xcf =  MC.xcf,
      Tdc = MC.Tdc,
      Kp_u = MC.Kp_u,
      Ki_u = MC.Ki_u,
      Kp_i = MC.Kp_i,
      Ki_i = MC.Ki_i,
      imax_csa = MC.imax_csa,
      imax_dc = MC.imax_dc,
      p_red = MC.p_red,
      p_ind = MC.p_ind,
      )
  end

    return [v_d_temp, v_q_temp,θ,udc,idc0,abs(U0),e_ud,e_uq,e_id,e_iq,dP,abs(i1),0.0], MC_new #P_before,Q_before,p,q,abs(idq)
end

function InitNode(VOC::dVOC,ind::Int64,I_c::Vector{Complex{Float64}},ic_lf::Array{Float64,1},ind_offset::Int64)
   v_d_temp = ic_lf[ind_offset]
   v_q_temp = ic_lf[ind_offset+1]
   U0 = v_d_temp+1im*v_q_temp
   θ = angle(U0)

   s = U0 * conj(I_c[ind]) / (VOC.Srated/VOC.Sbase)
   p = real(s)
   q = imag(s)

   #The current of the capacitor has to be related, since rf,xlf and xcf are related to Sbase!!!
   i1 = I_c[ind] / (VOC.Srated/VOC.Sbase) + U0/(-1im*VOC.xcf) #/(VOC.Srated/VOC.Sbase)
   E = U0 + (VOC.rf + 1im*VOC.xlf) * i1

   idqmeas = I_c[ind]*(cos(-θ)+1im*sin(-θ)) / (VOC.Srated/VOC.Sbase) #1im*
   idmeas = real(idqmeas)
   iqmeas = imag(idqmeas)

   idq = i1*(cos(-θ)+1im*sin(-θ)) #1im*
   id = real(idq)
   iq = imag(idq)

   P_before = real(conj(i1) * E)
   Q_before = imag(conj(i1) * E)
   dP =  P_before - p
   idc0 = VOC.gdc + VOC.p0set/ (VOC.Srated/VOC.Sbase)  + dP
   p0_new = idc0 - VOC.gdc - dP
   udc = 0.0 #ist hier nur das delta

   U0 = U0*(cos(-θ)+1im*sin(-θ))
   udmeas = real(U0) #should be equal to abs(U0)
   uqmeas = imag(U0) #should be zero

   qset = VOC.q0set
   x1 = (qset + VOC.u0set^2) * VOC.alpha/ (2 * VOC.alpha) + sqrt(((-qset - VOC.u0set^2*VOC.alpha) / (2 * VOC.alpha))^2 - qset * VOC.u0set / VOC.alpha)
   vd_int = sqrt(x1) - VOC.ϵ
   E0 = E*(cos(-θ)+1im*sin(-θ))
   umd = real(E0)
   umq = imag(E0)


   e_id = (umd - udmeas + iq * VOC.xlf ) #- id*VOC.rf
   e_iq = (umq - uqmeas - id * VOC.xlf )#- iq*VOC.rf

   e_ud = (id - idmeas + uqmeas / VOC.xcf)#/ VOC.Ki_u #hier müsste es ohne idmeas und iqmeas sein
   e_uq = (iq - iqmeas - udmeas / VOC.xcf) #/ VOC.Ki_u #passt das überhaupt mit dem Srated/Sbase???

   VOC_new = dVOC(
          Sbase = VOC.Sbase,
          Srated = VOC.Srated,
          p0set = p0_new, #new
          q0set = q,
          u0set = VOC.u0set,
          eta = VOC.eta,
          alpha = VOC.alpha,
          Kdc = VOC.Kdc,
          gdc = VOC.gdc,
          cdc = VOC.cdc,
          xlf = VOC.xlf,
          rf = VOC.rf,
          xcf =  VOC.xcf,
          Tdc = VOC.Tdc,
          Kp_u = VOC.Kp_u,
          Ki_u = VOC.Ki_u,
          Kp_i = VOC.Kp_i,
          Ki_i = VOC.Ki_i,
          imax_csa = VOC.imax_csa,
          imax_dc = VOC.imax_dc,
          p_red = VOC.p_red,
          ϵ = VOC.ϵ,
          LVRT_on = VOC.LVRT_on,
          p_ind = VOC.p_ind,
          )

    return [v_d_temp, v_q_temp,θ,udc,idc0,vd_int,e_ud,e_uq,e_id,e_iq,p,q,dP,abs(i1),0.0,0.0], VOC_new #,idmeas,iqmeas,id,iq
end

function InitNode(DR::droop,ind::Int64,I_c::Vector{Complex{Float64}},ic_lf::Array{Float64,1},ind_offset::Int64)
   v_d_temp = ic_lf[ind_offset]
   v_q_temp = ic_lf[ind_offset+1]
   U0 = v_d_temp+1im*v_q_temp
   θ = angle(U0)

   s = U0 * conj(I_c[ind]) / (DR.Srated/DR.Sbase)
   p = real(s)
   q = imag(s)

   #The current of the capacitor has to be related, since rf,xlf and xcf are related to Sbase!!!
   i1 = I_c[ind] / (DR.Srated/DR.Sbase) + U0/(-1im*DR.xcf) #/(DR.Srated/DR.Sbase)
   E = U0 + (DR.rf + 1im*DR.xlf) * i1

   idqmeas = I_c[ind]*(cos(-θ)+1im*sin(-θ)) / (DR.Srated/DR.Sbase) #1im*
   idmeas = real(idqmeas)
   iqmeas = imag(idqmeas)

   idq = i1*(cos(-θ)+1im*sin(-θ)) #1im*
   id = real(idq)
   iq = imag(idq)

   P_before = real(conj(i1) * E)
   Q_before = imag(conj(i1) * E)
   dP =  P_before - p
   idc0 = DR.gdc + DR.p0set/ (DR.Srated/DR.Sbase)  + dP
   p0_new = idc0 - DR.gdc - dP
   udc = 0.0 #ist hier nur das delta

   U0 = U0*(cos(-θ)+1im*sin(-θ))
   udmeas = real(U0) #should be equal to abs(U0)
   uqmeas = imag(U0) #should be zero

   E0 = E*(cos(-θ)+1im*sin(-θ))
   umd = real(E0)
   umq = imag(E0)

   e_id = (umd - udmeas + iq * DR.xlf ) #- id*DR.rf
   e_iq = (umq - uqmeas - id * DR.xlf )#- iq*DR.rf

   e_ud = (id - idmeas + uqmeas / DR.xcf)#/ DR.Ki_u #hier müsste es ohne idmeas und iqmeas sein
   e_uq = (iq - iqmeas - udmeas / DR.xcf) #/ DR.Ki_u #passt das überhaupt mit dem Srated/Sbase???

   droop_new = droop(
          Sbase = DR.Sbase,
          Srated = DR.Srated,
          p0set = p0_new, #new
          u0set = DR.u0set,
          Kp_droop = DR.Kp_droop,
          Kp_uset = DR.Kp_uset,
          Ki_uset = DR.Ki_uset,
          Kdc = DR.Kdc,
          gdc = DR.gdc,
          cdc = DR.cdc,
          xlf = DR.xlf,
          rf = DR.rf,
          xcf =  DR.xcf,
          Tdc = DR.Tdc,
          Kp_u = DR.Kp_u,
          Ki_u = DR.Ki_u,
          Kp_i = DR.Kp_i,
          Ki_i = DR.Ki_i,
          imax_csa = DR.imax_csa,
          imax_dc = DR.imax_dc,
          p_red = DR.p_red,
          LVRT_on = DR.LVRT_on,
          p_ind = DR.p_ind,
          )

    return [v_d_temp, v_q_temp,θ,udc,idc0,abs(U0),e_ud,e_uq,e_id,e_iq,p,dP,abs(i1),abs(i1),0.0], droop_new #,idc0,abs(E0),P_before,Q_before,p,q
end

function InitNode(VSM0::VSM,ind::Int64,I_c::Vector{Complex{Float64}},ic_lf::Array{Float64,1},ind_offset::Int64)
   v_d_temp = ic_lf[ind_offset]
   v_q_temp = ic_lf[ind_offset+1]
   U0 = v_d_temp+1im*v_q_temp
   θ = angle(U0)

   s = U0 * conj(I_c[ind]) / (VSM0.Srated/VSM0.Sbase)
   p = real(s)
   q = imag(s)

   #The current of the capacitor has to be related, since rf,xlf and xcf are related to Sbase!!!
   i1 = I_c[ind] / (VSM0.Srated/VSM0.Sbase) + U0/(-1im*VSM0.xcf) #/(VSM0.Srated/VSM0.Sbase)
   E = U0 + (VSM0.rf + 1im*VSM0.xlf) * i1

   idqmeas = I_c[ind]*(cos(-θ)+1im*sin(-θ)) / (VSM0.Srated/VSM0.Sbase) #1im*
   idmeas = real(idqmeas)
   iqmeas = imag(idqmeas)

   idq = i1*(cos(-θ)+1im*sin(-θ)) #1im*
   id = real(idq)
   iq = imag(idq)

   P_before = real(conj(i1) * E)
   Q_before = imag(conj(i1) * E)
   
   dP =  P_before - p
   idc0 = VSM0.gdc + VSM0.p0set/ (VSM0.Srated/VSM0.Sbase) + dP
   p0_new = idc0 - VSM0.gdc - dP
   udc = 0.0 #ist hier nur das delta

   U0 = U0*(cos(-θ)+1im*sin(-θ))
   udmeas = real(U0) #should be equal to abs(U0)
   uqmeas = imag(U0) #should be zero

   E0 = E*(cos(-θ)+1im*sin(-θ))
   umd = real(E0)
   umq = imag(E0)

   e_id = (umd - udmeas + iq * VSM0.xlf ) #- id*VSM0.rf
   e_iq = (umq - uqmeas - id * VSM0.xlf )#- iq*VSM0.rf

   e_ud = (id - idmeas + uqmeas / VSM0.xcf)#/ VSM0.Ki_u #hier müsste es ohne idmeas und iqmeas sein
   e_uq = (iq - iqmeas - udmeas / VSM0.xcf) #/ VSM0.Ki_u #passt das überhaupt mit dem Srated/Sbase???

   vsm_new = VSM(
          Sbase = VSM0.Sbase,
          Srated = VSM0.Srated,
          p0set = p0_new, #new
          u0set = VSM0.u0set,
          J = VSM0.J,
          Dp = VSM0.Dp,
          Kp_uset = VSM0.Kp_uset,
          Ki_uset = VSM0.Ki_uset,
          Kdc = VSM0.Kdc,
          gdc = VSM0.gdc,
          cdc = VSM0.cdc,
          xlf = VSM0.xlf,
          rf = VSM0.rf,
          xcf =  VSM0.xcf,
          Tdc = VSM0.Tdc,
          Kp_u = VSM0.Kp_u,
          Ki_u = VSM0.Ki_u,
          Kp_i = VSM0.Kp_i,
          Ki_i = VSM0.Ki_i,
          imax_csa = VSM0.imax_csa,
          imax_dc = VSM0.imax_dc,
          p_red = VSM0.p_red,
          LVRT_on = VSM0.LVRT_on,
          p_ind = VSM0.p_ind,
          )

    return [v_d_temp, v_q_temp,θ,0.0,udc,idc0,abs(U0),e_ud,e_uq,e_id,e_iq,p,dP,abs(i1),0.0], vsm_new #,idmeas,iqmeas,id,iq
end

##gentpj
function InitNode(GP::gentpj,ind::Int64,I_c::Vector{Complex{Float64}},ic_lf::Array{Float64,1},ind_offset::Int64)

   #saturation function
   function sat_scaled_quadratic(x1)
       A = (1.2 - sqrt(1.2*GP.S_12/GP.S_10)) / (1.0 - sqrt(1.2*GP.S_12/GP.S_10))
       B = GP.S_10 / (1 - A)^2
       if x1 > A
           return B*(x1-A)^2/x1
       else
           return 0.0
       end
   end
   function sat_exponential(x1)
      q = log(GP.S_12/GP.S_10)/log(1.2)
      return GP.S_10*x1^q
  end
   
   #get inital values
   v_r_term = ic_lf[ind_offset]
   v_i_term = ic_lf[ind_offset+1]
   i_r = real(I_c[ind])/(GP.Srated/GP.Sbase)
   i_i = imag(I_c[ind])/(GP.Srated/GP.Sbase)

   #Calculate saturation
   E_l = norm([v_i_term + i_i * GP.R_a + i_r * GP.X_l,v_r_term + i_r * GP.R_a - i_i * GP.X_l]) #magnitude not changed by dq-transformation
   S_d = sat_exponential(E_l)
   S_q = (GP.X_q)/(GP.X_d) * S_d

   #Calulate Polradspannung
   Xdsatss = (GP.X_dss - GP.X_l) /  (1 + S_d) + GP.X_l
   Xqsatss = (GP.X_qss - GP.X_l) /  (1 + S_q) + GP.X_l   
   Vr = v_r_term + GP.R_a *i_r  - Xqsatss * i_i
   Vi = v_i_term + GP.R_a *i_i  + Xdsatss * i_r

   #calculate δ
   K_temp = (GP.X_q-GP.X_l)/(1 + S_q) + GP.X_l
   δ = atan((v_i_term+i_r*K_temp + GP.R_a *i_i)/(v_r_term - i_i*K_temp +GP.R_a *i_r))

   #terminal voltage transformation to local dq-system
   v = v_r_term+1im*v_i_term
   v = 1im*v*exp(-1im*δ)
   v_d = real(v)
   v_q = imag(v)

   #current transformation to dq-system
   i_c = 1im*I_c[ind]*exp(-1im*δ)/(GP.Srated/GP.Sbase)
   i_d = real(i_c)
   i_q = imag(i_c)

   #initial conditions as derivative has to be 0
   e_q2 = 0;
   e_d2 = 0;
   e_d1 = 0;

   #equations 
   e_q1 = v_q + i_q * GP.R_a + i_d * ((GP.X_d - GP.X_l) / (1 + S_d) + GP.X_l)
   e_qss = e_q1 + e_q2 - i_d * ((GP.X_d - GP.X_dss) / (1 + S_d))  

   e_dss = e_d1 + e_d2 + i_q * ((GP.X_q - GP.X_qss) / (1 + S_q))

   e_qs = e_qss + i_d * ((GP.X_ds - GP.X_dss) / (1 + S_d))
   e_ds = e_dss - i_q * ((GP.X_qs - GP.X_qss) / (1 + S_q))

   #field voltage
   v_f = (1+S_d)*e_q1

   #Pm also needs to be initialized
   Pm = (v_q + GP.R_a * i_q) * i_q + (v_d + GP.R_a * i_d) * i_d

   #Create new bus
   node_temp = gentpj(Sbase=GP.Sbase,Srated=GP.Srated,H=GP.H, P=Pm, D=GP.D, Ω=GP.Ω,
                                           E_fd=v_f, R_a=GP.R_a, T_d0s=GP.T_d0s, T_q0s=GP.T_q0s, T_d0ss=GP.T_d0ss,
                                           T_q0ss=GP.T_q0ss, X_d=GP.X_d, X_q=GP.X_q, X_ds=GP.X_ds, X_qs=GP.X_qs,
                                           X_dss=GP.X_dss, X_qss=GP.X_qss, X_l=GP.X_l, S_10=GP.S_10, S_12=GP.S_12, K_is = GP.K_is);
return [v_r_term, v_i_term, δ, 0., e_ds, e_qs, e_dss, e_qss,E_l], node_temp
end

function InitNode(GP::gentpjAVROEL,ind::Int64,I_c::Vector{Complex{Float64}},ic_lf::Array{Float64,1},ind_offset::Int64)

   #saturation function
   function sat_scaled_quadratic(x1)
       A = (1.2 - sqrt(1.2*GP.S_12/GP.S_10)) / (1.0 - sqrt(1.2*GP.S_12/GP.S_10))
       B = GP.S_10 / (1 - A)^2
       if x1 > A
           return B*(x1-A)^2/x1
       else
           return 0.0
       end
   end
   function sat_exponential(x1)
      q = log(GP.S_12/GP.S_10)/log(1.2)
      return GP.S_10*x1^q
  end
   
   #get inital values
   v_r_term = ic_lf[ind_offset]
   v_i_term = ic_lf[ind_offset+1]
   i_r = real(I_c[ind])/(GP.Srated/GP.Sbase)
   i_i = imag(I_c[ind])/(GP.Srated/GP.Sbase)

   #Calculate saturation
   E_l = norm([v_i_term + i_i * GP.R_a + i_r * GP.X_l,v_r_term + i_r * GP.R_a - i_i * GP.X_l]) #magnitude not changed by dq-transformation
   S_d = sat_exponential(E_l)
   S_q = (GP.X_q)/(GP.X_d) * S_d

   #Calulate Polradspannung
   Xdsatss = (GP.X_dss - GP.X_l) /  (1 + S_d) + GP.X_l
   Xqsatss = (GP.X_qss - GP.X_l) /  (1 + S_q) + GP.X_l   
   Vr = v_r_term + GP.R_a *i_r  - Xqsatss * i_i
   Vi = v_i_term + GP.R_a *i_i  + Xdsatss * i_r

   #calculate δ
   K_temp = (GP.X_q-GP.X_l)/(1 + S_q) + GP.X_l
   δ = atan((v_i_term+i_r*K_temp + GP.R_a *i_i)/(v_r_term - i_i*K_temp +GP.R_a *i_r))

   #terminal voltage transformation to local dq-system
   v = v_r_term+1im*v_i_term
   v = 1im*v*exp(-1im*δ)
   v_d = real(v)
   v_q = imag(v)

   #current transformation to dq-system
   i_c = 1im*I_c[ind]*exp(-1im*δ)/(GP.Srated/GP.Sbase)
   i_d = real(i_c)
   i_q = imag(i_c)

   #initial conditions as derivative has to be 0
   e_q2 = 0;
   e_d2 = 0;
   e_d1 = 0;

   #equations
   edqss = 1im*(Vr+1im*Vi)*exp(-1im*δ)
   
   e_q1 = v_q + i_q * GP.R_a + i_d * ((GP.X_d - GP.X_l) / (1 + S_d) + GP.X_l)
   e_qss = e_q1 + e_q2 - i_d * ((GP.X_d - GP.X_dss) / (1 + S_d))
  
   e_dss = e_d1 + e_d2 + i_q * ((GP.X_q - GP.X_qss) / (1 + S_q))

   e_qs = e_qss + i_d * ((GP.X_ds - GP.X_dss) / (1 + S_d))
   e_ds = e_dss - i_q * ((GP.X_qs - GP.X_qss) / (1 + S_q))

   #field voltage
   v_f = (1+S_d)*e_q1

   #AVR & OEL
   #Field current in a non-reciprocal system, otherwise would be: ifd = (E_f - T_d0s * de_qs) / (X_d - X-l)
   ifd = v_f  #due to: ifd = (E_f - T_d0s * de_qs) de_qs = (E_f - (1 + S_d) * e_q1) * 1 / T_d0s
   Vref = abs(v) + v_f/GP.G1

   #Pm also needs to be initialized
   Pm = (v_q + GP.R_a * i_q) * i_q + (v_d + GP.R_a * i_d) * i_d

   #Create new bus
   node_temp = gentpjAVROEL(Sbase=GP.Sbase,Srated=GP.Srated,H=GP.H, P=Pm, D=GP.D, Ω=GP.Ω,
                                           R_a=GP.R_a, T_d0s=GP.T_d0s, T_q0s=GP.T_q0s, T_d0ss=GP.T_d0ss,
                                           T_q0ss=GP.T_q0ss, X_d=GP.X_d, X_q=GP.X_q, X_ds=GP.X_ds, X_qs=GP.X_qs,
                                           X_dss=GP.X_dss, X_qss=GP.X_qss, X_l=GP.X_l, S_10=GP.S_10, S_12=GP.S_12, K_is = GP.K_is,
                                           V0 = Vref, Ifdlim = GP.Ifdlim, L1 = GP.L1, G1 = GP.G1, Ta = GP.Ta, Tb = GP.Tb,
                                           G2 = GP.G2, L2 = GP.L2);
return [v_r_term, v_i_term, δ, 0., e_ds, e_qs, e_dss, e_qss,E_l,ifd,GP.L1,v_f,v_f], node_temp
end

function InitNode(F::ThreePhaseFaultContinouos,ind::Int64,I_c::Vector{Complex{Float64}},ic_lf::Array{Float64,1},ind_offset::Int64)
   v_d_temp = ic_lf[ind_offset]
   v_q_temp = ic_lf[ind_offset+1]
   return [v_d_temp, v_q_temp,F.rfault,F.xfault], F
end