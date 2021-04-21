using PowerDynamics
using PowerDynamics: _find_operationpoint_rootfind
using PowerDynamics: rhs, SixOrderMarconatoMachine, symbolsof
using NLsolve: nlsolve, converged
using IfElse
#include("PowerFlow.jl") # for NodalAdmittanceMatrice

function InitializeInternalDynamics(pg::PowerGrid,I_c::Array{Complex{Float64},2},ic_lf::Array{Float64,1})
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
function InitNode(SM::FourthOrderEq,ind::Int64,I_c::Array{Complex{Float64},2},ic_lf::Array{Float64,1},ind_offset::Int64)
   v_d_temp = ic_lf[ind_offset]
   v_q_temp = ic_lf[ind_offset+1]
   #Rotor angle
   δ = angle(v_d_temp+1im*v_q_temp+(+1im*SM.X_q)*I_c[ind]/(SM.S_r/Sbase))

   v   = v_d_temp +1im*v_q_temp
   v   = 1im*v*exp(-1im*δ)
   v_d = real(v)
   v_q = imag(v)
   i   = 1im*I_c[ind]*exp(-1im*δ)/(SM.S_r/Sbase)
   i_d = real(i)
   i_q = imag(i)

   V_f = v_q + (SM.X_d - SM.X_d_dash) * i_d
   pe = SM.P + (SM.X_q_dash - SM.X_d_dash)*i_d*i_q
   node_temp = FourthOrderEq(H=SM.H, P=pe, D=SM.D, Ω=SM.Ω, E_f=V_f, T_d_dash=SM.T_d_dash ,T_q_dash=SM.T_q_dash ,X_q_dash =SM.X_q_dash ,X_d_dash=SM.X_d_dash,X_d=SM.X_d, X_q=SM.X_q)
   return [v_d_temp, v_q_temp, δ, 0.], node_temp
end

InitNode(SM::SixOrderMarconatoMachineSin,ind::Int64,I_c::Array{Complex{Float64},2},ic_lf::Array{Float64,1},ind_offset::Int64) = InitNodeSM(SM,ind,I_c,ic_lf,ind_offset)
InitNode(SM::SixOrderMarconatoMachine,ind::Int64,I_c::Array{Complex{Float64},2},ic_lf::Array{Float64,1},ind_offset::Int64)    = InitNodeSM(SM,ind,I_c,ic_lf,ind_offset)
function InitNodeSM(SM,ind::Int64,I_c::Array{Complex{Float64},2},ic_lf::Array{Float64,1},ind_offset::Int64)
   v_d_temp = ic_lf[ind_offset]
   v_q_temp = ic_lf[ind_offset+1]
   #Rotor angle
   δ = angle(v_d_temp+1im*v_q_temp+(SM.R_a+1im*SM.X_q)*I_c[ind]) #+angle(v_d_temp + 1im*v_q_temp)

   v = v_d_temp +1im*v_q_temp
   v = 1im*v*exp(-1im*δ)
   v_d = real(v)
   v_q = imag(v)
   i   = 1im*I_c[ind]*exp(-1im*δ)
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
   node_temp = typeof(SM)== SixOrderMarconatoMachine ? SixOrderMarconatoMachine : SixOrderMarconatoMachineSin
   if typeof(SM) == SixOrderMarconatoMachine
      node_temp = SixOrderMarconatoMachine(H=SM.H, P=Pm, D=SM.D, Ω=SM.Ω, E_f=v_f, R_a=SM.R_a, T_ds=SM.T_ds, T_qs=SM.T_qs, T_dss=SM.T_dss, T_qss=SM.T_qss, X_d=SM.X_d, X_q=SM.X_q, X_ds=SM.X_ds, X_qs=SM.X_qs, X_dss=SM.X_dss, X_qss=SM.X_qss, T_AA=SM.T_AA);
   elseif typeof(SM) == SixOrderMarconatoMachineSin
      node_temp = SixOrderMarconatoMachineSin(H=SM.H, P=Pm, D=SM.D, Ω=SM.Ω, E_f=v_f, R_a=SM.R_a, T_ds=SM.T_ds, T_qs=SM.T_qs, T_dss=SM.T_dss, T_qss=SM.T_qss, X_d=SM.X_d, X_q=SM.X_q, X_ds=SM.X_ds, X_qs=SM.X_qs, X_dss=SM.X_dss, X_qss=SM.X_qss, T_AA=SM.T_AA);
   end
   # Structure from node: u_r, u_i, θ, ω, e_ds, e_qs, e_dss,e_qss
   return [v_d_temp, v_q_temp, δ, 0., e_ds, e_qs, e_dss, e_qss], node_temp
end

function InitNode(SM::SixOrderMarconatoMachineAVROEL,ind::Int64,I_c::Array{Complex{Float64},2},ic_lf::Array{Float64,1},ind_offset::Int64)
   v_d_temp = ic_lf[ind_offset]
   v_q_temp = ic_lf[ind_offset+1]#*15/380
   #Rotor angle
   δ = angle(v_d_temp+1im*v_q_temp+(SM.R_a+1im*SM.X_q)*I_c[ind]/(SM.Sr/SM.Sb))

   v = v_d_temp +1im*v_q_temp
   v = 1im*v*exp(-1im*δ)
   v_d = real(v)
   v_q = imag(v)
   i   = 1im*I_c[ind]*exp(-1im*δ)/6.0
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
   node_temp = SixOrderMarconatoMachineAVROEL(Sb=SM.Sb,Sr=SM.Sr,H=SM.H, P=Pm, D=SM.D, Ω=SM.Ω, R_a=SM.R_a, T_ds=SM.T_ds, T_qs=SM.T_qs,
                                        T_dss=SM.T_dss, T_qss=SM.T_qss, X_d=SM.X_d, X_q=SM.X_q, X_ds=SM.X_ds,
                                        X_qs=SM.X_qs, X_dss=SM.X_dss, X_qss=SM.X_qss, T_AA=SM.T_AA, V0 = Vref,
                                        Ifdlim = SM.Ifdlim, L1 = SM.L1, G1 = SM.G1, Ta = SM.Ta, Tb = SM.Tb,
                                        G2 = SM.G2, L2 = SM.L2);

   # Structure from node: u_r, u_i, θ, ω, e_ds, e_qs, e_dss,e_qss,ifd,timer,x1,E_f
   return [v_d_temp, v_q_temp, δ, 0., e_ds, e_qs, e_dss, e_qss,ifd,SM.L1,v_f,v_f], node_temp
end
