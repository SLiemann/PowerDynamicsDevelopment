using ModelingToolkit

@variables i_d i_q i_ld i_lq v_fd v_fq v_dref v_qref ζ_id ζ_iq ζ_vd ζ_vq  dζ_id dζ_iq dζ_vd dζ_vq 
@parameters K_pv K_iv K_ii K_pi xcf xlf rf 

i_cdref = i_ld - v_fq / xcf + K_pv * (v_dref*0 - v_fd) + ζ_vd * K_iv
i_cqref = i_lq + v_fd / xcf + K_pv * (v_qref*0 - v_fq) + ζ_vq * K_iv

v_cd = v_fd - i_q * xlf + K_pi * (i_cdref - i_d) + ζ_id * K_ii
v_cq = v_fq + i_d * xlf + K_pi * (i_cqref - i_q) + ζ_iq * K_ii

di_d = v_cd - i_d * rf + i_q * xlf - v_fd
di_q = v_cq - i_q * rf - i_d * xlf - v_fq