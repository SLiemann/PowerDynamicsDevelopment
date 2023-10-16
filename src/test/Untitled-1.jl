

eqs, aeqs, D_states, A_states = GetSymbolicEquationsAndStates(odesys);


Fx = Array{Num}(undef, 1, size(states(odesys))[1])
Diff_states = Differential.(states(odesys))
for (ind, val) in enumerate(Diff_states)
  Fx[ind] = Num.(expand_derivatives.(map(val, h)))
end

Fx*(my_rhs.(equations(odesys)))