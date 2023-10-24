Fx, Fy, Gx, Gy = GetSymbolicFactorizedJacobian(odesys);
Fxf, Fyf, Gxf, Gyf = [
  substitute(f, (params,)) for
  f in [Fx, Fy, Gx, Gy]
];
Af = Fxf - Fyf * inv(Gyf) * Gxf;
Af = substitute.(Symbolics.value.(Af),([],)) 