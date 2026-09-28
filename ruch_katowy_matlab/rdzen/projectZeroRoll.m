function [y, rollCorr] = projectZeroRoll(y)
%PROJECTZEROROLL Po kroku: normowanie λ̄ i przywrócenie φ = 0 układu całkowania (obrót składowych Ū, ϖ̄).
%   rollCorr – |kąt przechylenia| usunięty z układu całkowania [rad].
T = quatToT(y(10:13));
[psi, theta, phi] = eulerFromT(T);
rollCorr = abs(phi);
if rollCorr > 0
  c = cos(phi); s = sin(phi);
  a = y(5); b = y(6); y(5) = c * a - s * b; y(6) = s * a + c * b;
  a = y(8); b = y(9); y(8) = c * a - s * b; y(9) = s * a + c * b;
  y(14) = y(14) + phi;
end
y(10:13) = quatFromEuler(psi, theta, 0);
end
