function c = aeroCoefVec(P, M, alphaTot)
%AEROCOEFVEC Współczynniki aerodynamiczne w chwili (ℳ, α_tot) – wersja szybka (wektor).
%   c = [CD CLa CMa CMq CNq Cmagf Cmagm Cspin CNa CPN clamp] – znaczenie jak w aeroCoef.
%   Korzysta z macierzy P.coefTab przygotowanej w makeParams (kolumny: CDo CD2 CLa CLa3 Cma Cma3
%   Cma5 CNa CNq Cmq Cmq2 CYpa Clp CPN).
mach = P.mach; n = numel(mach);
clamp = 0; Mc = M;
if ~(Mc >= mach(1)), Mc = mach(1); clamp = -1;          % obsługuje też NaN
elseif Mc > mach(n), Mc = mach(n); clamp = 1;
end
i0 = 1 + sum(mach(2:n - 1) < Mc); i1 = i0 + 1;
w = (Mc - mach(i0)) / (mach(i1) - mach(i0));
t = P.coefTab(i0, :) + w * (P.coefTab(i1, :) - P.coefTab(i0, :));
s = sin(alphaTot); s2 = s * s; s4 = s2 * s2;
kp = P.kp; dX = P.dXcg; f = P.fv;
CNa = t(8); Cma = t(5); CNq = kp * t(9);
CMq = kp * (t(10) + t(11) * s2);
if P.fullCmqTransfer, CMq = CMq + dX * CNq - dX * Cma - dX * dX * CNa; end   % pełny transfer (normowanie qd/V)
% Cnpa(ℳ, α_tot): interpolacja biliniowa, α obcięte do zakresu tabeli
al = P.cnpaAlpha; na = numel(al); cn = P.cnpaData;
aD = min(max(alphaTot * (180 / pi), al(1)), al(na));
ja = 1 + sum(al(2:na - 1) < aD);
wa = (aD - al(ja)) / (al(ja + 1) - al(ja));
cnpa = (1 - w) * ((1 - wa) * cn(i0, ja) + wa * cn(i0, ja + 1)) + w * ((1 - wa) * cn(i1, ja) + wa * cn(i1, ja + 1));
c = [f(1) * t(1) + f(2) * t(2) * s2, ...                  % C_D (B.36), C_Dα⁴ = 0
     f(3) * (t(3) + t(4) * s2), ...                       % C_Lα
     f(4) * (Cma + t(6) * s2 + t(7) * s4 + CNa * dX), ... % C_Mα
     f(5) * CMq, f(6) * CNq, ...                          % C_Mq + C_Mα̇, C_Nq + C_Nα̇
     f(7) * kp * t(12), ...                               % C_mag-f
     f(8) * kp * (cnpa + t(12) * dX), ...                 % C_mag-m
     f(9) * kp * t(13), ...                               % C_spin
     CNa, t(14), clamp];
if P.anyFd
  A = P.A; fd = P.fd;
  c(1) = c(1) + fd.CD0 * sigmaAt(A, 'CX', Mc);
  c(2) = c(2) + fd.CLa * sigmaAt(A, 'CNa', Mc);
  c(3) = c(3) + fd.CMa * sigmaAt(A, 'Cma', Mc);
  c(4) = c(4) + fd.CMq * kp * sigmaAt(A, 'Cmq', Mc);
  c(7) = c(7) + fd.Cmagm * kp * sigmaCnpaAt(A, Mc, alphaTot * (180 / pi));
  c(8) = c(8) + fd.Cspin * kp * sigmaAt(A, 'Clp', Mc);
end
end
