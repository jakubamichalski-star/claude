function [T, p, rho, a] = atmAt(P, hAbs)
%ATMAT Atmosfera ISO 2533 z odchyłkami warstw METCM: T = T_wz + ΔT, ρ = ρ_wz·(1 + Δρ), a = √(γRT), p = ρRT.
[T, p, rho, a] = atmosphere(hAbs);
if ~isempty(P.met)
  R_AIR = 287.05287; GAMMA_AIR = 1.4;
  hr = hAbs - P.met.hMDP;
  T = max(150, T + pchipVal(P.met.fT, hr));
  rho = rho * (1 + pchipVal(P.met.fR, hr) / 100);
  p = rho * R_AIR * T; a = sqrt(GAMMA_AIR * R_AIR * T);
end
end
