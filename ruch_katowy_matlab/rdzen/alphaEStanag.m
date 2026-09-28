function ae = alphaEStanag(P, M, rho, p, VF, UF)
%ALPHAESTANAG Kąt równowagi wg STANAG 4355 (oznaczenie ᾱ_e).
%   ᾱ_e = −8·I_x·p·(V̄ × dŪ/dt)/(π·ρ·d³·C_Mα·V⁴), z przyspieszeniem jak w MPM:
%   siły aerodynamiczne liczone przy odchyleniu ᾱ_e (iteracja), + ḡ + Ā_C.
%   VF, UF – wektory w układzie ogniowym; wynik: wektor 1×3 w układzie ogniowym [rad].
VF = VF(:).'; UF = UF(:).';
V = norm(VF); Vh = VF / V;
Ac = [0 0 0];
if P.coriolis, Ac = -2 * cross(P.LamF, UF); end
qS = 0.5 * rho * P.S; d = P.d;
ae = [0 0 0];
for it = 1:4
  a = norm(ae);
  if a > 0, xe = Vh * cos(a) + ae / a * sin(a); else, xe = Vh; end
  c = aeroCoefVec(P, M, a);
  L3 = cross(VF, cross(xe, VF)); Vxe = cross(VF, xe);
  ud = (-qS * c(1) * V * VF + qS * c(2) * L3 + qS * p * d * c(6) * Vxe) / P.m + [0, -P.g, 0] + Ac;
  kk = -8 * P.Ix * p / (pi * rho * d^3 * c(3) * V^4);
  ae = kk * cross(VF, ud);
end
end
