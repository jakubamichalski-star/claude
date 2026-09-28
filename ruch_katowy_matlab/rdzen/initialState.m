function y = initialState(P)
%INITIALSTATE Stan początkowy: oś odchylona od lufy o α_tot0 w kierunku φ_A0, Ū wzdłuż lufy.
%   Stan y (15×1): X₁ X₂ X₃ | U_x U_y U_z | p q r | λ₀ λ₁ λ₂ λ₃ | φ | s
DEG = pi / 180;
F = P.cfg.fire;
ax = zrAxesF(0, P.QE);
a0 = F.alphaTot0_deg * DEG; f0 = F.phiA0_deg * DEG;
xF = cos(a0) * ax.x + sin(a0) * (cos(f0) * (-ax.z) + sin(f0) * ax.y);
theta = asin(max(-1, min(1, xF(2)))); psi = atan2(xF(3), xF(1));
lam = quatFromEuler(psi, theta, 0);
T = quatToT(lam);
UI = [ax.x(1) * P.u0, ax.x(3) * P.u0, -ax.x(2) * P.u0];
y = zeros(15, 1);
y(4) = T(1) * UI(1) + T(2) * UI(2) + T(3) * UI(3);
y(5) = T(4) * UI(1) + T(5) * UI(2) + T(6) * UI(3);
y(6) = T(7) * UI(1) + T(8) * UI(2) + T(9) * UI(3);
y(7) = P.p0; y(8) = F.q0; y(9) = F.r0;
y(10:13) = lam(:);
end
