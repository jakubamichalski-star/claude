function [dy, clamp] = rhs(P, wind, t, y)
%RHS Prawa strona równań ruchu (układ zero-roll).
%   [dy, clamp] = rhs(P, wind, t, y)
%   Stan y: X₁ X₂ X₃ | U_x U_y U_z | p q r | λ₀ λ₁ λ₂ λ₃ | φ | s
%   dX̄/dt = Tᵀ·Ū (7); dŪ/dt = F̄/m + ḡ + Ā_C − Ω̄ × Ū (4), (9);
%   I·dϖ̄/dt = M̄ − Ω̄ × (I·ϖ̄) (5), (10)–(11); dλ̄/dt wg (8);
%   Ω̄ = (−r·tan θ, q, r) (A.10); dφ/dt = p + r·tan θ.
l0 = y(10); l1 = y(11); l2 = y(12); l3 = y(13);
kq = 1 / (l0 * l0 + l1 * l1 + l2 * l2 + l3 * l3);
T00 = (l0 * l0 + l1 * l1 - l2 * l2 - l3 * l3) * kq; T01 = 2 * (l1 * l2 + l0 * l3) * kq; T02 = 2 * (l1 * l3 - l0 * l2) * kq;
T10 = 2 * (l1 * l2 - l0 * l3) * kq; T11 = (l0 * l0 - l1 * l1 + l2 * l2 - l3 * l3) * kq; T12 = 2 * (l2 * l3 + l0 * l1) * kq;
T20 = 2 * (l1 * l3 + l0 * l2) * kq; T21 = 2 * (l2 * l3 - l0 * l1) * kq; T22 = (l0 * l0 - l1 * l1 - l2 * l2 + l3 * l3) * kq;
tth = -T02 / sqrt(T00 * T00 + T01 * T01);        % tan θ
X2 = y(2);
Ux = y(4); Uy = y(5); Uz = y(6); p = y(7); q = y(8); r = y(9);
if P.frozen, hA = P.hM; else, hA = P.hM + X2; end
[~, ~, rho, aS] = atmAt(P, hA);
W = wind.W(t, X2);                                 % W̄ w układzie ogniowym
WI0 = W(1); WI1 = W(3); WI2 = -W(2);               % → układ I
Vx = Ux - (T00 * WI0 + T01 * WI1 + T02 * WI2);     % V̄ = Ū − W̄ (B.1), składowe zero-roll
Vy = Uy - (T10 * WI0 + T11 * WI1 + T12 * WI2);
Vz = Uz - (T20 * WI0 + T21 * WI1 + T22 * WI2);
V = sqrt(Vx * Vx + Vy * Vy + Vz * Vz);
alphaTot = atan2(sqrt(Vy * Vy + Vz * Vz), Vx);     % (B.5) = acos(V_x/V), postać atan2 dokładna przy małych kątach
c = aeroCoefVec(P, V / aS, alphaTot);   % [CD CLa CMa CMq CNq Cmagf Cmagm Cspin CNa CPN clamp]
qS = 0.5 * rho * P.S; d = P.d;
% siły: DF̄ (B.10), LF̄ (B.11), MF̄ (B.15), PDF̄ (B.13); x̄ = (1, 0, 0)
kD = -qS * c(1) * V; kL = qS * c(2); kMF = qS * p * d * c(6); kPF = qS * d * c(5) * V;
Fx = kD * Vx + kL * (Vy * Vy + Vz * Vz);
Fy = kD * Vy - kL * Vx * Vy + kMF * Vz + kPF * r;
Fz = kD * Vz - kL * Vx * Vz - kMF * Vy - kPF * q;
% momenty: OM̄ (B.12), PDM̄ (B.14), MM̄ (B.16), SDM̄ (B.18); składowe L, M, N
kOM = qS * d * c(3) * V; kPM = qS * d * d * c(4) * V; kMM = qS * d * d * p * c(7);
Lm = qS * d * d * p * V * c(8);
Mm = kOM * Vz + kPM * q + kMM * Vy;
Nm = -kOM * Vy + kPM * r + kMM * Vz;
% ḡ: w układzie I (0, 0, g)
g = P.g; gx = T02 * g; gy = T12 * g; gz = T22 * g;
Ax = 0; Ay = 0; Az = 0;                            % Ā_C = −2·Λ̄ × Ū (opcja)
if P.coriolis
  L = P.LamI;
  Lx = T00 * L(1) + T01 * L(2) + T02 * L(3); Ly = T10 * L(1) + T11 * L(2) + T12 * L(3); Lz = T20 * L(1) + T21 * L(2) + T22 * L(3);
  Ax = -2 * (Ly * Uz - Lz * Uy); Ay = -2 * (Lz * Ux - Lx * Uz); Az = -2 * (Lx * Uy - Ly * Ux);
end
Ox = -r * tth;                                     % Ω_x układu zero-roll
m = P.m; Ix = P.Ix; IT = P.IT;
dy = zeros(15, 1);
dy(4) = Fx / m + gx + Ax - (q * Uz - r * Uy);
dy(5) = Fy / m + gy + Ay - (r * Ux - Ox * Uz);
dy(6) = Fz / m + gz + Az - (Ox * Uy - q * Ux);
dy(7) = Lm / Ix;
dy(8) = (Mm - Ix * p * r + IT * Ox * r) / IT;      % I_T·q̇ = M − I_x·p·r − I_T·r²·tan θ
dy(9) = (Nm + Ix * p * q - IT * Ox * q) / IT;      % I_T·ṙ = N + I_x·p·q + I_T·q·r·tan θ
dy(10) = 0.5 * (-l1 * Ox - l2 * q - l3 * r);       % (8)
dy(11) = 0.5 * (l0 * Ox - l3 * q + l2 * r);
dy(12) = 0.5 * (l3 * Ox + l0 * q - l1 * r);
dy(13) = 0.5 * (-l2 * Ox + l1 * q + l0 * r);
XI0 = T00 * Ux + T10 * Uy + T20 * Uz; XI1 = T01 * Ux + T11 * Uy + T21 * Uz; XI2 = T02 * Ux + T12 * Uy + T22 * Uz;
dy(1) = XI0; dy(2) = -XI2; dy(3) = XI1;            % I → układ ogniowy
dy(14) = p - Ox;                                   % dφ/dt = p + r·tan θ
dy(15) = sqrt(Ux * Ux + Uy * Uy + Uz * Uz);
clamp = c(11);
end
