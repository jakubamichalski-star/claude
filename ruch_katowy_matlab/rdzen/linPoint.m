function L = linPoint(P, M, V, rho, p, alphaTot)
%LINPOINT Teoria liniowa (quasi-stacjonarna) w chwili lotu – sekcja 4.
%   ξ̈ + (η − i·ω_g)·ξ̇ − (ω_α² + i·ω_g·ζ)·ξ = wymuszenie
%   ω_g = I_x·p/I_T,  ω_α² = ρ·S·d·V²·C_Mα/(2·I_T),
%   η = ρSV/(2m)·[C_Lα − C_D0 − (md²/I_T)·(C_Mq + C_Mα̇)],
%   ζ = ρSV/(2m)·[C_Lα + (md²/I_x)·C_mag-m],
%   S_g = ω_g²/(4ω_α²) i S_d = 2ζ/η – tożsame ze wzorami sekcji 4.
%   Pierwiastki dokładne s_F, s_S (F: większe |Im s|); wzory zamknięte ν_F,S, σ_F,S.
%   Sufiks 0 – C_mag-m przy α_tot → 0, sufiks A – przy α_tot (sieczna).
c0 = aeroCoef(P, M, 0);
CLa = c0.CLa; CD0 = c0.CD; CMa = c0.CMa; CMq = c0.CMq; Cmagm0 = c0.Cmagm;
if alphaTot > 0, cA = aeroCoef(P, M, alphaTot); CmagmA = cA.Cmagm; else, CmagmA = Cmagm0; end
m = P.m; Ix = P.Ix; IT = P.IT; d = P.d; S = P.S; md2 = m * d * d;
wg = Ix * p / IT;
wa2 = rho * S * d * V * V * CMa / (2 * IT);
base = rho * S * V / (2 * m);
Hs = CLa - CD0 - (md2 / IT) * CMq;
eta = base * Hs; etaQ = -base * (md2 / IT) * CMq;
zeta0 = base * (CLa + (md2 / Ix) * Cmagm0); zetaA = base * (CLa + (md2 / Ix) * CmagmA);
if wg ~= 0, invSg = 4 * wa2 / (wg * wg);
elseif wa2 > 0, invSg = Inf;
elseif wa2 < 0, invSg = -Inf;
else, invSg = 0;
end
Sd0 = 2 * zeta0 / eta; SdA = 2 * zetaA / eta;
k = Ix * p / (4 * pi * IT);
if invSg < 1, rt = sqrt(1 - invSg); else, rt = NaN; end
B = complex(eta, -wg);
r0 = orderRoots(quadRoots(B, complex(-wa2, -wg * zeta0)));
rA = orderRoots(quadRoots(B, complex(-wa2, -wg * zetaA)));
L = struct('CLa', CLa, 'CD0', CD0, 'CMa', CMa, 'CMq', CMq, 'Cmagm0', Cmagm0, 'CmagmA', CmagmA, 'Cmagf', c0.Cmagf, ...
  'Cspin', c0.Cspin, 'CNq', c0.CNq, 'CNa', c0.CNa, 'CPN', c0.CPN, 'clamp', c0.clamp, ...
  'wg', wg, 'wa2', wa2, 'base', base, 'Hs', Hs, 'eta', eta, 'etaQ', etaQ, 'zeta0', zeta0, 'zetaA', zetaA, ...
  'invSg', invSg, 'Sg', 1 / invSg, 'Sd0', Sd0, 'SdA', SdA, ...
  'nuF0', k * (1 + rt), 'nuS0', k * (1 - rt), ...
  'sigF0', -(eta / 2) * (1 - (Sd0 - 1) / rt), 'sigS0', -(eta / 2) * (1 + (Sd0 - 1) / rt), ...
  'nuFA', k * (1 + rt), 'nuSA', k * (1 - rt), ...
  'sigFA', -(eta / 2) * (1 - (SdA - 1) / rt), 'sigSA', -(eta / 2) * (1 + (SdA - 1) / rt), ...
  'sF0', r0(1), 'sS0', r0(2), 'sFA', rA(1), 'sSA', rA(2));
end
