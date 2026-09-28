function L = pokazWylot(A, cfg)
%POKAZWYLOT Wielkości u wylotu (teoria liniowa, α_tot → 0) – liczone bez całkowania toru.
%   L = pokazWylot(A, cfg) – wypisuje tabelę i zwraca wynik linPoint.
P = makeParams(A, cfg);
[~, ~, rho, aS] = atmAt(P, P.hM);
M0 = P.u0 / aS;
L = linPoint(P, M0, P.u0, rho, P.p0, 0);
B = stabBands(P, M0, 0, L);
md2 = P.m * P.d^2;
cond = L.Sg > 1 && 1 / L.Sg < L.Sd0 * (2 - L.Sd0);
wmax = max([abs(L.sF0), abs(L.sS0), abs(L.wg), 1e-3]);
dt = min(max(2 * pi / (P.cfg.num.stepsPerPeriod * wmax), 1e-6), 5e-3);
if isempty(P.met), zm = ''; else, zm = ' (z warstw meteo)'; end
if cond, cs = 'spełniony'; else, cs = 'NIESPEŁNIONY'; end
fprintf('\nWielkości u wylotu (teoria liniowa, α_tot → 0)\n');
rows = {
  'p0', sprintf('%s rad/s (%s obr/s)', fmtPL(P.p0, 1), fmtPL(P.p0 / (2 * pi), 1))
  'ℳ0, ρ0', sprintf('%s, %s kg/m³%s', fmtPL(M0, 3), fmtPL(rho, 4), zm)
  'md²/Ix, md²/IT', sprintf('%s, %s', fmtPL(md2 / P.Ix, 3), fmtPL(md2 / P.IT, 4))
  'Ix/IT', fmtPL(P.Ix / P.IT, 5)
  'ΔX_CG', sprintf('%s kal. (%s mm)', fmtPL(P.dXcg, 4), fmtPL(P.dXcg * P.d * 1000, 2))
  'S_g', sprintf('%s ± %s', fmtPL(L.Sg, 3), fmtPL(B.SgSig, 3))
  'S_d (α_tot → 0)', sprintf('%s ± %s', fmtPL(L.Sd0, 3), fmtPL(B.SdSig, 3))
  'warunek 1/S_g < S_d(2 − S_d)', sprintf('%s < %s: %s', fmtPL(1 / L.Sg, 4), fmtPL(L.Sd0 * (2 - L.Sd0), 4), cs)
  'ν_F, ν_S', sprintf('%s Hz, %s Hz', fmtPL(L.nuF0, 2), fmtPL(L.nuS0, 2))
  'okres nutacji 1/ν_F', sprintf('%s ms', fmtPL(1000 / L.nuF0, 1))
  'σ_F, σ_S', sprintf('%s 1/s, %s 1/s', fmtPL(L.sigF0, 3), fmtPL(L.sigS0, 3))
  'krok całkowania dt', sprintf('%s ms', fmtPL(dt * 1000, 4))};
for k = 1:size(rows, 1), fprintf('  %-30s %s\n', rows{k, 1}, rows{k, 2}); end
end
