function fig = plotTimeSeries(R, tRange)
%PLOTTIMESERIES Przebiegi czasowe: pasek stanu, α_tot, amplitudy ramion, składowe ξ, S_g, S_d, σ, ν, ℳ, wiatr.
%   fig = plotTimeSeries(R, tRange) – tRange = [t0 t1] (domyślnie cały lot); osie czasu połączone.
DEG = pi / 180;
ch = R.ch; po = R.post; t = ch.t;
if nargin < 2 || isempty(tRange), tRange = [t(1), t(end)]; end
if tRange(2) <= tRange(1), tRange(2) = tRange(1) + 1e-3; end
cR = [0.72 0.18 0.18]; cB = [0.12 0.37 0.55]; cO = [0.85 0.50 0.10]; cG = [0.20 0.55 0.33]; cP = [0.45 0.25 0.60];
cI = [0.11 0.15 0.22]; cT = [0.10 0.55 0.55]; cGr = [0.45 0.47 0.50];
fig = figure('Name', 'Przebiegi czasowe', 'Color', 'w', 'Units', 'normalized', 'Position', [0.05 0.05 0.85 0.85]);
% --- pasek stanu ---
ax(1) = subplot(6, 2, [1 2]); hold(ax(1), 'on');
pasek(ax(1), R.summary.runsA, 0.55, 1.0); pasek(ax(1), R.summary.runs0, 0.0, 0.45);
set(ax(1), 'YTick', [0.225 0.775], 'YTickLabel', {'α → 0', 'ampl.'}, 'YLim', [0 1], 'Box', 'on');
title(ax(1), 'Stan wzdłuż toru: zielony – stateczny, żółty – słaba nadążność, pomarańczowy – cykl graniczny, czerwony – niestateczny dynamicznie, bordowy – giroskopowo');
for e = R.meta.events
  if strcmp(e.type, 'apex'), c = cG; else, c = cT; end
  plot(ax(1), [e.t e.t], [0 1], '-', 'Color', c, 'LineWidth', 1.2);
end
% --- α_tot ---
ax(2) = subplot(6, 2, 3); hold(ax(2), 'on');
pasmo(ax(2), t, po.envM / DEG, po.envP / DEG, cB, 0.13);
plot(ax(2), t, ch.alphaTot / DEG, 'Color', cI, 'LineWidth', 1.1, 'DisplayName', 'α_{tot} (6DOF)');
plot(ax(2), t, po.devE / DEG, 'Color', cB, 'LineWidth', 1.0, 'DisplayName', '|ξ − ξ_e|');
plot(ax(2), t, po.aeAbs / DEG, 'Color', cG, 'LineWidth', 1.5, 'DisplayName', '|ᾱ_e| STANAG');
plot(ax(2), t, abs(complex(po.xieRe, po.xieIm)) / DEG, '--', 'Color', cP, 'LineWidth', 1.1, 'DisplayName', '|ξ_e|');
title(ax(2), 'Kąt nutacji α_{tot} [°] (pasmo: obwiednie ||K_F| ∓ |K_S||)'); legend(ax(2), 'Location', 'best');
% --- amplitudy ramion ---
ax(3) = subplot(6, 2, 4); hold(ax(3), 'on');
fl = @(v) max(v / DEG, 1e-3);
semilogy(ax(3), t, fl(po.KFm), 'Color', cR, 'LineWidth', 1.4, 'DisplayName', '|K_F| z 6DOF');
semilogy(ax(3), t, fl(po.KSm), 'Color', cB, 'LineWidth', 1.4, 'DisplayName', '|K_S| z 6DOF');
semilogy(ax(3), t, fl(po.KFabs), '--', 'Color', cR, 'LineWidth', 1.1, 'DisplayName', '|K_F| teoria liniowa');
semilogy(ax(3), t, fl(po.KSabs), '--', 'Color', cB, 'LineWidth', 1.1, 'DisplayName', '|K_S| teoria liniowa');
set(ax(3), 'YScale', 'log'); title(ax(3), 'Amplitudy ramion |K_F|, |K_S| [°]'); legend(ax(3), 'Location', 'best');
% --- składowe ξ ---
ax(4) = subplot(6, 2, 5); hold(ax(4), 'on');
plot(ax(4), t, po.xiRe / DEG, 'Color', cI, 'DisplayName', 'Re ξ');
plot(ax(4), t, po.xiIm / DEG, 'Color', cB, 'DisplayName', 'Im ξ');
title(ax(4), 'Składowe ξ [°]: Re – nos w górę, Im – nos w prawo'); legend(ax(4), 'Location', 'best');
% --- S_g ---
ax(5) = subplot(6, 2, 6); hold(ax(5), 'on');
pasmo(ax(5), t, max(po.Sg - po.SgSig, 1e-3), po.Sg + po.SgSig, cI, 0.12);
semilogy(ax(5), t, po.Sg, 'Color', cI, 'LineWidth', 1.5, 'DisplayName', 'S_g');
set(ax(5), 'YScale', 'log'); title(ax(5), 'Stateczność giroskopowa S_g (pasmo ±1σ z pliku)');
% --- S_d ---
ax(6) = subplot(6, 2, 7); hold(ax(6), 'on');
pasmo(ax(6), t, po.Sd0 - po.SdSig, po.Sd0 + po.SdSig, cB, 0.13);
plot(ax(6), t, po.Sd0, 'Color', cB, 'LineWidth', 1.5, 'DisplayName', 'S_d (α_{tot} → 0)');
plot(ax(6), t, po.SdFA, 'Color', cR, 'LineWidth', 1.3, 'DisplayName', 'S_{d,F} (ampl.)');
plot(ax(6), t, po.SdSA, 'Color', cO, 'LineWidth', 1.3, 'DisplayName', 'S_{d,S} (ampl.)');
title(ax(6), 'Stateczność dynamiczna S_d'); legend(ax(6), 'Location', 'best');
yl = [min([po.Sd0; po.SdFA; po.SdSA; 0]), max([po.Sd0; po.SdFA; po.SdSA; 2])];
ylim(ax(6), [max(yl(1), -4) - 0.1, min(yl(2), 6) + 0.1]);
% --- σ ---
ax(7) = subplot(6, 2, 8); hold(ax(7), 'on');
plot(ax(7), t, po.sigF0, 'Color', cR, 'LineWidth', 1.4, 'DisplayName', 'σ_F (α → 0)');
plot(ax(7), t, po.sigS0, 'Color', cB, 'LineWidth', 1.4, 'DisplayName', 'σ_S (α → 0)');
plot(ax(7), t, po.sFreA, '--', 'Color', cR, 'LineWidth', 1.2, 'DisplayName', 'σ_F (ampl.)');
plot(ax(7), t, po.sSreA, '--', 'Color', cB, 'LineWidth', 1.2, 'DisplayName', 'σ_S (ampl.)');
title(ax(7), 'Tłumienie σ [1/s]'); legend(ax(7), 'Location', 'best');
% --- ν ---
ax(8) = subplot(6, 2, 9); hold(ax(8), 'on');
semilogy(ax(8), t, po.nuF0, 'Color', cR, 'LineWidth', 1.4, 'DisplayName', 'ν_F');
semilogy(ax(8), t, po.nuS0, 'Color', cB, 'LineWidth', 1.4, 'DisplayName', 'ν_S');
set(ax(8), 'YScale', 'log'); title(ax(8), 'Częstotliwości ν [Hz]'); legend(ax(8), 'Location', 'best');
% --- ℳ ---
ax(9) = subplot(6, 2, 10); hold(ax(9), 'on');
plot(ax(9), t, ch.M, 'Color', cI, 'LineWidth', 1.4);
title(ax(9), 'Liczba Macha ℳ');
% --- wiatr ---
ax(10) = subplot(6, 2, 11); hold(ax(10), 'on');
plot(ax(10), t, ch.W1, 'Color', cGr, 'LineWidth', 1.3, 'DisplayName', 'W_1');
plot(ax(10), t, ch.W2, 'Color', cG, 'LineWidth', 1.3, 'DisplayName', 'W_2');
plot(ax(10), t, ch.W3, 'Color', cT, 'LineWidth', 1.5, 'DisplayName', 'W_3');
title(ax(10), 'Wiatr W̄ [m/s] (układ ogniowy)'); legend(ax(10), 'Location', 'best'); xlabel(ax(10), 't [s]');
% --- tor h(t) ---
ax(11) = subplot(6, 2, 12); hold(ax(11), 'on');
plot(ax(11), t, ch.X2, 'Color', cI, 'LineWidth', 1.4);
title(ax(11), 'Wysokość X_2 [m]'); xlabel(ax(11), 't [s]');
for k = 1:numel(ax), grid(ax(k), 'on'); set(ax(k), 'XLim', tRange); end
linkaxes(ax, 'x'); xlim(ax(1), tRange);
liniaPozioma(ax(5), 1); liniaPozioma(ax(6), 0); liniaPozioma(ax(6), 2); liniaPozioma(ax(7), 0); liniaPozioma(ax(9), 1, [0.56 0.6 0.65]);
end

function pasek(ax, rs, y0, y1)
for r = rs
  t1 = r.t1; if t1 <= r.t0, t1 = r.t0 + 1e-6; end
  patch(ax, [r.t0 t1 t1 r.t0], [y0 y0 y1 y1], kolorStanu(r.key), 'EdgeColor', 'none');
end
end
