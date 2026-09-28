function fig = plotRosette(R, tc, win, pokazLin)
%PLOTROSETTE Rozeta ξ (widok od tyłu): ślad końca osi względem V̄ w oknie [tc − win, tc].
%   fig = plotRosette(R, tc, win, pokazLin)
%   tc – chwila (domyślnie koniec lotu), win – okno [s] (0 lub Inf – od wylotu, domyślnie 1 s),
%   pokazLin – także rozwiązanie liniowe ξ_lin (domyślnie true).
%   Oznaczenia: × – ξ_e, strzałki – ramię wolne K_S (niebieska) i szybkie K_F (czerwona),
%   kwadrat – ᾱ_e (STANAG), kółko – ξ(tc) z 6DOF.
DEG = pi / 180; ch = R.ch; po = R.post; t = ch.t;
if nargin < 2 || isempty(tc), tc = t(end); end
if nargin < 3 || isempty(win), win = 1; end
if nargin < 4, pokazLin = true; end
if win > 0 && isfinite(win), ta = max(t(1), tc - win); else, ta = t(1); end
i0 = find(t <= ta, 1, 'last'); i1 = find(t <= tc, 1, 'last');
fig = figure('Name', 'Rozeta ξ', 'Color', 'w'); ax = axes('Parent', fig); hold(ax, 'on'); axis(ax, 'equal'); grid(ax, 'on');
re = po.xiRe / DEG; im = po.xiIm / DEG;
if pokazLin
  plot(ax, po.xiLinIm(i0:i1) / DEG, po.xiLinRe(i0:i1) / DEG, '--', 'Color', [0.85 0.5 0.1], 'LineWidth', 1.0, 'DisplayName', 'ξ_{lin} (teoria liniowa)');
end
% 6DOF z zanikającą historią (starsze fragmenty jaśniejsze)
nseg = 16; span = max(1, i1 - i0);
for s = 0:nseg - 1
  a = i0 + floor(span * s / nseg); b = min(i1, i0 + floor(span * (s + 1) / nseg));
  al = 0.16 + 0.84 * (s + 1) / nseg;
  c = 1 - al * (1 - [0.11 0.15 0.22]);
  if s == nseg - 1, hv = 'on'; else, hv = 'off'; end
  plot(ax, im(a:b), re(a:b), 'Color', c, 'LineWidth', 1.3, 'HandleVisibility', hv, 'DisplayName', 'ξ z 6DOF');
end
xe = [po.xieIm(i1), po.xieRe(i1)] / DEG; KS = [po.KSim(i1), po.KSre(i1)] / DEG; KF = [po.KFim(i1), po.KFre(i1)] / DEG;
if all(isfinite(xe))
  plot(ax, xe(1), xe(2), 'x', 'Color', [0.45 0.25 0.6], 'MarkerSize', 10, 'LineWidth', 2, 'DisplayName', 'ξ_e');
  if all(isfinite(KS))
    quiver(ax, xe(1), xe(2), KS(1), KS(2), 0, 'Color', [0.12 0.37 0.55], 'LineWidth', 2, 'MaxHeadSize', 0.5, 'DisplayName', 'K_S (precesja)');
    quiver(ax, xe(1) + KS(1), xe(2) + KS(2), KF(1), KF(2), 0, 'Color', [0.72 0.18 0.18], 'LineWidth', 2, 'MaxHeadSize', 0.5, 'DisplayName', 'K_F (nutacja)');
  end
end
plot(ax, po.aeIm(i1) / DEG, po.aeRe(i1) / DEG, 's', 'Color', [0.2 0.55 0.33], 'MarkerSize', 9, 'LineWidth', 1.6, 'DisplayName', 'ᾱ_e STANAG');
plot(ax, im(i1), re(i1), 'o', 'MarkerFaceColor', 'w', 'MarkerEdgeColor', [0.11 0.15 0.22], 'MarkerSize', 7, 'LineWidth', 1.5, 'DisplayName', 'ξ(t)');
% okręgi stałego α_tot wokół początku układu
xl = xlim(ax); yl = ylim(ax);
rmax = hypot(max(abs(xl)), max(abs(yl))); st = niceStepR(rmax / 5); ph = linspace(0, 2 * pi, 181);
for r = st:st:rmax + st
  plot(ax, r * cos(ph), r * sin(ph), '-', 'Color', [0.89 0.91 0.93], 'HandleVisibility', 'off');
end
xlim(ax, xl); ylim(ax, yl);
xlabel(ax, 'Im ξ – nos w prawo [°]'); ylabel(ax, 'Re ξ – nos w górę [°]');
if win > 0 && isfinite(win), okno = sprintf('okno %.2f s', min(win, tc - t(1))); else, okno = 'od wylotu'; end
title(ax, sprintf('Rozeta ξ (widok od tyłu), t = %.3f s, %s', tc, okno));
legend(ax, 'Location', 'bestoutside');
end

function s = niceStepR(x)
if ~(x > 0), s = 1; return; end
e = 10^floor(log10(x)); m = x / e;
if m <= 1, s = e; elseif m <= 2, s = 2 * e; elseif m <= 5, s = 5 * e; else, s = 10 * e; end
end
