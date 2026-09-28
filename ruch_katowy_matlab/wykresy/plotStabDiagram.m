function fig = plotStabDiagram(R, tc)
%PLOTSTABDIAGRAM Diagram stateczności S_d – 1/S_g wzdłuż toru.
%   Obszar stateczności: 0 < 1/S_g < S_d(2 − S_d); lewa gałąź paraboli S_d = 1 − √(1 − 1/S_g) –
%   granica składowej wolnej (S), prawa S_d = 1 + √(1 − 1/S_g) – szybkiej (F).
%   Krzywe: S_d (α_tot → 0), S_d,F i S_d,S (przy bieżącej amplitudzie); kółka – wylot i chwila tc.
po = R.post; t = R.ch.t;
if nargin < 2 || isempty(tc), tc = t(end); end
i = find(t <= tc, 1, 'last');
fig = figure('Name', 'Diagram S_d – 1/S_g', 'Color', 'w'); ax = axes('Parent', fig); hold(ax, 'on'); grid(ax, 'on');
v = [po.Sd0; po.SdFA; po.SdSA]; v = v(isfinite(v));
x0 = max(min([-0.5; v]) - 0.1, -4); x1 = min(max([2.5; v]) + 0.1, 6);
y1 = min(max([1.15; 1.05 * po.invSg(isfinite(po.invSg))]), 3);
s = linspace(0, 2, 201);
fill(ax, [s, 0], [s .* (2 - s), 0], [0.24 0.55 0.35], 'FaceAlpha', 0.10, 'EdgeColor', 'none', 'HandleVisibility', 'off');
fill(ax, [x0 x1 x1 x0], [1 1 y1 y1], [0.43 0.1 0.16], 'FaceAlpha', 0.07, 'EdgeColor', 'none', 'HandleVisibility', 'off');
plot(ax, s, s .* (2 - s), 'Color', [0.24 0.55 0.35], 'LineWidth', 1.5, 'DisplayName', 'granica stateczności');
plot(ax, [x0 x1], [1 1], '--', 'Color', [0.55 0.23 0.29], 'HandleVisibility', 'off');
text(ax, 1, 0.42, 'stateczny', 'Color', [0.18 0.43 0.28], 'HorizontalAlignment', 'center');
text(ax, 0.22, 0.45, 'granica S', 'Color', [0.36 0.39 0.45], 'Rotation', 60, 'HorizontalAlignment', 'center');
text(ax, 1.78, 0.45, 'granica F', 'Color', [0.36 0.39 0.45], 'Rotation', -60, 'HorizontalAlignment', 'center');
text(ax, x0 + 0.05, 1.03, 'S_g < 1', 'Color', [0.55 0.23 0.29], 'VerticalAlignment', 'bottom');
spec = {'SdSA', [0.85 0.5 0.1], 'S_{d,S} (ampl.)'; 'SdFA', [0.72 0.18 0.18], 'S_{d,F} (ampl.)'; 'Sd0', [0.12 0.37 0.55], 'S_d (α_{tot} → 0)'};
for k = 1:3
  x = po.(spec{k, 1});
  plot(ax, x, po.invSg, 'Color', spec{k, 2}, 'LineWidth', 1.4, 'DisplayName', spec{k, 3});
  plot(ax, x(1), po.invSg(1), 'o', 'MarkerFaceColor', 'w', 'MarkerEdgeColor', spec{k, 2}, 'HandleVisibility', 'off');
  plot(ax, x(i), po.invSg(i), 'o', 'MarkerFaceColor', spec{k, 2}, 'MarkerEdgeColor', 'w', 'MarkerSize', 8, 'HandleVisibility', 'off');
end
xlim(ax, [x0 x1]); ylim(ax, [-0.03 y1]);
xlabel(ax, 'S_d'); ylabel(ax, '1/S_g');
title(ax, sprintf('Diagram S_d – 1/S_g (pełne kółka: t = %.2f s, puste: wylot)', t(i)));
legend(ax, 'Location', 'northeast');
end
