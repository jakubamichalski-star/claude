function fig = plotMeteo(cfg)
%PLOTMETEO Profile warstw METCM po interpolacji (PCHIP): W1, W3 (układ ogniowy), ΔT, Δρ w funkcji wysokości nad MDP.
%   fig = plotMeteo(cfg) – cfg z polami met i earth.az_deg (jak w defaultConfig).
met = cfg.met; met.enabled = true;
hDef = 0; if isfield(cfg, 'fire'), hDef = cfg.fire.hM; end
M = metProfile(met, cfg.earth.az_deg, hDef);
h = linspace(0, max(M.xs) + 1000, 600);
W = zeros(numel(h), 3); dT = zeros(size(h)); dR = dT;
for k = 1:numel(h), W(k, :) = metW(M, h(k)).'; dT(k) = pchipVal(M.fT, h(k)); dR(k) = pchipVal(M.fR, h(k)); end
fig = figure('Name', 'Profile meteo', 'Color', 'w', 'Units', 'normalized', 'Position', [0.1 0.1 0.7 0.6]);
ax = subplot(1, 3, 1); hold(ax, 'on'); grid(ax, 'on');
plot(ax, W(:, 1), h / 1000, 'Color', [0.45 0.47 0.5], 'LineWidth', 1.4, 'DisplayName', 'W_1 (wzdłuż strzału)');
plot(ax, W(:, 3), h / 1000, 'Color', [0.1 0.55 0.55], 'LineWidth', 1.4, 'DisplayName', 'W_3 (w prawo)');
plot(ax, M.w1, M.xs / 1000, 'o', 'Color', [0.45 0.47 0.5], 'HandleVisibility', 'off');
plot(ax, M.w3, M.xs / 1000, 'o', 'Color', [0.1 0.55 0.55], 'HandleVisibility', 'off');
xlabel(ax, 'wiatr [m/s]'); ylabel(ax, 'wysokość nad MDP [km]'); legend(ax, 'Location', 'best'); title(ax, 'Wiatr (układ ogniowy)');
ax = subplot(1, 3, 2); hold(ax, 'on'); grid(ax, 'on');
plot(ax, dT, h / 1000, 'Color', [0.72 0.18 0.18], 'LineWidth', 1.4); plot(ax, M.dT, M.xs / 1000, 'o', 'Color', [0.72 0.18 0.18]);
xlabel(ax, 'ΔT [K]'); title(ax, 'ΔT względem ISO 2533');
ax = subplot(1, 3, 3); hold(ax, 'on'); grid(ax, 'on');
plot(ax, dR, h / 1000, 'Color', [0.12 0.37 0.55], 'LineWidth', 1.4); plot(ax, M.dR, M.xs / 1000, 'o', 'Color', [0.12 0.37 0.55]);
xlabel(ax, 'Δρ [%]'); title(ax, 'Δρ względem ISO 2533');
end
