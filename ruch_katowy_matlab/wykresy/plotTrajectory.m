function fig = plotTrajectory(R, tc)
%PLOTTRAJECTORY Tor: rzut 1–2 (wysokość) i rzut 1–3 (derywacja); wierzchołek, podmuchy i chwila tc.
ch = R.ch; t = ch.t;
if nargin < 2 || isempty(tc), tc = t(end); end
i = find(t <= tc, 1, 'last');
fig = figure('Name', 'Tor', 'Color', 'w');
for k = 1:2
  ax = subplot(2, 1, k); hold(ax, 'on'); grid(ax, 'on');
  if k == 1, y = ch.X2 / 1000; yl = 'X_2 [km]'; tt = 'Tor: rzut 1–2'; else, y = ch.X3; yl = 'X_3 [m]'; tt = 'Tor: rzut 1–3 (derywacja)'; end
  plot(ax, ch.X1 / 1000, y, 'Color', [0.11 0.15 0.22], 'LineWidth', 1.5);
  for e = R.meta.events
    j = find(t <= e.t, 1, 'last');
    if strcmp(e.type, 'apex')
      plot(ax, ch.X1(j) / 1000, y(j), 'o', 'MarkerFaceColor', 'w', 'MarkerEdgeColor', [0.2 0.55 0.33]);
      text(ax, ch.X1(j) / 1000, y(j), 'wierzchołek', 'Color', [0.2 0.55 0.33], 'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'center');
    else
      plot(ax, ch.X1(j) / 1000, y(j), '.', 'Color', [0.1 0.55 0.55], 'MarkerSize', 14);
    end
  end
  plot(ax, ch.X1(i) / 1000, y(i), 'o', 'MarkerFaceColor', [0.72 0.18 0.18], 'MarkerEdgeColor', 'w', 'MarkerSize', 8);
  xlabel(ax, 'X_1 [km]'); ylabel(ax, yl); title(ax, tt);
end
end
