function plotResults(R, opts)
%PLOTRESULTS Komplet wykresów wyniku obliczenia toru.
%   plotResults(R, opts) – opts (wszystkie pola opcjonalne):
%     tRange  – zakres czasu wykresów czasowych [t0 t1] (domyślnie cały lot)
%     cursor  – chwila dla rozety, diagramu stateczności i toru (domyślnie koniec zakresu)
%     rosWin  – okno rozety [s] (domyślnie 1 s; 0 – od wylotu)
if nargin < 2, opts = struct(); end
t = R.ch.t;
tRange = getf(opts, 'tRange', [t(1), t(end)]);
tc = getf(opts, 'cursor', tRange(2));
plotTimeSeries(R, tRange);
plotRosette(R, tc, getf(opts, 'rosWin', 1));
plotStabDiagram(R, tc);
plotTrajectory(R, tc);
end

function v = getf(s, f, d)
if isfield(s, f) && ~isempty(s.(f)), v = s.(f); else, v = d; end
end
