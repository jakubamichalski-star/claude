function R = uruchomScenariusz(A, id, wariant, cfgBiezacy, rysuj)
%URUCHOMSCENARIUSZ Uruchomienie scenariusza (a)–(f).
%   R = uruchomScenariusz(A, id, wariant, cfgBiezacy, rysuj)
%   id – 'a', 'b', 'c1', 'c2', 'd', 'e', 'f1', 'f2'; wariant – [] (podstawowy) albo numer wariantu (1, …).
%   Ustawienia: wartości domyślne plus zmiany scenariusza; masa, bezwładność, środek masy i normowanie
%   z cfgBiezacy (dane pocisku użytkownika), warstwy meteo wyłączone.
%   rysuj (domyślnie true) – wykresy z zakresem czasu właściwym dla scenariusza.
if nargin < 3, wariant = []; end
if nargin < 4 || isempty(cfgBiezacy), cfgBiezacy = defaultConfig(); end
if nargin < 5, rysuj = true; end
SCN = scenariusze();
s = SCN(strcmp({SCN.id}, id));
if isempty(s), error('Nieznany scenariusz „%s”. Dostępne: %s.', id, strjoin({SCN.id}, ', ')); end
over = s.over; widok = s.view; tytul = s.title;
if ~isempty(wariant)
  v = s.variants{wariant}; over = v.over; widok = v.view; tytul = [tytul ' – ' v.label];
end
cfg = defaultConfig();
cfg.mass = cfgBiezacy.mass; cfg.norm = cfgBiezacy.norm;
if isfield(cfgBiezacy, 'met'), cfg.met = cfgBiezacy.met; cfg.met.enabled = false; end   % tabela meteo zostaje, warstwy wyłączone
cfg = deepMerge(cfg, over);
fprintf('\n==== Scenariusz %s ====\n%s\nUstawienia: %s\n', tytul, s.goal, s.set);
R = runModel(A, cfg);
rows = s.res(R);
fprintf('\nWynik:\n');
for k = 1:size(rows, 1), fprintf('  %-60s %s\n', rows{k, 1}, rows{k, 2}); end
fprintf('\nCo obserwować:\n'); for k = 1:numel(s.look), fprintf('  • %s\n', s.look{k}); end
fprintf('Pytania:\n'); for k = 1:numel(s.ask), fprintf('  %d. %s\n', k, s.ask{k}); end
R.scenario = struct('id', id, 'title', tytul, 'rows', {rows});
if ~rysuj, return; end
% zakres wykresów
T = R.ch.t(end); g = [];
for e = R.meta.events
  if strcmp(e.type, 'gust') && ~isempty(regexp(e.text, '^początek podmuchu', 'once')), g = e; break; end
end
tR = [0 T];
if isfield(widok, 'gust') && ~isempty(g), tR = [max(0, g.t + widok.gust(1)), min(T, g.t + widok.gust(2))];
elseif isfield(widok, 't1'), tR = [widok.t0, min(T, widok.t1)];
end
c = widok.cursor;
if ischar(c)
  if strcmp(c, 'apex')
    if ~isempty(R.summary.apex), c = R.summary.apex.t; else, c = T / 2; end
  elseif strncmp(c, 'gust', 4)
    if ~isempty(g), c = g.t + str2double(c(5:end)); else, c = T / 2; end
  end
end
plotResults(R, struct('tRange', tR, 'cursor', min(c, T), 'rosWin', widok.ros));
end
