function P = parseSpin2004(text, fileName)
%PARSESPIN2004 Parser pliku wyników Spin2004 (ASCII, kolumny o stałej szerokości).
%   P = parseSpin2004(text, fileName)
%   Tabele rozpoznawane po nagłówkach kolumn, nie po numerach wierszy.
%   P.tables.<rodzaj>: cols (cellstr), rows (macierz), raw (cellstr – zapis w pliku),
%   alpha (kąty tabeli 2D [°] albo []), line (numer wiersza nagłówka).
%   Rodzaje: forces, moments, cdcl, mcdrag, cnpa2d, cnp2d.
%   P.errors – tabele szacunków błędów (regime, m1, m2, sigma, n).
if nargin < 2, fileName = ''; end
text = strrep(strrep(char(text), sprintf('\r\n'), sprintf('\n')), sprintf('\r'), sprintf('\n'));
lines = regexp(text, '\n', 'split');
P = struct('fileName', fileName, 'project', '', 'date', '', 'version', '', 'geom', [], 'tables', struct(), ...
  'errors', struct('regime', {}, 'm1', {}, 'm2', {}, 'sigma', {}, 'n', {}), 'xcgHeaderCal', NaN);
P.warnings = {};
for k = 1:numel(lines)
  s = strtrim(lines{k});
  if isempty(P.project) && ~isempty(regexpi(s, '\.pr3(?![A-Za-z0-9_])', 'once')), P.project = regexprep(s, '\s+-\s+\d+\s*$', ''); end
  if isempty(P.date) && ~isempty(regexp(s, '^\d{1,2}\.\d{1,2}\.\d{4}\s+\d{1,2}:\d{2}', 'once')), P.date = s; end
  if isempty(P.version) && ~isempty(regexpi(s, '^Spin\s*\d+.*Version', 'once')), P.version = s; end
end
% geometria: wiersz liczb zakończony „mm”, pod nim wiersz zakończony „Caliber”
GEOM_NAMES = {'refDia', 'refCG', 'totalLength', 'ogiveLength', 'ogiveRadius', 'boattailLength', ...
  'boattailDia', 'meplatDia', 'bandLength', 'bandDia', 'boomLength', 'boomDia'};
for i = 1:numel(lines)
  tk = splitTok(lines{i});
  if numel(tk) >= 3 && strcmp(tk{end}, 'mm') && allNum(tk(1:end - 1))
    mm = str2double(tk(1:end - 1));
    if i < numel(lines), tk2 = splitTok(lines{i + 1}); else, tk2 = {}; end
    cal = [];
    if ~isempty(tk2) && strcmp(tk2{end}, 'Caliber') && allNum(tk2(1:end - 1)), cal = str2double(tk2(1:end - 1)); end
    if numel(mm) ~= numel(GEOM_NAMES)
      P.warnings{end + 1} = sprintf('Geometria: oczekiwano %d liczb, jest %d.', numel(GEOM_NAMES), numel(mm));
    end
    g = struct('cal', struct());
    for k = 1:numel(GEOM_NAMES)
      if k <= numel(mm), g.(GEOM_NAMES{k}) = mm(k); else, g.(GEOM_NAMES{k}) = NaN; end
      if ~isempty(cal) && k <= numel(cal), g.cal.(GEOM_NAMES{k}) = cal(k); end
    end
    P.geom = g; break;
  end
end
if isempty(P.geom), P.warnings{end + 1} = 'Nie znaleziono wiersza geometrii (liczby zakończone „mm”).'; end

i = 1;
while i <= numel(lines)
  L = lines{i};
  mX = regexpi(L, 'Moments at Baseline CG\s*=\s*([\d.]+)\s*Calibers', 'tokens', 'once');
  if ~isempty(mX) && isnan(P.xcgHeaderCal), P.xcgHeaderCal = str2double(mX{1}); end
  [tt, tc] = tokensWithPos(L);
  if isempty(tt) || ~strcmp(tt{1}, 'Mach'), i = i + 1; continue; end
  if ~isempty(regexpi(L, 'Regime', 'once'))
    [P, i] = parseErrorTable(lines, i, P); i = i + 1; continue;
  end
  cols = tt; alpha = [];
  if numel(cols) > 1 && allNum(cols(2:end))
    alpha = str2double(cols(2:end));            % tabela 2D w funkcji (ℳ, α[°])
    title = '';
    for k = i - 1:-1:max(1, i - 4)
      if ~isempty(regexpi(lines{k}, 'Derivative\s*-\s*Cnpa|Coefficient\s*-\s*Cnp(?![A-Za-z0-9_])', 'once')), title = lines{k}; break; end
    end
    if ~isempty(regexpi(title, 'Derivative', 'once')), kind = 'cnpa2d';
    elseif ~isempty(regexpi(title, 'Coefficient', 'once')), kind = 'cnp2d';
    else, kind = sprintf('table2d_%d', i);
    end
  else
    % powtórzone nazwy kolumn (np. trzy „CXo” w tabeli McDrag) rozróżnia wiersz nad nagłówkiem
    [~, ~, ic] = unique(cols); cnt = accumarray(ic(:), 1);
    aboveT = {}; aboveC = [];
    for k = i - 1:-1:1
      if ~isempty(strtrim(lines{k})), [aboveT, aboveC] = tokensWithPos(lines{k}); break; end
    end
    for k = 1:numel(tt)
      if cnt(ic(k)) > 1 && ~isempty(aboveT)
        [~, jb] = min(abs(aboveC - tc(k)));
        cols{k} = [aboveT{jb} ':' tt{k}];
      end
    end
    if any(strcmp(cols, 'CYpa')) || (any(strcmp(cols, 'CX0')) && any(strcmp(cols, 'CNa')) && any(strcmp(cols, 'CPN')))
      kind = 'forces';
    elseif any(strcmp(cols, 'Cma')) && any(strcmp(cols, 'Cmq')), kind = 'moments';
    elseif any(strcmp(cols, 'CDo')) && any(strcmp(cols, 'CLa')), kind = 'cdcl';
    elseif any(strcmp(cols, 'CXH')) || any(strcmp(cols, 'CXSF')), kind = 'mcdrag';
    else, kind = sprintf('table_%d', i);
    end
  end
  rows = zeros(0, numel(cols)); raw = cell(0, numel(cols));
  j = i + 1;
  while j <= numel(lines)
    s = strtrim(lines{j});
    if isempty(s)
      if ~isempty(rows), break; end
      j = j + 1; continue;
    end
    if ~isempty(regexp(s, '^-{3,}$', 'once')), j = j + 1; continue; end    % linie kreskowe
    tk = splitTok(s);
    if all(~cellfun(@isempty, regexpi(tk, '^Spinner$', 'once'))), j = j + 1; continue; end   % adnotacje „Spinner”
    if s(1) == '{' || ~allNum(tk), break; end
    if numel(tk) ~= numel(cols)
      P.warnings{end + 1} = sprintf('Tabela %s: wiersz %d ma %d liczb zamiast %d – pominięty.', kind, j, numel(tk), numel(cols));
      j = j + 1; continue;
    end
    rows(end + 1, :) = str2double(tk); raw(end + 1, :) = tk; %#ok<AGROW>
    j = j + 1;
  end
  if isfield(P.tables, kind)
    P.warnings{end + 1} = sprintf('Tabela %s wystąpiła ponownie (wiersz %d) – użyto pierwszej.', kind, i);
  else
    P.tables.(kind) = struct('cols', {cols}, 'rows', rows, 'raw', {raw}, 'alpha', alpha, 'line', i);
  end
  i = j;
end
need = {'forces', 'moments'};
for k = 1:2
  if ~isfield(P.tables, need{k}), P.warnings{end + 1} = sprintf('Brak tabeli „%s”.', need{k}); end
end
end

% ---------------------------------------------------------------------
function [P, iEnd] = parseErrorTable(lines, i, P)
% Tabele szacunków błędów: „Subsonic (0.010 - 0.850) = n σ n σ …”.
hdr = splitTok(lines{i});
ir = find(strcmp(hdr, 'Regime'), 1);
names = hdr(ir + 1:end); names = names(~strcmp(names, 'n'));
j = i + 1;
if j <= numel(lines) && ~isempty(regexpi(lines{j}, 'deg', 'once'))
  dt = splitTok(lines{j});
  for k = 1:numel(names)
    if k <= numel(dt)
      v = str2double(regexp(dt{k}, '^[+-]?(\d+\.?\d*|\.\d+)', 'match', 'once'));
      if isfinite(v), names{k} = sprintf('%s@%g', names{k}, v); end
    end
  end
  j = j + 1;
end
got = 0;
while j <= numel(lines)
  s = strtrim(lines{j});
  if isempty(s)
    if got, break; end
    j = j + 1; continue;
  end
  if ~isempty(regexp(s, '^-{3,}$', 'once')), j = j + 1; continue; end
  m = regexp(s, '^(.*?)\(\s*([\d.]+)\s*-\s*([\d.]+)\s*\)\s*=\s*(.+)$', 'tokens', 'once');
  if isempty(m), break; end
  nums = str2double(splitTok(m{4}));
  m1 = str2double(m{2}); m2 = str2double(m{3});
  ie = find([P.errors.m1] == m1 & [P.errors.m2] == m2, 1);
  if isempty(ie)
    ie = numel(P.errors) + 1;
    P.errors(ie).regime = strtrim(m{1}); P.errors(ie).m1 = m1; P.errors(ie).m2 = m2;
    P.errors(ie).sigma = struct(); P.errors(ie).n = struct();
  end
  for k = 1:numel(names)
    key = sigmaKey(names{k});
    if 2 * k <= numel(nums)
      P.errors(ie).n.(key) = nums(2 * k - 1); P.errors(ie).sigma.(key) = nums(2 * k);
    end
  end
  got = got + 1; j = j + 1;
end
if ~isempty(P.errors)
  [~, o] = sort([P.errors.m1]); P.errors = P.errors(o);
end
iEnd = j - 1;
end

function tk = splitTok(s)
tk = regexp(strtrim(s), '\s+', 'split');
if numel(tk) == 1 && isempty(tk{1}), tk = {}; end
end

function [t, c] = tokensWithPos(line)
[t, st] = regexp(line, '\S+', 'match', 'start');
c = st - 1 + cellfun(@numel, t) / 2;     % środek słowa (indeks od 0, jak w oryginale)
end

function ok = allNum(tk)
ok = all(~cellfun(@isempty, regexp(tk, '^[+-]?(\d+\.?\d*|\.\d+)([eE][+-]?\d+)?$', 'once')));
end
