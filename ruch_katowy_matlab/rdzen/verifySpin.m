function res = verifySpin(P)
%VERIFYSPIN Testy wczytania pliku Spin2004 (sekcja 6 specyfikacji).
%   res = verifySpin(P) – tablica struktur: id, name, pass, info, maxRatio, n, fails.
%   Tolerancja z zaokrągleń w pliku: suma połówek jednostki ostatniej cyfry każdej użytej wartości.
DEG = pi / 180;
T = P.tables;
F = tab(T, 'forces'); Mo = tab(T, 'moments'); DL = tab(T, 'cdcl'); MD = tab(T, 'mcdrag');
A2 = tab(T, 'cnpa2d'); AC = tab(T, 'cnp2d');
g = P.geom;
res = struct('id', {}, 'name', {}, 'pass', {}, 'info', {}, 'maxRatio', {}, 'n', {}, 'fails', {});
% V1: jednakowa siatka ℳ
all_ = {'siły', F; 'momenty', Mo; 'Cnpa(ℳ,α)', A2; 'Cnp(ℳ,α)', AC; 'CD/CL', DL; 'McDrag', MD};
if isempty(F), ref = []; else, ref = F.rows(:, 1); end
bad = {}; nt = 0;
for k = 1:size(all_, 1)
  t = all_{k, 2}; if isempty(t), continue; end
  nt = nt + 1; m = t.rows(:, 1);
  if isempty(ref) || numel(m) ~= numel(ref) || any(abs(m - ref) > 1e-9), bad{end + 1} = all_{k, 1}; end %#ok<AGROW>
end
if ~isempty(bad)
  info = ['inna siatka: ' strjoin(bad, ', ')];
elseif ~isempty(ref)
  info = sprintf('%d tabel, %d punktów, ℳ = %g…%g', nt, numel(ref), ref(1), ref(end));
else
  info = sprintf('%d tabel, 0 punktów', nt);
end
res(end + 1) = struct('id', 'V1', 'name', 'Jednakowa siatka ℳ we wszystkich tabelach', 'pass', isempty(bad), ...
  'info', info, 'maxRatio', NaN, 'n', nt, 'fails', {{}});

if ~isempty(g), xcgRef = g.refCG / g.refDia; xcgTol = 0.005 / g.refDia; else, xcgRef = NaN; xcgTol = 0; end
% V2: Cma ≈ CNa·(X_CG,ref − CPN)
d = []; tl = []; wh = {};
if ~isempty(F) && ~isempty(Mo) && ~isempty(g)
  kN = tcol(F, 'CNa'); kP = tcol(F, 'CPN'); kM = tcol(Mo, 'Cma');
  for r = 1:min(size(F.rows, 1), size(Mo.rows, 1))
    CNa = F.rows(r, kN); CPN = F.rows(r, kP); Cma = Mo.rows(r, kM);
    tl(end + 1) = halfUlp(Mo.raw{r, kM}) + halfUlp(F.raw{r, kN}) * abs(xcgRef - CPN) + CNa * (halfUlp(F.raw{r, kP}) + xcgTol); %#ok<AGROW>
    d(end + 1) = Cma - CNa * (xcgRef - CPN); wh{end + 1} = sprintf('ℳ = %g', F.rows(r, 1)); %#ok<AGROW>
  end
end
res(end + 1) = cmp('V2', 'Cma ≈ CNa·(X_CG,ref − CPN)  [tabele siły/momenty]', d, tl, wh);
% V3: Cnp ≈ Cnpa·sin α
d = []; tl = []; wh = {};
if ~isempty(A2) && ~isempty(AC)
  for r = 1:min(size(A2.rows, 1), size(AC.rows, 1))
    for k = 2:numel(A2.cols)
      a = A2.alpha(k - 1) * DEG;
      tl(end + 1) = halfUlp(AC.raw{r, k}) + halfUlp(A2.raw{r, k}) * sin(a); %#ok<AGROW>
      d(end + 1) = AC.rows(r, k) - A2.rows(r, k) * sin(a); %#ok<AGROW>
      wh{end + 1} = sprintf('ℳ = %g, α = %g°', A2.rows(r, 1), A2.alpha(k - 1)); %#ok<AGROW>
    end
  end
end
res(end + 1) = cmp('V3', 'Cnp ≈ Cnpa·sin α  [tabele 2D]', d, tl, wh);
% V4a: CDo = CXo
[d, tl, wh] = combo(DL, {'CDo', 1; 'CXo', -1}, DL);
res(end + 1) = cmp('V4a', 'CDo = CXo  [tabela CD/CL]', d, tl, wh);
% V4b: CD2 ≈ CX2 + CNa
d = []; tl = []; wh = {};
if ~isempty(DL) && ~isempty(F)
  a = tcol(DL, 'CD2'); b = tcol(F, 'CX2'); c = tcol(F, 'CNa');
  for r = 1:min(size(DL.rows, 1), size(F.rows, 1))
    d(end + 1) = DL.rows(r, a) - (F.rows(r, b) + F.rows(r, c)); %#ok<AGROW>
    tl(end + 1) = halfUlp(DL.raw{r, a}) + halfUlp(F.raw{r, b}) + halfUlp(F.raw{r, c}); %#ok<AGROW>
    wh{end + 1} = sprintf('ℳ = %g', DL.rows(r, 1)); %#ok<AGROW>
  end
end
res(end + 1) = cmp('V4b', 'CD2 ≈ CX2 + CNa  [CD/CL vs siły]', d, tl, wh);
% V4c: CLa ≈ CNa − CXo
[d, tl, wh] = combo(DL, {'CLa', 1; 'CNa', -1; 'CXo', 1}, DL);
res(end + 1) = cmp('V4c', 'CLa ≈ CNa − CXo  [tabela CD/CL]', d, tl, wh);
% V4d: CX0 = CXf + CXb
[d, tl, wh] = combo(F, {'CX0', 1; 'CXf', -1; 'CXb', -1}, F);
res(end + 1) = cmp('V4d', 'CX0 = CXf + CXb  [tabela sił]', d, tl, wh);
% V5: Cnpa@1, @3, @5 = kolumny 1°, 3°, 5° tabeli Cnpa
d = []; tl = []; wh = {};
if ~isempty(Mo) && ~isempty(A2)
  for dd = [1 3 5]
    km = tcol(Mo, sprintf('Cnpa@%d', dd)); ka = find(A2.alpha == dd, 1) + 1;
    if isempty(km) || isempty(ka), continue; end
    for r = 1:min(size(Mo.rows, 1), size(A2.rows, 1))
      d(end + 1) = Mo.rows(r, km) - A2.rows(r, ka); %#ok<AGROW>
      tl(end + 1) = halfUlp(Mo.raw{r, km}) + halfUlp(A2.raw{r, ka}); %#ok<AGROW>
      wh{end + 1} = sprintf('ℳ = %g, α = %d°', Mo.rows(r, 1), dd); %#ok<AGROW>
    end
  end
end
res(end + 1) = cmp('V5', 'Cnpa@1, @3, @5 = kolumny 1°, 3°, 5° tabeli Cnpa', d, tl, wh);
% V6: suma składników McDrag = CXo McDrag
[d, tl, wh] = combo(MD, {'McDrag:CXo', 1; 'CXH', -1; 'CXSF', -1; 'CXBND', -1; 'CXBT', -1; 'CXB', -1}, MD);
res(end + 1) = cmp('V6', 'Suma składników McDrag = CXo McDrag (dodatkowy)', d, tl, wh);
end

% ---------------------------------------------------------------------
function t = tab(T, name)
if isfield(T, name), t = T.(name); else, t = []; end
end

function k = tcol(t, name)
if isempty(t), k = []; else, k = find(strcmp(t.cols, name), 1); end
end

function [d, tl, wh] = combo(t, spec, ~)
% różnica Σ znak·kolumna w każdym wierszu tabeli (spec: {nazwa, znak}); pusta, gdy brak kolumny
d = []; tl = []; wh = {};
if isempty(t), return; end
ks = zeros(1, size(spec, 1));
for k = 1:size(spec, 1)
  c = tcol(t, spec{k, 1}); if isempty(c), return; end
  ks(k) = c;
end
sg = [spec{:, 2}];
for r = 1:size(t.rows, 1)
  d(end + 1) = sum(sg .* t.rows(r, ks)); %#ok<AGROW>
  tl(end + 1) = sum(cellfun(@halfUlp, t.raw(r, ks))); %#ok<AGROW>
  wh{end + 1} = sprintf('ℳ = %g', t.rows(r, 1)); %#ok<AGROW>
end
end

function r = cmp(id, name, d, tl, wh)
% porównanie z tolerancją z zaokrągleń
n = numel(d); ratio = abs(d) ./ tl;
if n == 0
  r = struct('id', id, 'name', name, 'pass', false, 'info', 'brak danych do porównania', 'maxRatio', 0, 'n', 0, 'fails', {{}});
  return;
end
[maxR, iw] = max(ratio);
fi = find(ratio > 1);
fails = arrayfun(@(k) sprintf('%s: różnica %.3g > tol. %.2g', wh{k}, d(k), tl(k)), fi(1:min(8, end)), 'UniformOutput', false);
info = sprintf('maks. |różnica| = %.2g (%s), %.0f %% tolerancji', abs(d(iw)), wh{iw}, 100 * maxR);
r = struct('id', id, 'name', name, 'pass', isempty(fi), 'info', info, 'maxRatio', maxR, 'n', n, 'fails', {fails});
end

function v = halfUlp(tok)
% ½ jednostki ostatniej cyfry zapisu – podstawa tolerancji testów wczytania
m = regexp(tok, '\.(\d*)', 'tokens', 'once');
if isempty(m), v = 0.5; else, v = 0.5 * 10^(-numel(m{1})); end
end
