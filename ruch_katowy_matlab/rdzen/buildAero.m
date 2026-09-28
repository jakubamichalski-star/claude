function A = buildAero(P)
%BUILDAERO Model aerodynamiczny: przypisanie kolumn pliku Spin2004 do oznaczeń AEP-96.
%   A = buildAero(P) – P z parseSpin2004. Wszystkie współczynniki na siatce ℳ tabeli sił
%   (kolumny o innej siatce – interpolacja liniowa).
if ~isfield(P.tables, 'forces') || ~isfield(P.tables, 'moments')
  error('W pliku brakuje tabeli sił lub momentów – nie mogę zbudować modelu.');
end
if isempty(P.geom), error('W pliku brakuje geometrii (Ref Dia, Ref CG).'); end
F = P.tables.forces; Mo = P.tables.moments;
if isfield(P.tables, 'cdcl'), DL = P.tables.cdcl; else, DL = []; end
warn = {};
mach = F.rows(:, 1);
A = struct('mach', mach);
for n = {'CX0', 'CX2', 'CNa', 'CPN', 'CYpa'}, A.(n{1}) = col(F, n{1}, true); end
for n = {'CX4', 'CNa3', 'CXf', 'CXb', 'CNq'}, A.(n{1}) = col(F, n{1}, false); end
for n = {'Cma', 'Cmq', 'Clp'}, A.(n{1}) = col(Mo, n{1}, true); end
for n = {'Cma3', 'Cma5', 'Cmq2', 'Cld'}, A.(n{1}) = col(Mo, n{1}, false); end
for d = [1 3 5], A.(sprintf('Cnpa%d', d)) = col(Mo, sprintf('Cnpa@%d', d), false); end
if ~isempty(DL)
  A.CDo = col(DL, 'CDo', true); A.CD2 = col(DL, 'CD2', true); A.CLa = col(DL, 'CLa', true); A.CLa3 = col(DL, 'CLa3', false);
  A.cdclSource = 'tabela „CD and CL Convention”';
else   % przeliczenie jak w Spin2004 (cos α ≈ 1)
  A.CDo = A.CX0; A.CD2 = A.CX2 + A.CNa; A.CLa = A.CNa - A.CX0; A.CLa3 = A.CNa3 - A.CX2;
  A.cdclSource = 'przeliczenie z CX, CN (cos α ≈ 1)';
  warn{end + 1} = 'Brak tabeli CD/CL – C_D i C_L przeliczone z CX, CN.';
end
% C_mag-m(ℳ, α): tabela 2D, a gdy jej brak – Cnpa@1, @3, @5
if isfield(P.tables, 'cnpa2d') && ~isempty(P.tables.cnpa2d.rows)
  T2 = P.tables.cnpa2d; al = T2.alpha(:).'; na = numel(al); m2 = T2.rows(:, 1);
  data = zeros(numel(mach), na);
  for j = 1:na
    for i = 1:numel(mach), data(i, j) = interp1c(m2, T2.rows(:, j + 1), mach(i)); end
  end
  A.cnpaAlpha = al; A.cnpaData = data; A.cnpaSource = 'tabela Cnpa(ℳ, α), interpolacja biliniowa';
else
  A.cnpaAlpha = [1 3 5];
  A.cnpaData = [A.Cnpa1, A.Cnpa3, A.Cnpa5];
  A.cnpaSource = 'Cnpa@1, @3, @5 (brak tabeli 2D; poza 1–5° wartości skrajne)';
  warn{end + 1} = 'Brak tabeli Cnpa(ℳ, α) – użyto Cnpa@1, @3, @5.';
end
g = P.geom;
A.refDia_mm = g.refDia; A.refCG_mm = g.refCG;
A.d = g.refDia / 1000; A.S = pi * A.d^2 / 4;
A.xcgRef_cal = g.refCG / g.refDia;
A.geom = g;
A.errors = P.errors;
A.warn = warn;

  function v = col(tab, name, required)
    k = find(strcmp(tab.cols, name), 1);
    if isempty(k)
      if required, warn{end + 1} = sprintf('Brak kolumny %s – przyjęto 0.', name); end
      v = zeros(numel(mach), 1); return;
    end
    m = tab.rows(:, 1); v = tab.rows(:, k);
    same = numel(m) == numel(mach) && all(abs(m - mach) < 1e-9);
    if ~same
      warn{end + 1} = sprintf('Kolumna %s: inna siatka ℳ – przeinterpolowano liniowo.', name);
      v = arrayfun(@(x) interp1c(m, tab.rows(:, k), x), mach);
    end
  end
end

function y = interp1c(xs, ys, x)
% interpolacja liniowa, poza zakresem – wartości skrajne
if x <= xs(1), y = ys(1); return; end
if x >= xs(end), y = ys(end); return; end
i = find(xs < x, 1, 'last');
y = ys(i) + (x - xs(i)) / (xs(i + 1) - xs(i)) * (ys(i + 1) - ys(i));
end
