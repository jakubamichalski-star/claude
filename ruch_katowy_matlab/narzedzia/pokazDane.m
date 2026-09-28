function pokazDane(Spin, A, V)
%POKAZDANE Opis wczytanego pliku Spin2004: nagłówek, geometria, rozpoznane tabele, weryfikacja wczytania.
%   pokazDane(Spin, A, V) – Spin z parseSpin2004, A z buildAero, V z verifySpin.
fprintf('\nPlik Spin2004: %s\n  projekt: %s; %s; %s\n', Spin.fileName, Spin.project, Spin.date, Spin.version);
g = Spin.geom;
fprintf('  geometria: d = %s mm, X_CG,ref = %s mm (%s kal.), długość %s mm\n', fmtPL(g.refDia, 2), fmtPL(g.refCG, 2), ...
  fmtPL(A.xcgRef_cal, 4), fmtPL(g.totalLength, 2));
tn = fieldnames(Spin.tables);
for k = 1:numel(tn)
  t = Spin.tables.(tn{k});
  fprintf('  tabela %-8s wiersz %4d: %2d kolumn, %2d wierszy (ℳ = %g…%g)\n', tn{k}, t.line, numel(t.cols), size(t.rows, 1), t.rows(1, 1), t.rows(end, 1));
end
fprintf('  C_D, C_L: %s\n  C_mag-m: %s\n', A.cdclSource, A.cnpaSource);
for k = 1:numel(Spin.errors)
  e = Spin.errors(k);
  fprintf('  błędy %-10s (ℳ %g–%g): σ_Cma = %g, σ_Cmq = %g, σ_Cnpa@1 = %g\n', e.regime, e.m1, e.m2, getf(e.sigma, 'Cma'), getf(e.sigma, 'Cmq'), getf(e.sigma, 'Cnpa_1'));
end
w = [Spin.warnings, A.warn];
for k = 1:numel(w), fprintf('  UWAGA: %s\n', w{k}); end
if nargin > 2 && ~isempty(V)
  fprintf('\nWeryfikacja wczytania (tolerancje z zaokrągleń zapisu w pliku)\n');
  for k = 1:numel(V)
    if V(k).pass, s = 'OK  '; else, s = 'BŁĄD'; end
    fprintf('  %-4s %s %s\n        %s\n', V(k).id, s, V(k).name, V(k).info);
    for q = 1:numel(V(k).fails), fprintf('        %s\n', V(k).fails{q}); end
  end
end
end

function v = getf(s, f)
if isfield(s, f), v = s.(f); else, v = NaN; end
end
