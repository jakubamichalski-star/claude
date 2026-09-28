function out = deepMerge(a, b)
%DEEPMERGE Scalenie struktur: pola b nadpisują pola a, podstruktury (skalarne) scalane rekurencyjnie.
%   Tablice (także tablice struktur, np. met.zones) są zastępowane w całości.
out = a;
if ~isScalarStruct(b), return; end
fn = fieldnames(b);
for k = 1:numel(fn)
  f = fn{k};
  if isScalarStruct(b.(f)) && isfield(a, f) && isScalarStruct(a.(f))
    out.(f) = deepMerge(a.(f), b.(f));
  else
    out.(f) = b.(f);
  end
end
end

function t = isScalarStruct(x)
t = isstruct(x) && isscalar(x);
end
