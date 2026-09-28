function v = sigmaAt(A, key, M)
%SIGMAAT Odchylenie standardowe współczynnika z pliku dla danego ℳ (najbliższy zakres).
%   key – nazwa jak w tabeli błędów ('CX', 'CNa', 'Cma', 'CPN', 'Cmq', 'Clp', 'Cnpa@1', …).
E = A.errors; v = 0;
if isempty(E), return; end
best = 1; bd = Inf;
for k = 1:numel(E)
  if M < E(k).m1, d = E(k).m1 - M; elseif M > E(k).m2, d = M - E(k).m2; else, d = 0; end
  if d < bd, bd = d; best = k; end
end
f = sigmaKey(key);
if isfield(E(best).sigma, f) && isfinite(E(best).sigma.(f)), v = E(best).sigma.(f); end
end
