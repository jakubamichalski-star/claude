function v = scenMemo(key, val)
%SCENMEMO Pamięć wyników między scenariuszami (np. |K_F| wariantu skokowego (f) do porównania z łagodnym).
%   scenMemo(key, val) – zapis; v = scenMemo(key) – odczyt (NaN, gdy brak).
persistent M
if isempty(M), M = struct(); end
if nargin > 1, M.(key) = val; end
if isfield(M, key), v = M.(key); else, v = NaN; end
end
