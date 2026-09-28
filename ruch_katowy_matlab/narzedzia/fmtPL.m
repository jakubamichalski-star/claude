function s = fmtPL(x, d)
%FMTPL Liczba z przecinkiem dziesiętnym (zapis polski), d miejsc po przecinku; '–' dla NaN/Inf.
if ~isfinite(x), s = '–'; return; end
if abs(x) < 0.5 * 10^(-d), x = 0; end
s = strrep(sprintf('%.*f', d, x), '.', ',');
end
