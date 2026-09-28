function pasmo(ax, x, lo, hi, kolor, alfa)
%PASMO Wypełnione pasmo między lo i hi (odcinki z wartościami skończonymi, bez legendy).
if nargin < 6, alfa = 0.15; end
x = x(:); lo = lo(:); hi = hi(:);
ok = isfinite(x) & isfinite(lo) & isfinite(hi);
d = diff([0; ok; 0]); s = find(d == 1); e = find(d == -1) - 1;
for k = 1:numel(s)
  j = s(k):e(k);
  fill(ax, [x(j); flipud(x(j))], [lo(j); flipud(hi(j))], kolor, 'EdgeColor', 'none', 'FaceAlpha', alfa, 'HandleVisibility', 'off');
end
end
