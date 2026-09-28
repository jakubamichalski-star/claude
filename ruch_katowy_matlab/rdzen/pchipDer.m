function v = pchipDer(pc, x)
%PCHIPDER Pochodna interpolanta PCHIP (pchipMake) w punkcie x; poza zakresem – 0.
xs = pc.xs; n = pc.n;
if ~(x > xs(1)) || x >= xs(n), v = 0; return; end
k = find(xs < x, 1, 'last');
h = pc.h(k); t = (x - xs(k)) / h; t2 = t * t;
v = ((6 * t2 - 6 * t) * pc.ys(k) + (3 * t2 - 4 * t + 1) * h * pc.m(k) ...
   + (-6 * t2 + 6 * t) * pc.ys(k + 1) + (3 * t2 - 2 * t) * h * pc.m(k + 1)) / h;
end
