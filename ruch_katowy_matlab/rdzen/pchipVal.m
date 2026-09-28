function v = pchipVal(pc, x)
%PCHIPVAL Wartość interpolanta PCHIP (pchipMake) w punkcie x; poza zakresem – wartości skrajne.
xs = pc.xs; n = pc.n;
if ~(x > xs(1)), v = pc.ys(1); return; end
if x >= xs(n), v = pc.ys(n); return; end
k = find(xs < x, 1, 'last');
h = pc.h(k); t = (x - xs(k)) / h; t2 = t * t; t3 = t2 * t;
v = (2 * t3 - 3 * t2 + 1) * pc.ys(k) + (t3 - 2 * t2 + t) * h * pc.m(k) ...
  + (-2 * t3 + 3 * t2) * pc.ys(k + 1) + (t3 - t2) * h * pc.m(k + 1);
end
