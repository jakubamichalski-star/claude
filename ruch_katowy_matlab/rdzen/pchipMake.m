function pc = pchipMake(xs, ys)
%PCHIPMAKE Interpolacja PCHIP (Fritsch–Carlson, warunki brzegowe jak w MATLAB-owym pchip).
%   pc = pchipMake(xs, ys) – struktura do pchipVal (f(x)) i pchipDer (df/dx).
%   Poza zakresem węzłów: wartości skrajne i pochodna 0.
xs = xs(:).'; ys = ys(:).'; n = numel(xs);
pc = struct('xs', xs, 'ys', ys, 'n', n, 'h', [], 'm', zeros(1, n));
if n == 1, return; end
h = diff(xs); d = diff(ys) ./ h; m = zeros(1, n);
if n == 2
  m(:) = d(1);
else
  for k = 2:n - 1
    if d(k - 1) * d(k) <= 0
      m(k) = 0;
    else
      w1 = 2 * h(k) + h(k - 1); w2 = h(k) + 2 * h(k - 1);
      m(k) = (w1 + w2) / (w1 / d(k - 1) + w2 / d(k));
    end
  end
  m(1) = endSlope(h(1), h(2), d(1), d(2));
  m(n) = endSlope(h(n - 1), h(n - 2), d(n - 1), d(n - 2));
end
pc.h = h; pc.m = m;
end

function v = endSlope(h0, h1, d0, d1)
v = ((2 * h0 + h1) * d0 - h0 * d1) / (h0 + h1);
if sign(v) ~= sign(d0)
  v = 0;
elseif sign(d0) ~= sign(d1) && abs(v) > abs(3 * d0)
  v = 3 * d0;
end
end
