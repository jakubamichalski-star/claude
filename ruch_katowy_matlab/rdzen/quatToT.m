function T = quatToT(l)
%QUATTOT Macierz T (wiersze: x̄, ȳ, z̄ w składowych układu I) z kwaternionu λ̄.
%   T – wektor 1×9 (macierz 3×3 zapisana wierszami: T(1..3) – wiersz x̄ itd.).
l0 = l(1); l1 = l(2); l2 = l(3); l3 = l(4);
k = 1 / (l0 * l0 + l1 * l1 + l2 * l2 + l3 * l3);
T = [(l0 * l0 + l1 * l1 - l2 * l2 - l3 * l3) * k, 2 * (l1 * l2 + l0 * l3) * k, 2 * (l1 * l3 - l0 * l2) * k, ...
     2 * (l1 * l2 - l0 * l3) * k, (l0 * l0 - l1 * l1 + l2 * l2 - l3 * l3) * k, 2 * (l2 * l3 + l0 * l1) * k, ...
     2 * (l1 * l3 + l0 * l2) * k, 2 * (l2 * l3 - l0 * l1) * k, (l0 * l0 - l1 * l1 - l2 * l2 + l3 * l3) * k];
end
