function v = gustExcitation(shape, tr, sF, N)
%GUSTEXCITATION Pobudzenie ramienia szybkiego przez podmuch narastający w czasie τ_r, względem skoku:
%   |∫₀^τr f′(τ)·e^{−s_F·τ} dτ| (f – kształt narastania: liniowy albo 1 − cos; Simpson, N przedziałów).
if nargin < 4, N = 200; end
h = tr / N; tau = (0:N) * h;
w = 2 * ones(1, N + 1); w(2:2:N) = 4; w([1 end]) = 1;
if strcmp(shape, 'cos'), fp = 0.5 * pi / tr * sin(pi * tau / tr); else, fp = ones(size(tau)) / tr; end
v = abs(sum(w .* fp .* exp(-sF * tau)) * h / 3);
end
