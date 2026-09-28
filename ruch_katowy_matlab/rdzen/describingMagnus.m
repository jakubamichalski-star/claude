function [CF, CS] = describingMagnus(Cfn, e, aF, aS, N)
%DESCRIBINGMAGNUS Funkcja opisująca (quasi-liniowa teoria Murphy'ego) członu f(ξ) = C(|ξ|)·ξ.
%   Dla ruchu ξ = e + a_F·e^{iφ} + a_S·e^{iψ} o niezależnych, jednostajnie rozłożonych fazach:
%   C_F = Re⟨f·e^{−iφ}⟩/a_F,  C_S = Re⟨f·e^{−iψ}⟩/a_S  (trapezy na siatce N×N, e = |ξ_e|).
%   Przypadki graniczne: samo ramię bez trymu → C(a) (sieczna); małe ramię wokół tła
%   o module b → ⟨C(b) + (b/2)·C′(b)⟩, czyli średnia siecznej i stycznej.
if nargin < 5 || isempty(N), N = 20; end
aFe = max(aF, 1e-4); aSe = max(aS, 1e-4);
ph = 2 * pi * (0:N - 1) / N; cs = cos(ph); sn = sin(ph);
[CSk, CSj] = meshgrid(cs, cs); [SNk, SNj] = meshgrid(sn, sn);    % j – faza ramienia F, k – S
xr = e + aFe * CSj + aSe * CSk; xi = aFe * SNj + aSe * SNk;
C = Cfn(sqrt(xr .^ 2 + xi .^ 2));
sF = sum(sum(C .* (xr .* CSj + xi .* SNj)));
sS = sum(sum(C .* (xr .* CSk + xi .* SNk)));
CF = sF / (N * N * aFe); CS = sS / (N * N * aSe);
end
