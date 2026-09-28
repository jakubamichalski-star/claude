function o = aeroCoef(P, M, alphaTot)
%AEROCOEF Współczynniki aerodynamiczne w chwili (ℳ, α_tot) – przypisanie wg sekcji 6.
%   o = aeroCoef(P, M, alphaTot) – P z makeParams, alphaTot [rad].
%   κ_p = 1 (pd/V) lub ½ (pd/2V) wg AEP-96 (B.37)–(B.44).
%   Ĉ = f_C·C + f_ΔC·ΔC (B.9), ΔC = σ z tabeli błędów pliku.
%   W ℳ interpolacja liniowa, poza zakresem tabel wartości skrajne; C_mag-m(ℳ, α_tot) – biliniowo.
%   Pola: CD, CLa, CMa, CMq, CNq, Cmagf, Cmagm, Cspin, CNa, CPN, clamp
%   (clamp = −1/+1, gdy ℳ poza zakresem tabel – użyto wartości skrajnych).
c = aeroCoefVec(P, M, alphaTot);
o = struct('CD', c(1), 'CLa', c(2), 'CMa', c(3), 'CMq', c(4), 'CNq', c(5), 'Cmagf', c(6), 'Cmagm', c(7), ...
  'Cspin', c(8), 'CNa', c(9), 'CPN', c(10), 'clamp', c(11));
end
