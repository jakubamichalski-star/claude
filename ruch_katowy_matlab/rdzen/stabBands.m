function B = stabBands(P, M, ~, L)
%STABBANDS Pasma 1σ dla S_g i S_d (α → 0): propagacja liniowa odchyleń z pliku (RSS).
%   B.SgSig, B.SdSig; L z linPoint.
A = P.A; kp = P.kp; md2 = P.m * P.d * P.d;
sCMa = sigmaAt(A, 'Cma', M); sCLa = sigmaAt(A, 'CNa', M); sCD = sigmaAt(A, 'CX', M);
sCMq = kp * sigmaAt(A, 'Cmq', M); sCm = kp * sigmaCnpaAt(A, M, 0);
B.SgSig = abs(L.Sg) * sCMa / abs(L.CMa);
Sd = @(CLa, CD0, CMq, Cm) 2 * (CLa + (md2 / P.Ix) * Cm) / (CLa - CD0 - (md2 / P.IT) * CMq);
S0 = Sd(L.CLa, L.CD0, L.CMq, L.Cmagm0);
d1 = Sd(L.CLa + sCLa, L.CD0, L.CMq, L.Cmagm0) - S0; d2 = Sd(L.CLa, L.CD0 + sCD, L.CMq, L.Cmagm0) - S0;
d3 = Sd(L.CLa, L.CD0, L.CMq + sCMq, L.Cmagm0) - S0; d4 = Sd(L.CLa, L.CD0, L.CMq, L.Cmagm0 + sCm) - S0;
B.SdSig = sqrt(d1 * d1 + d2 * d2 + d3 * d3 + d4 * d4);
end
