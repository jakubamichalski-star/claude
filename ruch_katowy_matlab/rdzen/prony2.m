function f = prony2(z, dt)
%PRONY2 Dwa mody zespolone z próbek równoodległych (Prony, najmniejsze kwadraty).
%   f = prony2(z, dt): f.sF, f.sS – wykładniki (F: większe |Im s|), f.KF, f.KS – amplitudy
%   w chwili pierwszej próbki, f.relRes – względna reszta dopasowania.
z = z(:);
u = z(2:end - 1); v = z(1:end - 2); w = z(3:end);
A11 = sum(abs(u) .^ 2); A22 = sum(abs(v) .^ 2);
A12 = sum(conj(u) .* v); b1 = sum(conj(u) .* w); b2 = sum(conj(v) .* w);
A21 = conj(A12);
det_ = A11 * A22 - abs(A12)^2;
a1 = (b1 * A22 - A12 * b2) / det_;
a2 = (b2 * A11 - A21 * b1) / det_;
s = orderRoots(log(quadRoots(-a1, -a2)) / dt);
k = (0:numel(z) - 1).';
eF = exp(s(1) * dt * k); eS = exp(s(2) * dt * k);
G11 = sum(abs(eF) .^ 2); G22 = sum(abs(eS) .^ 2); G12 = sum(conj(eF) .* eS);
r1 = sum(conj(eF) .* z); r2 = sum(conj(eS) .* z);
G21 = conj(G12); dg = G11 * G22 - abs(G12)^2;
KF = (r1 * G22 - G12 * r2) / dg;
KS = (r2 * G11 - G21 * r1) / dg;
res = sum(abs(z - (KF * eF + KS * eS)) .^ 2); nrm = sum(abs(z) .^ 2);
f = struct('sF', s(1), 'sS', s(2), 'KF', KF, 'KS', KS, 'relRes', sqrt(res / nrm));
end
