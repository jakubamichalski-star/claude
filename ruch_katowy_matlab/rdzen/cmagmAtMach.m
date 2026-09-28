function Cfn = cmagmAtMach(P, M)
%CMAGMATMACH C_mag-m(α) przy ustalonej ℳ – szybka ścieżka do uśredniania po fazach (zgodna z aeroCoef).
%   Cfn(alpha) – uchwyt funkcji (alpha [rad], dowolny wymiar).
mach = P.mach; n = numel(mach);
Mc = M; if ~(Mc >= mach(1)), Mc = mach(1); elseif Mc > mach(n), Mc = mach(n); end
i0 = 1 + sum(mach(2:n - 1) < Mc);
w = (Mc - mach(i0)) / (mach(i0 + 1) - mach(i0));
al = P.cnpaAlpha; na = numel(al);
row = (1 - w) * P.cnpaData(i0, :) + w * P.cnpaData(i0 + 1, :);
CYpa = P.coefTab(:, 12);
base = ((1 - w) * CYpa(i0) + w * CYpa(i0 + 1)) * P.dXcg; kf = P.f.Cmagm * P.kp; fdm = P.fd.Cmagm * P.kp;
A = P.A;
Cfn = @evalC;
  function c = evalC(alpha)
    aD = min(max(alpha * (180 / pi), al(1)), al(na));
    ja = ones(size(aD));
    for k = 2:na - 1, ja = ja + (aD > al(k)); end
    wa = (aD - al(ja)) ./ (al(ja + 1) - al(ja));
    c = kf * ((1 - wa) .* row(ja) + wa .* row(ja + 1) + base);
    if fdm, c = c + fdm * sigmaCnpaAt(A, Mc, alpha * (180 / pi)); end
  end
end
