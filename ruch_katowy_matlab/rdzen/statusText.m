function s = statusText(k)
%STATUSTEXT Opis klasy stanu (statusClass).
switch k
  case 'ok', s = 'STATECZNY';
  case 'gyro', s = 'NIESTATECZNY GIROSKOPOWO (S_g < 1)';
  case 'dynF', s = 'NIESTATECZNY DYNAMICZNIE – składowa szybka (σ_F > 0)';
  case 'dynS', s = 'NIESTATECZNY DYNAMICZNIE – składowa wolna (σ_S > 0)';
  case 'dynFS', s = 'NIESTATECZNY DYNAMICZNIE – obie składowe';
  case 'weak', s = 'SŁABA NADĄŻNOŚĆ (|ᾱ_e| > próg)';
  case 'limit', s = 'CYKL GRANICZNY – niestateczny przy α_tot → 0, wzrost zatrzymany przez nieliniowy moment Magnusa';
  otherwise, s = k;
end
end
