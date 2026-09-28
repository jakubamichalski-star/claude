function c = kolorStanu(k)
%KOLORSTANU Kolor klasy stanu (statusClass) na pasku stanu.
switch k
  case 'ok', c = [60 141 90] / 255;
  case 'weak', c = [217 169 58] / 255;
  case 'limit', c = [196 105 43] / 255;
  case 'dynF', c = [181 58 58] / 255;
  case 'dynS', c = [150 45 75] / 255;
  case 'dynFS', c = [130 30 40] / 255;
  case 'gyro', c = [110 26 42] / 255;
  otherwise, c = [205 211 219] / 255;
end
end
