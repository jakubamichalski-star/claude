function [T, p, rho, a] = atmosphere(h)
%ATMOSPHERE Atmosfera wzorcowa ISO 2533 (ICAO), jak w STANAG 4355.
%   [T, p, rho, a] = atmosphere(h)
%   h – wysokość geometryczna [m]; warstwy w wysokości geopotencjalnej H.
%   T [K], p [Pa], rho [kg/m^3], a [m/s].
persistent LAY
if isempty(LAY)
  LAY = [ 0     288.15 -0.0065 101325.0   % H_b [m], T_b [K], L_b [K/m], p_b [Pa]
          11000 216.65  0.0     22632.06
          20000 216.65  0.001    5474.889
          32000 228.65  0.0028    868.0187
          47000 270.65  0.0       110.9063
          51000 270.65 -0.0028     66.93887
          71000 214.65 -0.002       3.956420 ];
end
G0 = 9.80665; R_AIR = 287.05287; GAMMA_AIR = 1.4; R_EARTH_ISO = 6356766;
H = R_EARTH_ISO * h / (R_EARTH_ISO + h);
i = 1;
while i < size(LAY, 1) && H >= LAY(i + 1, 1), i = i + 1; end
L = LAY(i, :);
T = L(2) + L(3) * (H - L(1));
if L(3) ~= 0
  p = L(4) * (T / L(2))^(-G0 / (L(3) * R_AIR));
else
  p = L(4) * exp(-G0 * (H - L(1)) / (R_AIR * L(2)));
end
rho = p / (R_AIR * T);
a = sqrt(GAMMA_AIR * R_AIR * T);
end
