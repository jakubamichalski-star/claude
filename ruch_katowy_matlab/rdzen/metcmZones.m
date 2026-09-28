function Z = metcmZones()
%METCMZONES Strefy komunikatu METCM (STANAG 4082).
%   Linia 00 – powierzchnia MDP; strefy 01–31 (wysokości nad MDP): 0–200, 200–500, 500–1000 m,
%   dalej co 500 m do 5 km, co 1000 m do 20 km i co 2000 m do 30 km.
%   Z.n, Z.lo, Z.hi, Z.mid – wektory 1×32 (element k odpowiada strefie k−1).
b = [0 200 500 1000:500:5000 6000:1000:20000 22000:2000:30000];
Z.n = 0:numel(b) - 1;
Z.lo = [0 b(1:end - 1)];
Z.hi = [0 b(2:end)];
Z.mid = (Z.lo + Z.hi) / 2;
end
