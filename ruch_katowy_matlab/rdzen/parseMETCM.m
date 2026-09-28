function res = parseMETCM(text, hMDPfallback)
%PARSEMETCM Odczyt komunikatu METCM.
%   res = parseMETCM(text, hMDPfallback)
%   Nagłówek „METCMQ LaLaLaLoLoLo YYGoGoGoG hhhPdPdPd” (hhh – wysokość MDP w dekametrach)
%   i linie „ZZdddFFF TTTTPPPP”: ddd – kierunek, skąd wieje, w dziesiątkach tysięcznych (6400),
%   FFF – prędkość w węzłach, TTTT – średnia temperatura wirtualna strefy w 0,1 K,
%   PPPP – ciśnienie w środku strefy w hPa. Odchyłki liczone względem ISO 2533 na wysokości
%   środka strefy: ΔT = T_v − T_wz, Δρ = p/(R·T_v)/ρ_wz − 1.
%   Wynik: hMDP, hMDPsource, pMDP, zones (tablica struktur), warnings (cellstr).
R_AIR = 287.05287;
s = upper(char(text)); warn = {};
hMDP = []; pMDP = []; src = 'nagłówek komunikatu';
hd = regexp(s, 'METCM\s*(\d)\s*(\d{6})\s+(\d{6})\s+(\d{6})', 'tokens', 'once');
if ~isempty(hd)
  hMDP = str2double(hd{4}(1:3)) * 10;
  pd = str2double(hd{4}(4:6));
  if pd < 500, pMDP = 1000 + pd; else, pMDP = pd; end
else
  if nargin > 1 && isfinite(hMDPfallback), hMDP = hMDPfallback; else, hMDP = 0; end
  src = 'pole formularza (brak nagłówka)';
  warn{end + 1} = 'Nie rozpoznano nagłówka METCM – wysokość MDP z formularza.';
end
Z = metcmZones();
tk = regexp(s, '(?:^|\s)(\d{2})(\d{3})(\d{3})\s*(\d{4})(\d{4})(?=\s|$)', 'tokens');
found = cell(1, numel(Z.mid));
for k = 1:numel(tk)
  m = tk{k};
  zz = str2double(m{1}); ddd = str2double(m{2}); fff = str2double(m{3});
  tttt = str2double(m{4}); pppp = str2double(m{5});
  if zz > 31, warn{end + 1} = sprintf('Linia strefy %s – numer strefy poza 00–31, pominięta.', m{1}); continue; end %#ok<*AGROW>
  if ddd > 640, warn{end + 1} = sprintf('Strefa %s: kierunek %d > 640 dziesiątek tysięcznych.', m{1}, ddd); end
  if tttt < 1500 || tttt > 3500, warn{end + 1} = sprintf('Strefa %s: temperatura %g K poza zakresem 150–350 K.', m{1}, tttt / 10); end
  if pppp < 5 || pppp > 1100, warn{end + 1} = sprintf('Strefa %s: ciśnienie %d hPa poza zakresem.', m{1}, pppp); end
  hMid = Z.mid(zz + 1);
  [Tst, ~, rhoSt] = atmosphere(hMDP + hMid);
  Tv = tttt / 10; p = pppp * 100; rho = p / (R_AIR * Tv);
  found{zz + 1} = struct('n', zz, 'dir_mil', ddd * 10, 'spd_ms', fff * 0.514444, 'spd_kn', fff, 'Tv', Tv, ...
    'p_hPa', pppp, 'rho', rho, 'dT', Tv - Tst, 'dRho', 100 * (rho / rhoSt - 1), 'hMid', hMid);
end
found = found(~cellfun(@isempty, found));
if isempty(found)
  zones = struct('n', {}, 'dir_mil', {}, 'spd_ms', {}, 'spd_kn', {}, 'Tv', {}, 'p_hPa', {}, 'rho', {}, 'dT', {}, 'dRho', {}, 'hMid', {});
  warn{end + 1} = 'Nie znaleziono linii danych „ZZdddFFF TTTTPPPP”.';
else
  zones = [found{:}];
end
res = struct('hMDP', hMDP, 'hMDPsource', src, 'pMDP', pMDP, 'zones', zones);
res.warnings = warn;
end
