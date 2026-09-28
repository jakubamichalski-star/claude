function met = metcmToConfig(res, met)
%METCMTOCONFIG Przeniesienie odczytanego komunikatu METCM (parseMETCM) do tabeli stref konfiguracji.
%   met = metcmToConfig(res, met) – met: pole cfg.met (dirUnit decyduje o jednostce kierunku).
%   Strefy brakujące w komunikacie przyjmują wartości ostatniej wcześniejszej strefy
%   (przed pierwszą – zera). Włącza warstwy i ustawia wysokość MDP z komunikatu.
Z = metcmZones(); r2 = @(v) round(v * 100) / 100;
zones = repmat(struct('dir', 0, 'spd', 0, 'dT', 0, 'dRho', 0), 1, numel(Z.n));
nz = [res.zones.n]; last = [];
for k = 1:numel(Z.n)
  j = find(nz == Z.n(k), 1);
  if ~isempty(j), last = res.zones(j); end
  if ~isempty(last)
    if strcmp(met.dirUnit, 'deg'), dir = round(last.dir_mil * 3600 / 6400) / 10; else, dir = last.dir_mil; end
    zones(k) = struct('dir', dir, 'spd', r2(last.spd_ms), 'dT', r2(last.dT), 'dRho', r2(last.dRho));
  end
end
met.zones = zones; met.enabled = true; met.hMDP = res.hMDP;
end
