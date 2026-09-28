function met = metProfile(M, azDeg, hDef)
%METPROFILE Profil meteo z tabeli stref METCM.
%   met = metProfile(M, azDeg, hDef)
%   Wiatr (kierunek, skąd wieje – od północy geograficznej zgodnie ze wskazówkami zegara;
%   prędkość [m/s]) przeliczony na składowe w układzie ogniowym dla azymutu strzału azDeg,
%   odchyłki temperatury ΔT [K] i gęstości Δρ [%] od atmosfery ISO 2533.
%   Wartość strefy przypisana jest środkowi strefy (linia 00 – wysokości 0); między środkami
%   interpolacja PCHIP, poza zakresem wartości skrajne. [] – gdy warstwy wyłączone.
%   Pola: hMDP, xs, w1, w3, dT, dR, f1, f3, fT, fR (struktury pchipMake).
%   Wiatr na wysokości hr nad MDP: metW(met, hr), gradient: metdW(met, hr).
met = [];
if isempty(M) || ~isfield(M, 'enabled') || ~M.enabled || ~isfield(M, 'zones') || isempty(M.zones), return; end
DEG = pi / 180;
if isfield(M, 'hMDP') && ~isempty(M.hMDP) && isfinite(M.hMDP), hMDP = M.hMDP; else, hMDP = hDef; end
Z = metcmZones();
az = azDeg * DEG; n = min(numel(M.zones), numel(Z.mid));
if strcmp(M.dirUnit, 'deg'), k0 = DEG; else, k0 = 2 * pi / 6400; end
xs = zeros(1, n); w1 = xs; w3 = xs; dT = xs; dR = xs;
for k = 1:n
  z = M.zones(k); D = z.dir * k0; s = z.spd;
  xs(k) = Z.mid(k);
  w1(k) = -s * cos(D - az); w3(k) = -s * sin(D - az);     % wiatr „skąd” D → wektor ku D + π
  dT(k) = z.dT; dR(k) = z.dRho;
end
met = struct('hMDP', hMDP, 'xs', xs, 'w1', w1, 'w3', w3, 'dT', dT, 'dR', dR, ...
  'f1', pchipMake(xs, w1), 'f3', pchipMake(xs, w3), 'fT', pchipMake(xs, dT), 'fR', pchipMake(xs, dR));
end
