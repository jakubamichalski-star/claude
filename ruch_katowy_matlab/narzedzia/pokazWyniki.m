function tab = pokazWyniki(R)
%POKAZWYNIKI Wyniki obliczenia toru: tor, wielkości u wylotu, ekstrema, komunikaty
%   i przebieg w wybranych chwilach (ok. 24 chwile + próbki po zdarzeniach).
%   tab = pokazWyniki(R) – zwraca tabelę przebiegu w wybranych chwilach (table).
DEG = pi / 180;
S = R.summary; ch = R.ch; po = R.post; n = R.n; num = R.cfg.num;
switch S.stop
  case 'ground', st = 'poziom wylotu';
  case 'alpha', st = sprintf('α_tot > %g° (poza zakresem modelu)', num.alphaStop_deg);
  case 'theta', st = sprintf('|θ| > %g°', num.thetaStop_deg);
  case 'tmax', st = 't_max';
  case 'below', st = 'tor poniżej wylotu';
  case 'nan', st = 'błąd numeryczny';
  otherwise, st = S.stop;
end
I = S.impact; Ap = S.apex;
if S.X3 >= 0, strona = 'w prawo'; else, strona = 'w lewo'; end
fprintf('\nTOR\n');
kv('koniec obliczenia', st); kv('czas lotu', [fmtPL(S.tEnd, 2) ' s']); kv('donośność X1', [fmtPL(S.X1, 0) ' m']);
kv('derywacja X3', [fmtPL(S.X3, 1) ' m ' strona]);
if ~isempty(Ap), kv('wierzchołek', sprintf('h = %s m, t = %s s', fmtPL(Ap.h, 0), fmtPL(Ap.t, 2))); else, kv('wierzchołek', '–'); end
if ~isempty(I)
  kv('prędkość końcowa', sprintf('%s m/s (ℳ %s)', fmtPL(I.V, 1), fmtPL(I.M, 3))); kv('kąt upadku', [fmtPL(I.angle_deg, 2) '°']);
else
  kv('prędkość końcowa', [fmtPL(ch.V(n), 1) ' m/s']); kv('kąt upadku', '–');
end
kv('p końcowe', [fmtPL(ch.p(n), 1) ' rad/s']);
Mz = S.muzzle; cond = Mz.Sg > 1 && 1 / Mz.Sg < Mz.Sd0 * (2 - Mz.Sd0);
if cond, cs = 'spełniony'; else, cs = 'NIESPEŁNIONY'; end
fprintf('\nU WYLOTU (α_tot → 0)\n');
kv('ℳ0', fmtPL(Mz.M, 3)); kv('p0', [fmtPL(Mz.p, 1) ' rad/s']);
kv('S_g', sprintf('%s ± %s', fmtPL(Mz.Sg, 3), fmtPL(Mz.SgSig, 3))); kv('S_d', sprintf('%s ± %s', fmtPL(Mz.Sd0, 3), fmtPL(Mz.SdSig, 3)));
kv('1/S_g < S_d(2 − S_d)', cs);
kv('ν_F, ν_S', sprintf('%s, %s Hz', fmtPL(Mz.nuF, 2), fmtPL(Mz.nuS, 2))); kv('σ_F, σ_S', sprintf('%s, %s 1/s', fmtPL(Mz.sigF, 3), fmtPL(Mz.sigS, 3)));
fprintf('\nEKSTREMA I NUMERYKA\n');
kv('maks. α_tot', sprintf('%s° przy %s s', fmtPL(S.maxAlpha.v, 3), fmtPL(S.maxAlpha.t, 2)));
kv('maks. |ᾱ_e| (STANAG)', sprintf('%s° przy %s s', fmtPL(S.maxAe.v, 3), fmtPL(S.maxAe.t, 1)));
kv('min. S_g', sprintf('%s przy %s s', fmtPL(S.minSg.v, 3), fmtPL(S.minSg.t, 1)));
kv('krok dt', [fmtPL(S.dt * 1000, 4) ' ms']); kv('kroki, próbki', sprintf('%d, %d', S.steps, S.n));
kv('czas: 6DOF + teoria', sprintf('%s + %s s', fmtPL(S.cpu_ms / 1000, 2), fmtPL(S.post_ms / 1000, 2)));
fprintf('\nSTATECZNOŚĆ WZDŁUŻ TORU\n');
pr('  przy bieżącej amplitudzie (funkcja opisująca):', S.runsA);
pr('  teoria liniowa, α_tot → 0:', S.runs0);
fprintf('\nKOMUNIKATY\n');
for k = 1:numel(S.messages)
  switch S.messages(k).level, case 'stop', p = '[STOP] '; case 'warn', p = '[!] '; otherwise, p = ''; end
  fprintf('  - %s%s\n', p, S.messages(k).text);
end
fprintf('  - %s\n', S.stepInfo);
% przebieg w wybranych chwilach
span = ch.t(n) - ch.t(1); step = niceStep(span / 24);
idx = []; i = 1;
for target = 0:step:ch.t(n) + 1e-9
  while i < n && ch.t(i + 1) <= target, i = i + 1; end
  if i < n && abs(ch.t(i + 1) - target) < abs(ch.t(i) - target), idx(end + 1) = i + 1; else, idx(end + 1) = i; end %#ok<AGROW>
end
for k = 1:n
  if ch.ev(k) && (k == n || ch.t(k + 1) - ch.t(k) > 1e-12), idx(end + 1) = k; end %#ok<AGROW>
end
idx = unique(idx);
evTxt = repmat({''}, numel(idx), 1); stan = cell(numel(idx), 1);
for q = 1:numel(idx)
  k = idx(q);
  stan{q} = lower(strtok(strtok(statusText(statusClass(po.status(k))), '–'), '('));
  stan{q} = strtrim(stan{q});
  if ch.ev(k)
    e = R.meta.events(abs([R.meta.events.t] - ch.t(k)) < 1e-9);
    if ~isempty(e), evTxt{q} = strjoin({e.text}, '; ');
    elseif k == 1, evTxt{q} = 'wylot';
    elseif k == n, evTxt{q} = 'koniec obliczenia';
    end
  end
end
if nargout > 0
c = @(v) v(idx);
tab = table(c(ch.t), c(ch.X2), c(ch.V), c(ch.M), c(ch.p), c(ch.alphaTot) / DEG, c(ch.phiA) / DEG, c(po.KFm) / DEG, c(po.KSm) / DEG, ...
  c(po.aeAbs) / DEG, c(po.Sg), c(po.Sd0), c(po.SdFA), c(po.SdSA), c(po.nuF0), c(po.nuS0), c(po.sigF0), c(po.sigS0), ...
  c(po.sFreA), c(po.sSreA), stan, evTxt, 'VariableNames', {'t_s', 'h_m', 'V_ms', 'Mach', 'p_rads', 'alphaTot_deg', 'phiA_deg', ...
  'KF_deg', 'KS_deg', 'alphaE_deg', 'Sg', 'Sd0', 'SdF_A', 'SdS_A', 'nuF_Hz', 'nuS_Hz', 'sigF0', 'sigS0', 'sigF_A', 'sigS_A', 'stan', 'zdarzenie'});
else
  fprintf('\nPRZEBIEG W WYBRANYCH CHWILACH\n');
  fprintf('%8s %6s %7s %6s %7s %7s %6s %7s %7s %7s %8s %6s %6s %6s %6s %6s %7s %7s %7s %7s  %s\n', 't [s]', 'h [m]', 'V [m/s]', 'ℳ', 'p', 'α_tot°', 'φ_A°', ...
    '|K_F|°', '|K_S|°', '|ᾱ_e|°', 'S_g', 'S_d0', 'S_dF', 'S_dS', 'ν_F', 'ν_S', 'σ_F0', 'σ_S0', 'σ_FA', 'σ_SA', 'stan / zdarzenie');
  for q = 1:numel(idx)
    k = idx(q);
    fprintf('%8.3f %6.0f %7.1f %6.3f %7.1f %7.3f %6.1f %7.3f %7.3f %7.3f %8.3f %6.3f %6.3f %6.3f %6.2f %6.3f %7.3f %7.3f %7.3f %7.3f  %s%s\n', ...
      ch.t(k), ch.X2(k), ch.V(k), ch.M(k), ch.p(k), ch.alphaTot(k) / DEG, ch.phiA(k) / DEG, po.KFm(k) / DEG, po.KSm(k) / DEG, po.aeAbs(k) / DEG, ...
      po.Sg(k), po.Sd0(k), po.SdFA(k), po.SdSA(k), po.nuF0(k), po.nuS0(k), po.sigF0(k), po.sigS0(k), po.sFreA(k), po.sSreA(k), stan{q}, ...
      ifelse(isempty(evTxt{q}), '', [' – ' evTxt{q}]));
  end
end
end

function kv(k, v)
fprintf('  %-26s %s\n', k, v);
end

function pr(lbl, rs)
fprintf('%s\n', lbl);
for r = rs
  fprintf('    %8.2f – %8.2f s  %s\n', r.t0, r.t1, statusText(r.key));
end
end

function s = niceStep(x)
e = 10^floor(log10(x)); m = x / e;
if m <= 1, s = e; elseif m <= 2, s = 2 * e; elseif m <= 5, s = 5 * e; else, s = 10 * e; end
end

function v = ifelse(c, a, b)
if c, v = a; else, v = b; end
end
