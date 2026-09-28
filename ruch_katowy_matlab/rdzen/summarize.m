function S = summarize(R)
%SUMMARIZE Podsumowanie przebiegu: tor, ekstrema, wielkości u wylotu, przedziały stanów i komunikaty.
%   S = summarize(R) – R z simulate z polem R.post (postprocess).
%   S.messages – tablica struktur (level: 'info' | 'warn' | 'stop', text).
DEG = pi / 180; ST = statusBits();
ch = R.ch; n = R.n; meta = R.meta; P = R.P; po = R.post; num = P.cfg.num;
S = struct();
S.messages = struct('level', {}, 'text', {});
S.tEnd = ch.t(n); S.X1 = ch.X1(n); S.X2 = ch.X2(n); S.X3 = ch.X3(n);
S.stop = meta.stopReason; S.apex = meta.apex; S.impact = [];
if strcmp(meta.stopReason, 'ground')
  S.impact = struct('t', ch.t(n), 'X1', ch.X1(n), 'X3', ch.X3(n), 'V', ch.V(n), 'u', ch.u(n), 'M', ch.M(n), ...
    'angle_deg', atan2(-ch.U2(n), hypot(ch.U1(n), ch.U3(n))) / DEG, 'p', ch.p(n));
end
[~, iA] = max(ch.alphaTot); [~, iE] = max(po.aeAbs);
iSg = 1;
for i = 1:n
  if po.Sg(i) > 0 && (po.Sg(i) < po.Sg(iSg) || ~(po.Sg(iSg) > 0)), iSg = i; end
end
S.maxAlpha = struct('v', ch.alphaTot(iA) / DEG, 't', ch.t(iA));
S.maxAe = struct('v', po.aeAbs(iE) / DEG, 't', ch.t(iE));
S.minSg = struct('v', po.Sg(iSg), 't', ch.t(iSg));
S.muzzle = struct('M', ch.M(1), 'V', ch.V(1), 'p', ch.p(1), 'rho', ch.rho(1), 'Sg', po.Sg(1), 'Sd0', po.Sd0(1), ...
  'SdFA', po.SdFA(1), 'SdSA', po.SdSA(1), 'nuF', po.nuF0(1), 'nuS', po.nuS0(1), 'sigF', po.sigF0(1), 'sigS', po.sigS0(1), ...
  'sF', complex(po.sFre0(1), po.sFim0(1)), 'sS', complex(po.sSre0(1), po.sSim0(1)), 'SgSig', po.SgSig(1), 'SdSig', po.SdSig(1));
kA = arrayfun(@(s) statusClass(s), po.status, 'UniformOutput', false);
k0 = arrayfun(@(s) statusClass(s, true), po.status, 'UniformOutput', false);
S.runsA = runs(ch.t, kA);
S.runs0 = runs(ch.t, k0);
switch S.stop
  case 'ground', txt = sprintf('Lot zakończony na poziomie wylotu po %s s.', fmt(S.tEnd, 2));
  case 'alpha', txt = sprintf('POZA ZAKRESEM MODELU: α_tot przekroczył %g° przy t = %s s (koziołkowanie). Symulację zatrzymano.', num.alphaStop_deg, fmt(S.tEnd, 3));
  case 'theta', txt = sprintf('|θ| przekroczył %g° przy t = %s s (osobliwość układu zero-roll). Symulację zatrzymano.', num.thetaStop_deg, fmt(S.tEnd, 3));
  case 'tmax', txt = sprintf('Osiągnięto t_max = %s s.', fmt(num.tMax, 1));
  case 'below', txt = 'Tor od początku opada poniżej poziomu wylotu – sprawdź QE.';
  case 'nan', txt = 'Błąd numeryczny (wartość nieskończona) – symulację przerwano.';
  otherwise, txt = S.stop;
end
if any(strcmp(S.stop, {'ground', 'tmax'})), msg('info', txt); else, msg('stop', txt); end
for e = meta.events
  if strcmp(e.type, 'gust') && isempty(strfind(e.text, 'nieaktywny')), msg('info', sprintf('t = %s s: %s', fmt(e.t, 3), e.text)); end %#ok<STREMP>
end
% podmuch: czas narastania w okresach nutacji i pobudzenie ramienia szybkiego względem skoku (teoria liniowa)
G0 = meta.gust;
if ~isempty(G0) && G0.on
  starts = meta.events([]);
  for e = meta.events
    if strcmp(e.type, 'gust') && (~isempty(regexp(e.text, '^początek podmuchu', 'once')) || ...
        (strcmp(G0.mode, 'height') && ~isempty(regexp(e.text, '^granica warstwy', 'once')) && isempty(strfind(e.text, 'nieaktywny')))) %#ok<STREMP>
      starts(end + 1) = e; %#ok<AGROW>
    end
  end
  if strcmp(G0.mode, 'height'), starts = starts(1:min(1, end)); end
  for e = starts
    i = 1; while i < n && ch.t(i + 1) <= e.t + 1e-12, i = i + 1; end
    nuF = po.nuF0(i); sF = complex(po.sFreA(i), po.sFimA(i));
    if ~(nuF > 0), continue; end
    if strcmp(G0.shape, 'step'), tr = 0;
    elseif strcmp(G0.mode, 'time'), tr = G0.tr;
    else, tr = G0.dh / max(1e-6, abs(ch.U2(i)));
    end
    if ~(tr > 0)
      msg('info', sprintf('Podmuch w t = %s s jest skokowy: ramię szybkie (nutacja, ν_F = %s Hz) pobudzone w pełni.', fmt(e.t, 3), fmt(nuF, 2)));
      continue;
    end
    if strcmp(G0.mode, 'height'), zast = ' (zastępcze δh/|U₂|)'; else, zast = ''; end
    msg('info', sprintf(['Podmuch w t = %s s: narastanie τ_r = %s ms%s = %s okresu nutacji (ν_F = %s Hz); ' ...
      'pobudzenie ramienia szybkiego względem podmuchu skokowego wg teorii liniowej: %s.'], ...
      fmt(e.t, 3), fmt(tr * 1000, 1), zast, fmt(tr * nuF, 2), fmt(nuF, 2), fmt(gustExcitation(G0.shape, tr, sF), 3)));
  end
end
% przedziały stanów innych niż „stateczny”
pretty = struct('gyro', 'niestateczny giroskopowo (S_g < 1)', 'dynF', 'niestateczny dynamicznie – składowa szybka', ...
  'dynS', 'niestateczny dynamicznie – składowa wolna', 'dynFS', 'niestateczny dynamicznie – obie składowe', ...
  'limit', 'cykl graniczny', 'weak', sprintf('słaba nadążność (|ᾱ_e| > %g°)', num.alphaEWarn_deg));
lst = {S.runs0, 'Teoria liniowa (α_tot → 0)'; S.runsA, 'Przy bieżącej amplitudzie (funkcja opisująca)'};
for q = 1:2
  rs = lst{q, 1}; lbl = lst{q, 2};
  bad = rs(~strcmp({rs.key}, 'ok') & ([rs.t1] - [rs.t0]) >= 0.1);
  for r = bad(1:min(8, end))
    msg('warn', sprintf('%s: %s dla t = %s–%s s.', lbl, pretty.(r.key), fmt(r.t0, 2), fmt(r.t1, 2)));
  end
  if numel(bad) > 8, msg('warn', sprintf('%s: … i %d dalszych przedziałów (zob. pasek stanu).', lbl, numel(bad) - 8)); end
end
% przedziały słabej nadążności (|ᾱ_e| > próg) – w pasku stanu przysłaniane przez niestateczność
wk = runs(ch.t, num2cell(bitand(po.status, ST.WEAK) > 0));
wk = wk([wk.key] == 1 & ([wk.t1] - [wk.t0]) >= 0.1);
for r = wk(1:min(6, end))
  msg('warn', sprintf('Słaba nadążność: |ᾱ_e| > %g° dla t = %s–%s s (oś nie nadąża za zakrzywiającym się torem).', num.alphaEWarn_deg, fmt(r.t0, 1), fmt(r.t1, 1)));
end
% przedziały ostrzeżeń α_tot
aw = runs(ch.t, num2cell(bitand(po.status, ST.AWARN) > 0));
aw = aw([aw.key] == 1);
for r = aw
  msg('warn', sprintf('α_tot > %g° dla t = %s–%s s.', num.alphaWarn_deg, fmt(r.t0, 2), fmt(r.t1, 2)));
end
if meta.clampSteps > 0
  msg('warn', sprintf('ℳ poza zakresem tabel (%s–%s) dla t ≈ %s–%s s (ℳ_min = %s, ℳ_max = %s): użyto wartości skrajnych.', ...
    fmt(P.A.mach(1), 2), fmt(P.A.mach(end), 2), fmt(meta.clampFirst, 1), fmt(meta.clampLast, 1), fmt(meta.Mmin, 3), fmt(meta.Mmax, 3)));
end
G = meta.gust;
if ~isempty(G) && G.on
  reached = any(ch.piece ~= 0);
  if ~reached && strcmp(G.branch, 'apex')
    msg('warn', 'Podmuch nie wystąpił: lot zakończył się przed wierzchołkiem toru.');
  elseif ~reached
    if ~isempty(meta.apex), hmax = meta.apex.h; else, hmax = max(ch.X2); end
    if hmax < G.h1
      msg('stop', sprintf('Podmuch nie wystąpił: h₁ = %s m leży POWYŻEJ wierzchołka toru (h_max = %s m).', fmt(G.h1, 0), fmt(hmax, 0)));
    else
      msg('warn', sprintf('Podmuch nie wystąpił: wybrana gałąź toru nie przecięła h₁ = %s m przed końcem lotu.', fmt(G.h1, 0)));
    end
  end
  if strcmp(G.mode, 'height') && ~strcmp(G.branch, 'both') && ~isempty(meta.apex)
    hA = meta.apex.h; top = G.h1 + G.HL;
    if hA > G.h1 && hA < top, msg('warn', 'Wierzchołek toru leży w warstwie podmuchu: przełączenie gałęzi daje skokową zmianę wiatru na wierzchołku.'); end
  end
end
if ~isempty(P.met)
  if ~isempty(meta.apex), hA = meta.apex.h; else, hA = max(ch.X2); end
  hr = P.hM + hA - P.met.hMDP; W = metW(P.met, hr);
  msg('info', sprintf(['Meteo: warstwy METCM, MDP = %s m n.p.m., odchyłki od ISO 2533 interpolowane między środkami warstw (PCHIP). ' ...
    'Na wysokości wierzchołka (%s m nad MDP): ΔT = %s K, Δρ = %s %%, W₁ = %s m/s, W₃ = %s m/s.'], ...
    fmt(P.met.hMDP, 0), fmt(hr, 0), fmt(pchipVal(P.met.fT, hr), 2), fmt(pchipVal(P.met.fR, hr), 2), fmt(W(1), 1), fmt(W(3), 1)));
end
S.dt = meta.dt; S.kOut = meta.kOut; S.steps = meta.steps; S.cpu_ms = meta.cpu_ms; S.n = n;
S.stepInfo = sprintf('dt = %s ms (%g kroków na okres najszybszej składowej u wylotu: %s Hz), zapis co %d kroków, %d kroków, %d próbek.', ...
  fmt(meta.dt * 1000, 4), num.stepsPerPeriod, fmt(meta.wmax / (2 * pi), 2), meta.kOut, meta.steps, n);

  function msg(level, text)
    S.messages(end + 1) = struct('level', level, 'text', text);
  end
end

% ---------------------------------------------------------------------
function out = runs(t, key)
% przedziały stałej wartości klucza: key – cell (tekst albo liczba/logiczna) o długości t
n = numel(t); out = struct('key', {}, 't0', {}, 't1', {});
k0 = 1;
for i = 2:n + 1
  if i == n + 1 || ~isequal(key{i}, key{k0})
    out(end + 1) = struct('key', key{k0}, 't0', t(k0), 't1', t(i - 1)); %#ok<AGROW>
    k0 = i;
  end
end
if numel(out) > 1, out = out([out.t1] > [out.t0]); end
end

function s = fmt(x, d)
% liczba z przecinkiem dziesiętnym (zapis polski); '–' dla wartości nieskończonych
if ~isfinite(x), s = '–'; return; end
s = strrep(sprintf('%.*f', d, x), '.', ',');
end
