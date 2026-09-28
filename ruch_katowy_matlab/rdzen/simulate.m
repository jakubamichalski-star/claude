function R = simulate(A, cfgIn, onProgress)
%SIMULATE Całkowanie równań ruchu 6DOF: RK4 ze stałym krokiem, dokładna lokalizacja zdarzeń
%   (bisekcja długości kroku), restart po zdarzeniu.
%   R = simulate(A, cfg, onProgress)
%   R.n – liczba próbek, R.ch – przebiegi (pola jak w simChannels, wektory kolumnowe),
%   R.meta – dane numeryczne i zdarzenia, R.P – parametry (makeParams).
%   onProgress(t, h) – opcjonalnie, wywoływane co 4000 kroków.
if nargin < 2, cfgIn = struct(); end
if nargin < 3, onProgress = []; end
DEG = pi / 180;
P = makeParams(A, cfgIn);
cfg = P.cfg; num = cfg.num;
if isempty(P.met), hOff = 0; else, hOff = P.hM - P.met.hMDP; end
wind = Wiatr(cfg.wind, P.met, hOff);
wind.init(0);
N = 15;
y = initialState(P);
% krok z najszybszej składowej przy wylocie: co najmniej N kroków na okres max(|s_F|, |s_S|, ω_g)
[~, ~, rho0, a0] = atmAt(P, P.hM);
L0 = linPoint(P, P.u0 / a0, P.u0, rho0, P.p0, 0);
wmax = max([abs(L0.sF0), abs(L0.sS0), abs(L0.wg), 1e-3]);
dt = min(max(2 * pi / (num.stepsPerPeriod * wmax), 1e-6), 5e-3);
kOut = max(1, round(num.stepsPerPeriod / num.samplesPerPeriod));
meta = struct('dt', dt, 'kOut', kOut, 'wmax', wmax, 'steps', 0, 'maxRollCorr', 0, 'clampSteps', 0, ...
  'clampFirst', NaN, 'clampLast', NaN, 'Mmin', Inf, 'Mmax', -Inf, 'apex', [], 'impact', [], ...
  'stopReason', '', 'p0', P.p0, 'QE', P.QE, 'dXcg', P.dXcg, 'kp', P.kp, 'gustPulses', 0, 'cpu_ms', 0);
meta.events = struct('t', {}, 'type', {}, 'text', {});
CH = simChannels(); nch = numel(CH);
iEv = find(strcmp(CH, 'ev')); iRe = find(strcmp(CH, 'reinit')); iM = find(strcmp(CH, 'M')); iCl = find(strcmp(CH, 'clamp'));
ctxClamp = 0;
buf = zeros(4096, nch); nrec = 0;
aStop = num.alphaStop_deg * DEG; thStop = num.thetaStop_deg * DEG;
t = 0; stop = '';
hStart = 0;
reachedAbove = false;
record(t, y, true, true);
t0cpu = tic;
while isempty(stop)
  % zdarzenia aktywne w bieżącym stanie
  evs = struct('name', {}, 'dir', {}, 'level', {}, 'w', {});
  if ~P.frozen
    if ~wind.st.afterApex, evs(end + 1) = struct('name', 'apex', 'dir', -1, 'level', 0, 'w', []); end
    if reachedAbove || t > 0.5, evs(end + 1) = struct('name', 'ground', 'dir', -1, 'level', hStart, 'w', []); end
  end
  evs(end + 1) = struct('name', 'alpha', 'dir', 1, 'level', aStop, 'w', []);
  evs(end + 1) = struct('name', 'theta', 'dir', 1, 'level', thStop, 'w', []);
  wev = wind.events();
  for k = 1:numel(wev), evs(end + 1) = struct('name', 'wind', 'dir', wev(k).dir, 'level', wev(k).level, 'w', wev(k)); end %#ok<AGROW>
  h = dt; hitTB = false;
  tb = min(wind.nextTimeBreak(t), num.tMax);
  if t + h >= tb - 1e-12, h = max(tb - t, 0); hitTB = true; end
  if h > 0, y1 = rk4(t, y, h); else, y1 = y; end
  best = []; bestH = Inf; bestUp = false;
  if h > 0
    for k = 1:numel(evs)
      ev = evs(k);
      g0 = gval(ev, t, y); g1 = gval(ev, t + h, y1);
      if ev.dir > 0, crossed = g0 < 0 && g1 >= 0;
      elseif ev.dir < 0, crossed = g0 > 0 && g1 <= 0;
      else, crossed = (g0 < 0 && g1 >= 0) || (g0 > 0 && g1 <= 0);
      end
      if ~crossed, continue; end
      lo = 0; hi = h;
      for it = 1:64
        if hi - lo <= 1e-11, break; end
        mid = 0.5 * (lo + hi); gm = gval(ev, t + mid, rk4(t, y, mid));
        if ev.dir > 0, cr = gm >= 0;
        elseif ev.dir < 0, cr = gm <= 0;
        elseif g0 < 0, cr = gm >= 0;
        else, cr = gm <= 0;
        end
        if cr, hi = mid; else, lo = mid; end
      end
      if isempty(best) || hi < bestH, best = ev; bestH = hi; bestUp = g0 < 0; end
    end
  end
  if ~isempty(best)
    y = rk4(t, y, bestH); t = t + bestH; meta.steps = meta.steps + 1;
    [y, rc] = projectZeroRoll(y); meta.maxRollCorr = max(meta.maxRollCorr, rc);
    switch best.name
      case 'apex'
        wind.st.afterApex = false;
        meta.apex = struct('t', t, 'h', y(2), 'X1', y(1), 'X3', y(3));
        addEvent(t, 'apex', sprintf('wierzchołek toru: h = %.0f m', y(2)));
        record(t, y, true, false);
        ra = wind.onApex(t);
        if ischar(ra) || ra
          if ischar(ra), txt = 'początek podmuchu na wierzchołku toru'; else, txt = 'zmiana podmuchu na wierzchołku (wybrana gałąź)'; end
          addEvent(t, 'gust', txt); record(t, y, true, ~wind.inTransition());
        end
      case {'ground', 'alpha', 'theta'}
        stop = best.name; record(t, y, true, false);
      case 'wind'
        record(t, y, true, false);
        txt = wind.fire(best.w, t, bestUp);
        addEvent(t, 'gust', txt);
        record(t, y, true, ~wind.inTransition());
    end
  else
    y = y1; t = t + h; meta.steps = meta.steps + 1;
    [y, rc] = projectZeroRoll(y); meta.maxRollCorr = max(meta.maxRollCorr, rc);
    if hitTB
      if t >= num.tMax - 1e-12
        stop = 'tmax'; record(t, y, true, false);
      else
        record(t, y, true, false);
        pc0 = wind.st.piece;
        wind.onTimeBreak(t);
        pc1 = wind.st.piece;
        if pc0 == 0, txt = 'początek podmuchu (chwila zadana w teście)';
        elseif pc1 == 2, txt = 'pełny podmuch (koniec narastania)';
        elseif pc1 == 3, txt = 'początek zanikania podmuchu';
        else, txt = 'koniec podmuchu';
        end
        addEvent(t, 'gust', txt);
        record(t, y, true, ~wind.inTransition());
      end
    elseif mod(meta.steps, kOut) == 0
      record(t, y, false, false);
    end
  end
  if y(2) > 0.5, reachedAbove = true; end
  if ctxClamp
    meta.clampSteps = meta.clampSteps + 1;
    if isnan(meta.clampFirst), meta.clampFirst = t; end
    meta.clampLast = t;
  end
  if ~P.frozen && isempty(stop) && y(2) < -1 && t > 0.2 && ~reachedAbove, stop = 'below'; record(t, y, true, false); end
  if ~isfinite(y(1)) || ~isfinite(y(4)), stop = 'nan'; end
  if ~isempty(onProgress) && mod(meta.steps, 4000) == 0, onProgress(t, y(2)); end
end
meta.cpu_ms = 1000 * toc(t0cpu);
meta.stopReason = stop; meta.tEnd = t;
if strcmp(stop, 'ground'), meta.impact = struct('t', t, 'X1', y(1), 'X3', y(3)); end
meta.gustPulses = wind.st.pulses;
meta.gust = struct('on', wind.on, 'mode', wind.mode, 'h1', wind.h1, 'branch', wind.branch, 'dW', wind.dW, ...
  'shape', wind.shape, 'tr', wind.tr, 'td', wind.td, 'dh', wind.dh, 'HL', wind.HL);
buf = buf(1:nrec, :);
ch = struct();
for k = 1:nch, ch.(CH{k}) = buf(:, k); end
R = struct('n', nrec, 'ch', ch, 'meta', meta, 'P', P);

  % ---------------------------------------------------------------
  function yo = rk4(tt, yy, hh)
    k1 = rhs(P, wind, tt, yy);
    k2 = rhs(P, wind, tt + 0.5 * hh, yy + 0.5 * hh * k1);
    k3 = rhs(P, wind, tt + 0.5 * hh, yy + 0.5 * hh * k2);
    [k4, ctxClamp] = rhs(P, wind, tt + hh, yy + hh * k3);   % ctxClamp: ℳ poza tabelami w ostatnim wywołaniu
    yo = yy + hh / 6 * (k1 + 2 * k2 + 2 * k3 + k4);
  end

  function g = gval(ev, tt, yy)
    switch ev.name
      case 'apex'                          % U₂ = 0
        Tq = quatToT(yy(10:13)); g = -(Tq(3) * yy(4) + Tq(6) * yy(5) + Tq(9) * yy(6));
      case 'ground', g = yy(2) - ev.level;
      case 'alpha'                         % α_tot − α_stop
        Tq = quatToT(yy(10:13)); Wf = wind.W(tt, yy(2)); WIq = [Wf(1), Wf(3), -Wf(2)];
        Vx = yy(4) - (Tq(1) * WIq(1) + Tq(2) * WIq(2) + Tq(3) * WIq(3));
        Vy = yy(5) - (Tq(4) * WIq(1) + Tq(5) * WIq(2) + Tq(6) * WIq(3));
        Vz = yy(6) - (Tq(7) * WIq(1) + Tq(8) * WIq(2) + Tq(9) * WIq(3));
        g = atan2(sqrt(Vy * Vy + Vz * Vz), Vx) - ev.level;
      case 'theta'                         % |θ| − θ_stop
        Tq = quatToT(yy(10:13)); g = abs(asin(max(-1, min(1, -Tq(3))))) - ev.level;
      otherwise                            % przecięcie wysokości (podmuch)
        g = yy(2) - ev.level;
    end
  end

  function record(tt, yy, ev, reinit)
    o = diagnose(P, wind, tt, yy); o(iEv) = double(ev); o(iRe) = double(reinit); ctxClamp = o(iCl);
    nrec = nrec + 1;
    if nrec > size(buf, 1), buf = [buf; zeros(size(buf, 1), nch)]; end
    buf(nrec, :) = o;
    meta.Mmin = min(meta.Mmin, o(iM)); meta.Mmax = max(meta.Mmax, o(iM));
  end

  function addEvent(tt, type, text)
    meta.events(end + 1) = struct('t', tt, 'type', type, 'text', text);
  end
end
