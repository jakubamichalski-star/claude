function out = runTests(A, userCfg, verbose)
%RUNTESTS Testy modelu (sekcja 10 specyfikacji + testy dodatkowe).
%   out = runTests(A, userCfg, verbose)
%   Warunki testów stałe (u₀ = 800 m/s, QE = 30°, t_c = 20 kal., h_M = 0), współczynniki
%   nominalne (f = 1, f_Δ = 0); masa, bezwładność i normowanie – z userCfg (bieżące).
%   out – tablica struktur: id, name, pass, value, details (cellstr), ms.
%   verbose (domyślnie true) – wypisywanie wyników na bieżąco.
if nargin < 2 || isempty(userCfg), userCfg = struct(); end
if nargin < 3, verbose = true; end
DEG = pi / 180; GAMMA_AIR = 1.4; R_AIR = 287.05287;
d0 = defaultConfig();
base = d0;
if isfield(userCfg, 'norm'), base.norm = userCfg.norm; end
if isfield(userCfg, 'mass'), base.mass = deepMerge(d0.mass, userCfg.mass); end
if isfield(userCfg, 'mult'), base.mult = deepMerge(d0.mult, userCfg.mult); end
base.fire = struct('u0', 800, 'QE', 30, 'QEunit', 'deg', 'tc', 20, 'bore_mm', 155, 'hM', 0, 'q0', 1, 'r0', 0, ...
  'alphaTot0_deg', 0, 'phiA0_deg', 0);
base.wind = struct('base', struct('W', [0 0 0], 'table', []), 'gust', struct('enabled', false));
base.earth.coriolis = false;
if isfield(userCfg, 'num') && isfield(userCfg.num, 'stepsPerPeriod')
  base.num.stepsPerPeriod = max(100, userCfg.num.stepsPerPeriod);
end
sim = @(over) simulate(A, deepMerge(base, over));
out = struct('id', {}, 'name', {}, 'pass', {}, 'value', {}, 'details', {}, 'ms', {});
rel = @(a, b) abs(a - b) ./ abs(b);
frozen = @(u0, extra) deepMerge(st('test', st('frozen', true), 'fire', st('u0', u0, 'QE', 0, 'q0', 0.02, 'r0', 0), ...
  'num', st('tMax', 3, 'samplesPerPeriod', 25)), extra);
fits = struct();

% ---- T1 / T1b: ν, σ z 6DOF (Prony) vs sekcja 4 ----
freqTest('T1', 800, 0.02, 3, 'ℳ ≈ 2,35, naddźwiękowo');
[~, ~, ~, a0] = atmosphere(0);
freqTest('T1b', 0.6 * a0, 0.001, 4, 'ℳ = 0,6, poddźwiękowo');

% ---- T2: zerowe momenty ----
T('T2', 'Zerowe momenty, q₀ = r₀ = 0: oś x̄ zachowuje kierunek, α_tot rośnie z zakrzywieniem toru', @testT2);
% ---- T3: p = 0 ----
T('T3', 'p = 0 i C_Mα > 0: α_tot rośnie wykładniczo', @testT3);
% ---- T4: S_g, S_d vs obliczenie ręczne ----
T('T4', 'S_g i S_d: implementacja vs obliczenie ręczne (ℳ = 2,0, h = 0, p = 1500 rad/s)', @testT4);
% ---- pełny tor do T5, T7, T10 ----
Rf = sim(st('num', st('tMax', 25))); Rf.post = postprocess(Rf);
T('T5', 'p > 0: oba ramiona obracają się w stronę rosnącego φ_A, ᾱ_e wskazuje w prawo (φ_A ≈ 90°)', @testT5);
T('T6', 'Zbieżność względem kroku całkowania (dt, dt/2, dt/4; 3 s lotu z zaburzeniem wylotowym)', @testT6);
T('T7', 'Po wytłumieniu oscylacji średni kąt natarcia z 6DOF zgodny z ᾱ_e', @testT7);
% ---- T8, T9: podmuch skokowy i łagodny ----
stepArms = [];
T('T8', 'Podmuch skokowy (dodatkowy): skok ξ = i·ΔW/V oraz |K_F|/|K_S| = |s_S|/|s_F|', @testT8);
T('T9', 'Podmuch łagodny (dodatkowy): amplituda nutacji po narastaniu vs teoria liniowa', @testT9);
T('T10', 'Kwaternion ↔ kąty Eulera (A.8), kinematyka (8) i więzy zero-roll (dodatkowy)', @testT10);
T('T11', 'Atmosfera ISO 2533: wartości tablicowe (dodatkowy)', @testT11);
T('T12', 'Ramiona liniowe przy zmianie prędkości od grawitacji: amplituda ramienia F z 6DOF vs K_F (dodatkowy)', @testT12);
T('T13', 'Funkcja opisująca C_mag-m: przypadki graniczne i ramię F po podmuchu w wierzchołku (dodatkowy)', @testT13);
T('T14', 'Warstwy METCM: tożsamość przy zerowych odchyłkach, atmosfera z odchyłkami, wiatr, odczyt komunikatu (dodatkowy)', @testT14);
if verbose
  fprintf('Wynik: %d z %d testów zaliczonych.\n', sum([out.pass]), numel(out));
end

  % =================================================================
  function T(id, name, fn)
    t0 = tic;
    try
      r = fn();
    catch err
      r = struct('pass', false, 'value', ['błąd wykonania: ' err.message], 'details', {{getReport(err)}});
    end
    r.id = id; r.name = name; r.ms = 1000 * toc(t0);
    out(end + 1) = orderfields(r, out);
    if verbose
      if r.pass, ok = 'ZALICZONY'; else, ok = 'NIEZALICZONY'; end
      fprintf('%-4s %-12s %s\n     %s\n', id, ok, name, r.value);
      for q = 1:numel(r.details), if ~isempty(r.details{q}), fprintf('       %s\n', r.details{q}); end, end
    end
  end

  function L = muzzle(R, u0)
    [~, ~, rho, aS] = atmosphere(R.P.hM);
    L = linPoint(R.P, u0 / aS, u0, rho, R.P.p0, 0);
  end

  function freqTest(id, u0, q0, tMax, label)
    T(id, sprintf('Małe kąty, współczynniki zamrożone (%s): ν_F, ν_S, σ_F, σ_S z 6DOF vs wzory sekcji 4', label), ...
      @() freqRun(id, u0, q0, tMax));
  end

  function r = freqRun(id, u0, q0, tMax)
    R = sim(frozen(u0, st('fire', st('q0', q0), 'num', st('tMax', tMax))));
    L = muzzle(R, u0);
    u = uniformXi(R, 0, tMax);
    f = prony2(u.z, u.dt); fits.(id) = struct('f', f, 'L', L);
    fit = [imag(f.sF) / (2 * pi), imag(f.sS) / (2 * pi), real(f.sF), real(f.sS)];
    cf = [L.nuF0, L.nuS0, L.sigF0, L.sigS0];
    exr = [imag(L.sF0) / (2 * pi), imag(L.sS0) / (2 * pi), real(L.sF0), real(L.sS0)];
    e = rel(fit, cf); eEx = rel(fit, exr); emax = max(e);
    nm = {'ν_F [Hz]', 'ν_S [Hz]', 'σ_F [1/s]', 'σ_S [1/s]'};
    det = {sprintf('ℳ = %s, p = %s rad/s, S_g = %s, S_d = %s; g = 0, C_D = 0, C_spin = 0, ρ = const; dopasowanie Prony''ego: %d próbek, reszta %s', ...
      f4(R.ch.M(1), 3), f4(R.P.p0, 1), f4(L.Sg, 3), f4(L.Sd0, 3), numel(u.z), f4(f.relRes, 7))};
    for k = 1:4
      det{end + 1} = sprintf('%s: 6DOF %s | wzór sekcji 4 %s (%s) | pierwiastki dokładne %s (%s)', nm{k}, f4(fit(k), 4), f4(cf(k), 4), pct(e(k)), f4(exr(k), 4), pct(eEx(k))); %#ok<AGROW>
    end
    r = st('pass', emax < 0.02 && u.uniform, 'value', sprintf('maks. różnica %s (próg 2 %%)', pct(emax)), 'details', det);
  end

  function r = testT2()
    R = sim(st('f', st('CMa', 0, 'CMq', 0, 'Cmagm', 0, 'Cspin', 0), 'fire', st('q0', 0, 'r0', 0), 'num', st('tMax', 10)));
    ch = R.ch;
    dth = max(abs(ch.theta - ch.theta(1))); dps = max(abs(ch.psi - ch.psi(1)));
    k = 2:R.n; dec = ch.t(k) > ch.t(k - 1) & ch.alphaTot(k) < ch.alphaTot(k - 1) - 1e-13;
    mono = ~any(dec);
    dB5 = 0;
    for i = 1:R.n
      ax = zrAxesF(ch.psi(i), ch.theta(i));
      VF = [ch.U1(i) - ch.W1(i), ch.U2(i) - ch.W2(i), ch.U3(i) - ch.W3(i)];
      aGeo = atan2(norm(cross(ax.x, VF)), dot(ax.x, VF));      % kąt x̄–V̄ geometrycznie
      dB5 = max(dB5, abs(aGeo - ch.alphaTot(i)));
    end
    k = R.n; gam = atan2(ch.U2(k), hypot(ch.U1(k), ch.U3(k)));
    r = st('pass', dth < 1e-9 && dps < 1e-9 && mono && dB5 < 1e-9, ...
      'value', sprintf('α_tot(10 s) = %s°, |Δθ|, |Δψ| ≤ %s rad', f4(ch.alphaTot(k) / DEG, 3), ex(max(dth, dps), 1)), ...
      'details', {sprintf('maks. |θ − θ₀| = %s rad, maks. |ψ − ψ₀| = %s rad (próg 1e-9)', ex(dth, 2), ex(dps, 2)), ...
        sprintf('α_tot rośnie monotonicznie: %s', tn(mono)), ...
        sprintf('α_tot z (B.5) (składowe zero-roll) vs kąt między x̄ i V̄ w układzie ogniowym: maks. różnica %s rad', ex(dB5, 2)), ...
        sprintf('po 10 s: QE − γ = %s° (kąt obrotu wektora prędkości), α_tot = %s° (różnica ze składowej bocznej od siły Magnusa)', ...
          f4((R.P.QE - gam) / DEG, 3), f4(ch.alphaTot(k) / DEG, 3))});
  end

  function r = testT3()
    R = sim(st('test', st('frozen', true), 'mult', st('p0', 0), 'fire', st('u0', 800, 'QE', 0, 'q0', 1e-4, 'r0', 0), ...
      'num', st('tMax', 1.5, 'samplesPerPeriod', 25)));
    L = muzzle(R, 800);
    sPlus = (-L.eta + sqrt(L.eta^2 + 4 * L.wa2)) / 2;
    ch = R.ch; a = ch.alphaTot;
    planar = max([0; abs(ch.phiA(a > 1e-6))]);
    sel = a >= 1e-3 & a <= 5e-2 & ~ch.ev; nn = sum(sel);
    slope = linSlope(ch.t(sel), log(a(sel)));
    e = rel(slope, sPlus);
    r = st('pass', L.CMa > 0 && e < 0.02 && planar < 1e-6 && nn > 10, ...
      'value', sprintf('wykładnik %s 1/s vs teoria %s 1/s (%s)', f4(slope, 3), f4(sPlus, 3), pct(e)), ...
      'details', {sprintf('C_Mα = %s > 0, ω_α = %s 1/s, η = %s 1/s; teoria: s = [−η + √(η² + 4ω_α²)]/2', f4(L.CMa, 3), f4(sqrt(L.wa2), 3), f4(L.eta, 4)), ...
        sprintf('regresja ln α_tot na %d próbkach (α_tot od 0,057° do 2,9°)', nn), sprintf('ruch płaski: maks. |φ_A| = %s rad', ex(planar, 2))});
  end

  function r = testT4()
    P = makeParams(A, base);
    [~, i] = min(abs(A.mach - 2.0));
    M = A.mach(i); [~, ~, rho, aS] = atmosphere(0); V = M * aS; p = 1500;
    kp = P.kp; dX = P.dXcg; m = P.m; Ix = P.Ix; IT = P.IT; d = P.d;
    Cnpa0 = A.cnpaData(i, 1);
    Cma = A.Cma(i) + A.CNa(i) * dX; CLa = A.CLa(i); CDo = A.CDo(i); Cmq = kp * A.Cmq(i); Cmm = kp * (Cnpa0 + A.CYpa(i) * dX);
    SgH = 2 * Ix^2 * p^2 / (pi * rho * IT * d^3 * V^2 * Cma);
    SdH = 2 * (CLa + m * d^2 / Ix * Cmm) / (CLa - CDo - m * d^2 / IT * Cmq);
    L = linPoint(P, M, V, rho, p, 0);
    e = max(rel(L.Sg, SgH), rel(L.Sd0, SdH));
    if 1 / SgH < SdH * (2 - SdH), wn = 'stateczny'; else, wn = 'niestateczny'; end
    r = st('pass', e < 1e-9, 'value', sprintf('S_g = %s, S_d = %s; różnica %s', f4(SgH, 4), f4(SdH, 4), ex(e, 1)), ...
      'details', {sprintf('dane: ℳ = %g, ρ = %s kg/m³, V = %s m/s, p = %g rad/s, d = %g m, m = %g kg, I_x = %g, I_T = %g kg·m², ΔX_CG = %s kal., κ_p = %g', ...
          M, f4(rho, 4), f4(V, 3), p, d, m, Ix, IT, f4(dX, 4), kp), ...
        sprintf('z pliku (wiersz ℳ = %g): Cma = %g, CNa = %g, CLa = %g, CDo = %g, Cmq = %g, Cnpa(α = 0) = %g, CYpa = %g', ...
          M, A.Cma(i), A.CNa(i), CLa, CDo, A.Cmq(i), Cnpa0, A.CYpa(i)), ...
        sprintf('S_g = 2·I_x²·p²/(π·ρ·I_T·d³·V²·C_Mα) = %.10g | implementacja: %.10g', SgH, L.Sg), ...
        sprintf('S_d = 2·[C_Lα + (md²/I_x)·C_mag-m]/[C_Lα − C_D0 − (md²/I_T)·C_Mq] = %.10g | implementacja: %.10g', SdH, L.Sd0), ...
        sprintf('warunek 1/S_g < S_d·(2 − S_d): %s < %s → %s', f4(1 / SgH, 4), f4(SdH * (2 - SdH), 4), wn)});
  end

  function r = testT5()
    rot = isfield(fits, 'T1');
    if rot, f1 = fits.T1.f; rot = imag(f1.sF) > 0 && imag(f1.sS) > 0; end
    if rot && isfield(fits, 'T1b'), f1b = fits.T1b.f; rot = imag(f1b.sF) > 0 && imag(f1b.sS) > 0; end
    po = Rf.post;
    ae = meanOver(Rf, 9.5, 10.5, po.aeRe, po.aeIm); xm = meanOver(Rf, 9.5, 10.5, po.xiRe, po.xiIm);
    dirAe = angle(ae) / DEG; dirXi = angle(xm) / DEG;
    drift = Rf.ch.X3(end);
    s2 = Rf.ch.t <= 2;
    dl = max(abs(complex(po.xiLinRe(s2) - po.xiRe(s2), po.xiLinIm(s2) - po.xiIm(s2)))); dn = max(po.devE(s2));
    det = {};
    if isfield(fits, 'T1'), det{end + 1} = sprintf('T1: ν_F = %s Hz, ν_S = %s Hz (dodatnie → ruch w stronę rosnącego φ_A)', f4(imag(fits.T1.f.sF) / (2 * pi), 3), f4(imag(fits.T1.f.sS) / (2 * pi), 3)); end
    if isfield(fits, 'T1b'), det{end + 1} = sprintf('T1b: ν_F = %s Hz, ν_S = %s Hz', f4(imag(fits.T1b.f.sF) / (2 * pi), 3), f4(imag(fits.T1b.f.sS) / (2 * pi), 3)); end
    det = [det, {sprintf('t = 10 s (średnia z ±0,5 s): φ_A(ᾱ_e STANAG) = %s° (próg 90° ± 5°), φ_A(średnie ξ z 6DOF) = %s° (próg 90° ± 10°)', f4(dirAe, 2), f4(dirXi, 2)), ...
      sprintf('derywacja X₃(25 s) = %s m (w prawo)', f4(drift, 2)), ...
      sprintf('nakładka liniowa (informacyjnie): maks. |ξ_lin − ξ| w 0–2 s = %s° przy maks. |ξ − ξ_e| = %s° (%s)', f4(dl / DEG, 5), f4(dn / DEG, 4), pct(dl / dn))}];
    r = st('pass', rot && abs(dirAe - 90) < 5 && abs(dirXi - 90) < 10 && drift > 0, ...
      'value', sprintf('ν_F, ν_S > 0: %s; kierunek ᾱ_e %s°, 6DOF %s°; derywacja %s m', tn(rot), f4(dirAe, 1), f4(dirXi, 1), f4(drift, 1)), ...
      'details', det);
  end

  function r = testT6()
    Nn = [100 200 400]; Rs = cell(1, 3);
    for k = 1:3, Rs{k} = sim(st('num', st('tMax', 3, 'stepsPerPeriod', Nn(k), 'samplesPerPeriod', 20))); end
    xiE = @(R) R.ch.alphaTot(end) * exp(1i * R.ch.phiA(end));
    XE = @(R) [R.ch.X1(end), R.ch.X2(end), R.ch.X3(end)];
    e1 = abs(xiE(Rs{1}) - xiE(Rs{2})); e2 = abs(xiE(Rs{2}) - xiE(Rs{3}));
    x1 = norm(XE(Rs{1}) - XE(Rs{2})); x2 = norm(XE(Rs{2}) - XE(Rs{3}));
    order = log2(e1 / e2);
    amax = max(Rs{3}.ch.alphaTot);
    r = st('pass', order > 3.3 && order < 4.7 && e1 < 1e-4 * amax, ...
      'value', sprintf('rząd zbieżności ≈ %s; |Δξ(dt; dt/2)| = %s°', f4(order, 2), ex(e1 / DEG, 2)), ...
      'details', {sprintf('dt = %s ms; t_końc = %s s', f4(Rs{1}.meta.dt * 1e3, 5), f4(Rs{1}.ch.t(end), 6)), ...
        sprintf('|ξ(dt) − ξ(dt/2)| = %s rad, |ξ(dt/2) − ξ(dt/4)| = %s rad, maks. α_tot = %s°', ex(e1, 3), ex(e2, 3), f4(amax / DEG, 3)), ...
        sprintf('|X̄(dt) − X̄(dt/2)| = %s m, |X̄(dt/2) − X̄(dt/4)| = %s m', ex(x1, 3), ex(x2, 3))});
  end

  function r = testT7()
    po = Rf.post; det = {}; ok = true; worst = 0;
    for tc = [10 20]
      xm = meanOver(Rf, tc - 0.5, tc + 0.5, po.xiRe, po.xiIm);
      xe = meanOver(Rf, tc - 0.5, tc + 0.5, po.xieRe, po.xieIm);
      xe0 = meanOver(Rf, tc - 0.5, tc + 0.5, po.xie0Re, po.xie0Im);
      ae = meanOver(Rf, tc - 0.5, tc + 0.5, po.aeRe, po.aeIm);
      eFull = abs(xm - xe) / abs(xe); eSt = abs(abs(xm) - abs(ae)) / abs(xm);
      eQs = abs(xm - xe0) / abs(xe0);
      dDir = (angle(ae) - angle(xm)) / DEG;
      [~, iM] = min(abs(Rf.ch.t - tc));
      ok = ok && eFull < 0.03 && eSt < 0.03; worst = max([worst, eFull, eSt]);
      det = [det, {sprintf('t = %g s (ℳ = %s): |średnie ξ| = %s°, φ_A = %s°', tc, f4(Rf.ch.M(iM), 3), f4(abs(xm) / DEG, 4), f4(angle(xm) / DEG, 2)), ...
        sprintf('   ξ_e (rozwiązanie szczególne z poprawką na zmienność): %s°, %s° → różnica zespolona %s', f4(abs(xe) / DEG, 4), f4(angle(xe) / DEG, 2), pct(eFull)), ...
        sprintf('   ξ_e0 (quasi-stacjonarne, informacyjnie): %s°, %s° → %s', f4(abs(xe0) / DEG, 4), f4(angle(xe0) / DEG, 2), pct(eQs)), ...
        sprintf('   |ᾱ_e| STANAG: %s° (%s), kierunek różni się o %s°', f4(abs(ae) / DEG, 4), pct(eSt), f4(dDir, 2))}]; %#ok<AGROW>
    end
    r = st('pass', ok, 'value', sprintf('maks. różnica %s (próg 3 %%)', pct(worst)), 'details', det);
  end

  function R = gustRun(shape, tr)
    R = sim(st('test', st('frozen', true), 'fire', st('u0', 800, 'QE', 0, 'q0', 0, 'r0', 0), 'num', st('tMax', 3.5, 'samplesPerPeriod', 25), ...
      'wind', st('gust', st('enabled', true, 'mode', 'time', 'triggerAt', 0.5, 'A', 10, 'azimuth_deg', 90, 'W2', 0, 'shape', shape, 'tr', tr, 'td', 0))));
  end

  function a = armsAfter(R, te)
    u = uniformXi(R, te + 1e-9, 3.5); f = prony2(u.z, u.dt); back = u.ts(1) - te;
    a = struct('f', f, 'KF', f.KF * exp(-f.sF * back), 'KS', f.KS * exp(-f.sS * back), 'uniform', u.uniform);
  end

  function r = testT8()
    R = gustRun('step', 0);
    ch = R.ch; idx = find(ch.ev & abs(ch.t - 0.5) < 1e-9);
    xiAt = @(i) ch.alphaTot(i) * exp(1i * ch.phiA(i));
    dxi = xiAt(idx(end)) - xiAt(idx(1));
    Vb = ch.V(idx(1)); expd = 1i * atan2(10, Vb);
    eJ = abs(dxi - expd) / abs(expd);
    a = armsAfter(R, 0.5); stepArms = a;
    ip = idx(end);
    L = linPoint(R.P, ch.M(ip), ch.V(ip), ch.rho(ip), ch.p(ip), ch.alphaTot(ip));
    ratio = abs(a.KF) / abs(a.KS); th = abs(L.sSA) / abs(L.sFA);
    e = rel(ratio, th);
    r = st('pass', eJ < 0.01 && e < 0.02 && a.uniform, ...
      'value', sprintf('skok ξ: %s; |K_F|/|K_S| = %s vs %s (%s)', pct(eJ), f4(ratio, 4), f4(th, 4), pct(e)), ...
      'details', {sprintf('ΔW₃ = +10 m/s (wiatr ku osi 3) w t = 0,5 s, g = 0, C_D = 0: Δξ = %s + i·%s° vs i·atan(ΔW/V) = i·%s°', ...
          f4(real(dxi) / DEG, 4), f4(imag(dxi) / DEG, 4), f4(imag(expd) / DEG, 4)), ...
        sprintf('ramiona po skoku: |K_F| = %s°, |K_S| = %s°; ν_S/ν_F = %s, |s_S|/|s_F| = %s', f4(abs(a.KF) / DEG, 4), f4(abs(a.KS) / DEG, 4), f4(L.nuSA / L.nuFA, 4), f4(th, 4))});
  end

  function r = testT9()
    R0 = gustRun('step', 0); L = muzzle(R0, 800);
    sF = L.sFA; nuF = L.nuF0;
    if isempty(stepArms), stepArms = armsAfter(R0, 0.5); end
    cases = {'trapez', 1 / nuF, 'liniowe, τ_r = 1/ν_F'; 'cos', 1.5 / nuF, '1 − cos, τ_r = 1,5/ν_F'; 'trapez', 0.5 / nuF, 'liniowe, τ_r = 0,5/ν_F'};
    det = {sprintf('ν_F = %s Hz; odniesienie: podmuch skokowy w tej samej chwili; teoria |∫₀^τr f′(τ)·e^(−s_F·τ) dτ|', f4(nuF, 3))};
    ok = true; worst = 0;
    for k = 1:size(cases, 1)
      shape = cases{k, 1}; tr = cases{k, 2};
      R = gustRun(shape, tr);
      a = armsAfter(R, 0.5 + tr);
      KFstep = abs(stepArms.KF * exp(stepArms.f.sF * tr));
      meas = abs(a.KF) / KFstep; th = gustExcitation(shape, tr, sF, 400); dd = abs(meas - th);
      ok = ok && dd < 0.02; worst = max(worst, dd);
      det{end + 1} = sprintf('%s = %s ms: |K_F|/|K_F,skok| = %s (6DOF) vs %s (teoria) – różnica %s', cases{k, 3}, f4(tr * 1000, 2), f4(meas, 4), f4(th, 4), f4(dd, 4)); %#ok<AGROW>
    end
    r = st('pass', ok, 'value', sprintf('maks. różnica %s (próg 0,02)', f4(worst, 4)), 'details', det);
  end

  function r = testT10()
    maxD = 0; maxE = 0;
    wrap = @(x) atan2(sin(x), cos(x));
    for ps = [-2.5, -1, 0, 0.7, 2.9]
      for th = [-1.4, -0.5, 0, 0.3, 1.2]
        for ph = [-3, -1, 0, 1.5, 3]
          T1 = quatToT(quatFromEuler(ps, th, ph)); T2 = TfromEuler(ps, th, ph);
          maxD = max(maxD, max(abs(T1 - T2)));
          [e1, e2, e3] = eulerFromT(T2);
          maxE = max([maxE, abs(wrap(e1 - ps)), abs(e2 - th), abs(wrap(e3 - ph))]);
        end
      end
    end
    P = makeParams(A, base); wind = Wiatr(base.wind);
    y = initialState(P); y(8) = 0.3; y(9) = -0.2;
    dy = rhs(P, wind, 0, y);
    h = 1e-6; lp = y(10:13) + h * dy(10:13); lm = y(10:13) - h * dy(10:13);
    [pp, tp, fp] = eulerFromT(quatToT(lp)); [pm, tm, fm] = eulerFromT(quatToT(lm)); [~, t0] = eulerFromT(quatToT(y(10:13)));
    dpsi = (pp - pm) / (2 * h); dth = (tp - tm) / (2 * h); dph = (fp - fm) / (2 * h);
    kinE = max([abs(dpsi - y(9) / cos(t0)), abs(dth - y(8)), abs(dph)]);
    rc = Rf.meta.maxRollCorr;
    r = st('pass', maxD < 1e-12 && maxE < 1e-9 && kinE < 1e-6 && rc < 1e-10, ...
      'value', sprintf('maks. |ΔT| = %s, kinematyka %s 1/s, korekta φ ≤ %s rad', ex(maxD, 1), ex(kinE, 1), ex(rc, 1)), ...
      'details', {sprintf('T(λ̄(ψ,θ,φ)) vs (A.8) na 125 zestawach kątów: maks. różnica %s; powrót do kątów: %s rad', ex(maxD, 2), ex(maxE, 2)), ...
        sprintf('λ̇ z (8) przy Ω̄ = (−r·tan θ, q, r): ψ̇ − r/cos θ, θ̇ − q, φ̇ → maks. %s 1/s (φ̇ = 0: układ zero-roll)', ex(kinE, 2)), ...
        sprintf('maks. korekta kąta przechylenia układu całkowania po kroku RK4 (25 s lotu): %s rad', ex(rc, 2))});
  end

  function r = testT11()
    ref = [0, 288.15, 101325, 1.225, 340.294; 5000, 255.676, 54048, 0.736429, 320.545;
           11019.1, 216.65, 22632.1, 0.363918, 295.070; 20000, 216.65, 5529.31, 0.0889097, 295.070];
    worst = 0; det = {};
    for k = 1:size(ref, 1)
      [Tk, pk, rk, ak] = atmosphere(ref(k, 1));
      e = max([rel(Tk, ref(k, 2)), rel(pk, ref(k, 3)), rel(rk, ref(k, 4)), rel(ak, ref(k, 5))]); worst = max(worst, e);
      det{end + 1} = sprintf('h = %g m: T = %s K, p = %s Pa, ρ = %s kg/m³, a = %s m/s (maks. odchyłka %s)', ref(k, 1), f4(Tk, 3), f4(pk, 1), f4(rk, 6), f4(ak, 3), pct(e)); %#ok<AGROW>
    end
    r = st('pass', worst < 5e-4, 'value', sprintf('maks. odchyłka %s', pct(worst)), 'details', det);
  end

  function r = testT12()
    % współczynniki zamrożone (C_D = 0, C_spin = 0, ρ = const), grawitacja włączona, nurkowanie 60°:
    % V rośnie wyłącznie od grawitacji. Postać ramion po drodze (z H* i V̇ z 6DOF) dawałaby tu
    % fałszywy przyrost ramienia F w tempie ≈ −(g·sin θ/V)·Re[s_F/(2s_F + η − iω_g)].
    R = sim(st('test', st('frozen', true, 'gravity', true), 'fire', st('u0', 500, 'QE', -60, 'q0', 0.1, 'r0', 0), 'num', st('tMax', 6, 'samplesPerPeriod', 25)));
    R.post = postprocess(R);
    ch = R.ch; po = R.post; u = uniformXi(R, 0, 6); N = numel(u.z);
    idx = find(~ch.ev);
    Phi = zeros(N, 1);
    for k = 2:N, Phi(k) = Phi(k - 1) + 0.5 * (po.sFimA(idx(k)) + po.sFimA(idx(k - 1))) * u.dt; end
    L = bitor(round(8 * 2 * pi / po.sFimA(idx(1)) / u.dt), 1); h = (L - 1) / 2;
    w = 0.5 - 0.5 * cos(2 * pi * ((0:L - 1).' + 0.5) / L); ws = sum(w);
    worst = 0; det = {};
    for tc = [0.5, 1.5, 2.5, 3.5, 4.5, 5.5]
      c = 0; while c < N - 1 && u.ts(c + 1) < tc, c = c + 1; end      % indeks od 0
      if c < h || c >= N - h, continue; end
      jj = (c - h:c + h).' + 1;
      kf = sum(u.z(jj) .* exp(-1i * Phi(jj)) .* w);
      meas = abs(kf) / ws; lin = po.KFabs(idx(c + 1)); e = lin / meas - 1;
      worst = max(worst, abs(e));
      det{end + 1} = sprintf('t = %s s, V = %s m/s: |K_F| z 6DOF %s°, ramię liniowe %s° (%s)', f4(u.ts(c + 1), 1), f4(ch.V(idx(c + 1)), 1), f4(meas / DEG, 5), f4(lin / DEG, 5), pct(abs(e))); %#ok<AGROW>
    end
    extra = 0; cnt = 0;
    for k = 1:20:N      % średni fałszywy przyrost w postaci „po drodze”
      i = idx(k); s = complex(po.sFreA(i), po.sFimA(i));
      rr = real(s / complex(2 * real(s) + po.eta(i), 2 * imag(s) - po.wg(i)));
      extra = extra - R.P.g * ch.U2(i) / ch.V(i)^2 * rr; cnt = cnt + 1;
    end
    det = [{sprintf(['g ≠ 0, C_D = 0, C_spin = 0, ρ = const; u₀ = 500 m/s, QE = −60°, q₀ = 0,1 rad/s; V: %s → %s m/s; ' ...
      'amplituda ramienia F z 6DOF przez demodulację (okno Hanna, 8 okresów ν_F)'], f4(ch.V(1), 0), f4(ch.V(end), 0))}, det, ...
      {sprintf('postać po drodze dałaby dodatkowe tempo %s 1/s, czyli ok. %s %% po 5,5 s', f4(extra / cnt, 4), f4(100 * (exp(extra / cnt * 5.5) - 1), 1))}];
    r = st('pass', worst < 0.03 && u.uniform, 'value', sprintf('maks. różnica %s (próg 3 %%)', pct(worst)), 'details', det);
  end

  function r = testT13()
    % (a) C(α) = c₀ + c₁·α: samo ramię → sieczna C(a); małe ramię wokół trymu e → c₀ + 1,5·c₁·e
    c0 = -0.4; c1 = 8; Cf = @(a) c0 + c1 * a;
    r1 = describingMagnus(Cf, 0, 3 * DEG, 0); r2 = describingMagnus(Cf, 10 * DEG, 0, 0); [r3F, r3S] = describingMagnus(Cf, 0, 0, 6 * DEG);
    eA = max([abs(r1 - Cf(3 * DEG)), abs(r2 - (c0 + 1.5 * c1 * 10 * DEG)), abs(r3F - (c0 + 1.5 * c1 * 6 * DEG)), abs(r3S - Cf(6 * DEG))]);
    % (b) QE = 65°, skok wiatru bocznego 10 m/s w wierzchołku: ramię F z 6DOF (dopasowanie dwóch modów) vs ramię liniowe
    R = sim(st('fire', st('QE', 65), 'wind', st('gust', st('enabled', true, 'mode', 'time', 'branch', 'apex', 'A', 10, 'azimuth_deg', 90, 'W2', 0, 'shape', 'step', 'td', 0))));
    R.post = postprocess(R);
    ch = R.ch; po = R.post;
    ev = R.meta.events(strcmp({R.meta.events.type}, 'gust'));
    if isempty(ev), error('podmuch nie wystąpił'); end
    ta = ev(1).t;
    w0 = ta + 1; w1 = ta + 8;
    sMeas = logSlope(ch.t, po.KFm, w0, w1); sLin = logSlope(ch.t, po.KFabs, w0, w1);
    [pMv, pMt] = peak(ch.t, po.KFm, ta + 1, ta + 35); [pLv, pLt] = peak(ch.t, po.KFabs, ta + 1, ta + 35);
    md2 = R.P.m * R.P.d^2; aDF = 0; aSec = 0; a0_ = 0; m = 0;
    for i = find(ch.t >= w0 & ch.t <= w1).'
      zs = po.base(i) * (po.CLa(i) + md2 / R.P.Ix * po.CmagmSec(i));
      rr = orderRoots(quadRoots(complex(po.eta(i), -po.wg(i)), complex(-po.wa2(i), -po.wg(i) * zs)));
      aSec = aSec + real(rr(1)); aDF = aDF + po.sFreA(i); a0_ = a0_ + po.sFre0(i); m = m + 1;
    end
    ok = eA < 1e-5 && abs(sMeas - sLin) < 0.005 && abs(pLv / pMv - 1) < 0.15 && abs(pLt - pMt) < 3;
    iT = find(ch.t <= ta, 1, 'last');
    r = st('pass', ok, 'value', sprintf('przypadki graniczne: %s; tempo |K_F|: 6DOF %s, ramię liniowe %s 1/s; maksimum %s', ex(eA, 1), f4(sMeas, 4), f4(sLin, 4), pct(abs(pLv / pMv - 1))), ...
      'details', {sprintf('C(α) = −0,4 + 8·α: samo ramię 3° → C_F = %s (sieczna %s); małe ramię wokół trymu 10° → %s (c₀ + 1,5·c₁·e = %s); małe F wokół ramienia S 6° → %s, C_S = %s (sieczna)', ...
          f4(r1, 6), f4(Cf(3 * DEG), 6), f4(r2, 6), f4(c0 + 1.5 * c1 * 10 * DEG, 6), f4(r3F, 6), f4(r3S, 6)), ...
        sprintf('QE = 65°, podmuch w t = %s s (ℳ = %s): d ln|K_F|/dt w %s–%s s – 6DOF %s, ramię liniowe %s 1/s (próg 0,005)', f4(ta, 2), f4(ch.M(iT), 3), f4(w0, 1), f4(w1, 1), f4(sMeas, 4), f4(sLin, 4)), ...
        sprintf('σ_F średnio w tym oknie: funkcja opisująca %s, sieczna przy α_ref %s, α_tot → 0 %s 1/s', f4(aDF / m, 4), f4(aSec / m, 4), f4(a0_ / m, 4)), ...
        sprintf('maksimum |K_F|: 6DOF %s° przy %s s, ramię liniowe %s° przy %s s (progi 15 %%, 3 s)', f4(pMv / DEG, 4), f4(pMt, 1), f4(pLv / DEG, 4), f4(pLt, 1))});
  end

  function r = testT14()
    det = {};
    % (a) włączone, zerowe odchyłki i wiatr → wynik identyczny jak bez warstw
    on = deepMerge(base, st('num', st('tMax', 5))); on.met = d0.met; on.met.enabled = true;
    Ra = simulate(A, deepMerge(base, st('num', st('tMax', 5)))); Rb = simulate(A, on);
    dA = max(abs(Ra.ch.X1(end) - Rb.ch.X1(end)), abs(Ra.ch.alphaTot(end) - Rb.ch.alphaTot(end)));
    det{end + 1} = sprintf('(a) zerowe odchyłki: różnica X₁ i α_tot po 5 s = %s', ex(dA, 1));
    % (b) jednolite ΔT = +15 K, Δρ = −10 %: ρ = 0,9·ρ_wz, a² = γR(T_wz + 15)
    dv = deepMerge(base, st('num', st('tMax', 5))); dv.met = d0.met; dv.met.enabled = true;
    for k = 1:numel(dv.met.zones), dv.met.zones(k).dT = 15; dv.met.zones(k).dRho = -10; end
    Rc = simulate(A, dv); eB = 0;
    for i = 1:97:Rc.n
      [Ts, ~, rs] = atmosphere(Rc.P.hM + Rc.ch.X2(i));
      eB = max([eB, abs(Rc.ch.rho(i) / rs - 0.9), abs(Rc.ch.a(i)^2 / (GAMMA_AIR * R_AIR * (Ts + 15)) - 1)]);
    end
    det{end + 1} = sprintf('(b) ΔT = +15 K, Δρ = −10 %%: maks. błąd ρ i a² w próbkach = %s', ex(eB, 1));
    % (c) wiatr: azymut 30°; z 30° (czołowy) → W₁ = −v; z 300° (z lewej) → W₃ = +v
    cw = base; cw.earth.az_deg = 30; cw.met = d0.met; cw.met.enabled = true; cw.met.dirUnit = 'deg';
    cw.met.zones(1) = struct('dir', 30, 'spd', 10, 'dT', 0, 'dRho', 0); cw.met.zones(2) = struct('dir', 300, 'spd', 10, 'dT', 0, 'dRho', 0);
    Pw = makeParams(A, cw);
    eC = max([abs(Pw.met.w1(1) + 10), abs(Pw.met.w3(1)), abs(Pw.met.w1(2)), abs(Pw.met.w3(2) - 10)]);
    det{end + 1} = sprintf('(c) wiatr 10 m/s: z kierunku strzału W₁ = %s, z lewej W₃ = %s m/s', f4(Pw.met.w1(1), 3), f4(Pw.met.w3(2), 3));
    % (d) komunikat syntetyczny: T_v = T_wz + 2 K, p = 0,99·p_wz w środkach stref 00–15, MDP 120 m
    hMDP = 120; Z = metcmZones(); msgT = sprintf('METCM1 520210 281200 012013\n');
    for k = 0:15
      [Ts, ps] = atmosphere(hMDP + Z.mid(k + 1));
      msgT = [msgT, sprintf('%02d%03d%03d %04d%04d\n', k, 100 + 10 * k, 5 + k, round((Ts + 2) * 10), round(0.99 * ps / 100))]; %#ok<AGROW>
    end
    pm = parseMETCM(msgT, 0); eD = 0; eW = 0;
    for z = pm.zones
      [Ts, ps, rs] = atmosphere(hMDP + Z.mid(z.n + 1)); rhoEx = 0.99 * ps / (R_AIR * (Ts + 2));
      tolR = 100 * (0.5 / z.p_hPa) + 100 * 0.05 / z.Tv + 1e-9;             % zaokrąglenie PPPP do 1 hPa, TTTT do 0,1 K
      eD = max([eD, abs(z.dT - 2) / (0.05 + 1e-9), abs(z.dRho - 100 * (rhoEx / rs - 1)) / tolR]);
      eW = max([eW, abs(z.dir_mil - 10 * (100 + 10 * z.n)), abs(z.spd_ms - 0.514444 * (5 + z.n))]);
    end
    det{end + 1} = sprintf('(d) komunikat syntetyczny: %d stref, MDP = %g m (%s); ΔT, Δρ w granicach zaokrągleń zapisu: %s (wykorzystanie %s), wiatr: %s', ...
      numel(pm.zones), pm.hMDP, pm.hMDPsource, tn(eD <= 1), pct(eD), ex(eW, 1));
    % (e) PCHIP: dane monotoniczne bez przerzutów, ciągła pochodna w węźle
    ip = pchipMake([0 1 3 4 7], [0 1 1.2 4 4.5]); over = 0;
    for x = 0:0.01:7, v = pchipVal(ip, x); over = max([over, -v, v - 4.5]); end
    jump = abs(pchipDer(ip, 3 - 1e-7) - pchipDer(ip, 3 + 1e-7));
    det{end + 1} = sprintf('(e) PCHIP: przerzut ponad zakres danych %s, skok pochodnej w węźle %s', ex(max(0, over), 1), ex(jump, 1));
    ok = dA == 0 && eB < 1e-12 && eC < 1e-12 && numel(pm.zones) == 16 && pm.hMDP == 120 && eD <= 1 && eW < 1e-9 && over <= 1e-12 && jump < 1e-5;
    r = st('pass', ok, 'value', sprintf('zerowe odchyłki: %s; atmosfera: %s; komunikat: %d stref', ex(dA, 1), ex(eB, 1), numel(pm.zones)), 'details', det);
  end
end

% =====================================================================
function s = st(varargin)
% skrót do budowy struktur konfiguracji (wartości tablicowe opakowane, by nie tworzyć tablic struktur)
s = struct();
for k = 1:2:numel(varargin), s.(varargin{k}) = varargin{k + 1}; end
end

function m = meanOver(R, a, b, re, im)
sel = R.ch.t >= a & R.ch.t <= b & ~R.ch.ev;
m = complex(mean(re(sel)), mean(im(sel)));
end

function k = linSlope(x, y)
n = numel(x); k = (n * sum(x .* y) - sum(x) * sum(y)) / (n * sum(x .^ 2) - sum(x)^2);
end

function k = logSlope(t, y, a, b)
sel = t >= a & t <= b & y > 0;
k = linSlope(t(sel), log(y(sel)));
end

function [v, tt] = peak(t, y, a, b)
sel = find(t >= a & t <= b & y > 0);
[v, j] = max(y(sel)); tt = t(sel(j));
end

function s = f4(x, d)
if nargin < 2, d = 4; end
if ~isfinite(x), s = '–'; return; end
s = strrep(sprintf('%.*f', d, x), '.', ',');
end

function s = pct(x)
if x < 0.001, d = 4; else, d = 3; end
s = [strrep(sprintf('%.*f', d, 100 * x), '.', ','), ' %'];
end

function s = ex(x, d)
if x == 0, s = '0'; else, s = sprintf('%.*e', d, x); end
end

function s = tn(b)
if b, s = 'tak'; else, s = 'NIE'; end
end
