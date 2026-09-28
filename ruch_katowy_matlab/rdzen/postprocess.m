function o = postprocess(R)
%POSTPROCESS Przetwarzanie przebiegu 6DOF: wielkości sekcji 4 w każdej próbce,
%   ξ = α_tot·e^{iφ_A}, ξ_e (pełne rozwiązanie szczególne), ᾱ_e (STANAG),
%   ramiona K_F, K_S (inicjalizacja u wylotu i po każdej zmianie podmuchu),
%   całki ∫s dt + czynnik WKB, obwiednie |K_F| ± |K_S|, kod stanu.
%   o = postprocess(R) – R z simulate; pola jak w postChannels (wektory kolumnowe).
DEG = pi / 180;
P = R.P; ch = R.ch; n = R.n; num = P.cfg.num; ST = statusBits();
names = postChannels();
for k = 1:numel(names), o.(names{k}) = NaN(n, 1); end
md2 = P.m * P.d * P.d;
aWarn = num.alphaWarn_deg * DEG; aeWarn = num.alphaEWarn_deg * DEG;
wkb = ~isfield(num, 'wkb') || ~isequal(num.wkb, false);
if isfield(num, 'sigmaEps') && ~isempty(num.sigmaEps), sEps = num.sigmaEps; else, sEps = 0.02; end  % martwa strefa |σ| [1/s]
zetaE = NaN(n, 1);
t = ch.t;
isTrans = ch.piece == 1 | ch.piece == 3;

% ---- przebieg 0: α_ref = √⟨α_tot²⟩ z najdłuższego okresu fazy względnej (1/(ν_F − ν_S) albo 1/ν_S).
%      Służy już tylko jako wariant zapasowy (sieczna C_mag-m(α_ref)), gdy ramion nie da się rozdzielić. ----
aRef = NaN(n, 1);
I = zeros(n, 1);
for i = 2:n
  a0 = ch.alphaTot(i - 1); a1 = ch.alphaTot(i);
  I(i) = I(i - 1) + 0.5 * (a0 * a0 + a1 * a1) * (t(i) - t(i - 1));
end
j = 1;
for i = 1:n
  c = aeroCoefVec(P, ch.M(i), 0);
  wg = P.Ix * ch.p(i) / P.IT; wa2 = ch.rho(i) * P.S * P.d * ch.V(i)^2 * c(3) / (2 * P.IT);
  if wg ~= 0, inv = 4 * wa2 / (wg * wg); else, inv = Inf; end
  if inv < 1
    r = sqrt(1 - inv); dw = min(abs(wg) * r, abs(wg) * (1 - r) / 2);
  else
    dw = max([abs(wg) / 2, sqrt(abs(wa2)), 1e-3]);
  end
  Tw = min(max(2 * pi / max(dw, 1e-3), 0.02), 3);
  ta = t(i) - Tw;
  while j < i && t(j + 1) <= ta, j = j + 1; end
  if t(i) - t(1) < 1e-12, aRef(i) = ch.alphaTot(i); continue; end
  if t(j) >= ta || j == i
    Ia = I(j); tA = t(j);
  else
    w = (ta - t(j)) / (t(j + 1) - t(j)); Ia = I(j) + w * (I(j + 1) - I(j)); tA = ta;
  end
  span = t(i) - tA;
  if span > 0, aRef(i) = sqrt(max(0, (I(i) - Ia) / span)); else, aRef(i) = ch.alphaTot(i); end
end
o.alphaRef = aRef;

% ---- przebieg 1: wielkości sekcji 4 dla α_tot → 0, ᾱ_e (STANAG), ξ_e quasi-stacjonarne ----
prevF = []; prevS = [];
for i = 1:n
  M = ch.M(i); V = ch.V(i); rho = ch.rho(i); p = ch.p(i);
  L = linPoint(P, M, V, rho, p, aRef(i));
  sF0 = L.sF0; sS0 = L.sS0;
  if ~isempty(prevF) && abs(sF0 - prevF) + abs(sS0 - prevS) > abs(sF0 - prevS) + abs(sS0 - prevF)
    q_ = sF0; sF0 = sS0; sS0 = q_;
  end
  prevF = sF0; prevS = sS0;
  o.invSg(i) = L.invSg; o.Sg(i) = L.Sg; o.Sd0(i) = L.Sd0;
  o.nuF0(i) = L.nuF0; o.nuS0(i) = L.nuS0; o.sigF0(i) = L.sigF0; o.sigS0(i) = L.sigS0;
  o.sFre0(i) = real(sF0); o.sFim0(i) = imag(sF0); o.sSre0(i) = real(sS0); o.sSim0(i) = imag(sS0);
  o.wg(i) = L.wg; o.wa2(i) = L.wa2; o.eta(i) = L.eta; o.zeta0(i) = L.zeta0; o.base(i) = L.base;
  o.CLa(i) = L.CLa; o.CD0(i) = L.CD0; o.CMa(i) = L.CMa; o.CMq(i) = L.CMq; o.CNq(i) = L.CNq; o.Cmagf(i) = L.Cmagf;
  o.Cmagm0(i) = L.Cmagm0; o.CmagmSec(i) = L.CmagmA; o.Cspin(i) = L.Cspin; o.CPN(i) = L.CPN;
  B = stabBands(P, M, ch.alphaTot(i), L); o.SgSig(i) = B.SgSig; o.SdSig(i) = B.SdSig;
  o.xiRe(i) = ch.alphaTot(i) * cos(ch.phiA(i)); o.xiIm(i) = ch.alphaTot(i) * sin(ch.phiA(i));
  % ᾱ_e wg STANAG: wektor w układzie ogniowym → Re = nos w górę (−z̄), Im = nos w prawo (ȳ)
  ax = zrAxesF(ch.psi(i), ch.theta(i));
  UF = [ch.U1(i), ch.U2(i), ch.U3(i)]; VF = UF - [ch.W1(i), ch.W2(i), ch.W3(i)];
  ae = alphaEStanag(P, M, rho, p, VF, UF);
  o.aeRe(i) = -dot(ae, ax.z); o.aeIm(i) = dot(ae, ax.y); o.aeAbs(i) = norm(ae);
  % ξ_e0 = −(i·ω_g − η_q)·Γ/(ω_α² + i·ω_g·ζ), Γ = (−γ_z + i·γ_y)/V, γ̄ = ḡ + Ā_C − dW̄₀/dt;
  % rzut γ̄ na układ zero-roll wektora V̄ (bez nutacyjnego „drgania” osi), ζ przy α = |ξ_e| (jedna iteracja)
  if ~isTrans(i)
    Ac = [0 0 0]; if P.coriolis, Ac = -2 * cross(P.LamF, UF); end
    gam = [Ac(1) - ch.dW0dt1(i), -P.g + Ac(2) - ch.dW0dt2(i), Ac(3) - ch.dW0dt3(i)];
    axV = zrAxesF(atan2(VF(3), VF(1)), asin(max(-1, min(1, VF(2) / V))));
    Gam = complex(-dot(axV.z, gam) / V, dot(axV.y, gam) / V);
    num_ = complex(L.etaQ, -L.wg) * Gam;
    zeta = L.zeta0; xe = num_ / complex(L.wa2, L.wg * zeta);
    a1 = abs(xe);
    if a1 > 0 && isfinite(a1)
      cB = aeroCoefVec(P, M, a1);
      zeta = L.base * (L.CLa + (md2 / P.Ix) * cB(7));
      xe = num_ / complex(L.wa2, L.wg * zeta);
    end
    o.xie0Re(i) = real(xe); o.xie0Im(i) = imag(xe); zetaE(i) = zeta;
  end
end

% ---- przebieg 2: ξ̇_e0 (w odcinkach bez zdarzeń) i poprawka pierwszego rzędu
%      na powolną zmienność: ξ_e = ξ_e0 + (η_t − i·ω_g)·ξ̇_e0/(ω_α² + i·ω_g·ζ) ----
xedRe = zeros(n, 1); xedIm = zeros(n, 1);
segs = zeros(0, 2);                  % odcinki [a, b): t rosnące, bez stref przejściowych podmuchu
i = 1;
while i <= n
  if ~isfinite(o.xie0Re(i)), i = i + 1; continue; end
  j = i + 1;
  while j <= n && isfinite(o.xie0Re(j)) && t(j) - t(j - 1) > 1e-9, j = j + 1; end
  segs(end + 1, :) = [i, j]; %#ok<AGROW>
  i = j;
end
% ξ̇_e0: nachylenie lokalnej prostej MNK w oknie ±0,25 s (obciętym do odcinka). C_mag-m(α) jest
% odcinkami liniowa, więc ξ_e0 ma załamania w węzłach tabeli; różnica trójpunktowa dawałaby w tych
% miejscach skok poprawki ξ_e (a za nim fałszywe ramię w dopasowaniu modów).
Hd = 0.25;
fr = o.xie0Re; fi = o.xie0Im;
for s = 1:size(segs, 1)
  a = segs(s, 1); b = segs(s, 2);
  if b - a < 2, continue; end
  j0 = a; j1 = a;
  for i = a:b - 1
    while j0 < i && t(j0) < t(i) - Hd, j0 = j0 + 1; end
    while j1 + 1 < b && t(j1 + 1) <= t(i) + Hd, j1 = j1 + 1; end
    lo = j0; hi = j1;
    if hi - lo < 2, lo = max(a, min(lo, i - 1)); hi = min(b - 1, max(hi, i + 1)); end
    tt = t(lo:hi); tm = mean(tt); dtt = tt - tm;
    stt = sum(dtt .* dtt);
    if stt > 0
      xedRe(i) = sum(dtt .* (fr(lo:hi) - mean(fr(lo:hi)))) / stt;
      xedIm(i) = sum(dtt .* (fi(lo:hi) - mean(fi(lo:hi)))) / stt;
    end
  end
end
% η_t – współczynnik przy ξ̇ w równaniu w dziedzinie czasu: (ρSV/2m)[C_Lα − (md²/I_T)(C_Mq + C_Mα̇)],
% bez −C_D0 (opór działa wzdłuż V̄ i nie obraca wektora prędkości; −C_D0 w H* pojawia się dopiero
% po przejściu do drogi jako zmiennej niezależnej, przy założeniu, że V̇ pochodzi tylko od oporu)
etaT = o.eta + o.base .* o.CD0;
for i = 1:n
  if ~isfinite(o.xie0Re(i)), continue; end
  corr = complex(etaT(i), -o.wg(i)) * complex(xedRe(i), xedIm(i)) / complex(o.wa2(i), o.wg(i) * zetaE(i));
  o.xieRe(i) = o.xie0Re(i) + real(corr); o.xieIm(i) = o.xie0Im(i) + imag(corr);
end

% ---- przebieg 3: amplitudy ramion z 6DOF i funkcja opisująca C_mag-m.
%      Lokalne dopasowanie MNK dwóch modów: ξ − ξ_e ≈ K_F·e^{iΦ_F} + K_S·e^{iΦ_S} (Φ_j = ∫ω_j dt z pierwiastków
%      przy α → 0), waga Hanna, okno ±1,5 okresu dudnień 1/(ν_F − ν_S), obcięte do odcinka bez zdarzeń.
%      Dla |K_F|, |K_S| i |ξ_e| – współczynniki efektywne C_F, C_S (osobno dla każdego ramienia).
%      Wariant zapasowy (sieczna przy α_ref), gdy ramion nie da się rozdzielić (S_g ≤ 1, dudnienia wolniejsze
%      niż 1 Hz, za krótki odcinek, źle uwarunkowany układ) oraz w strefach przejściowych podmuchu. ----
PhF = zeros(n, 1); PhS = zeros(n, 1);
for i = 2:n
  dt = t(i) - t(i - 1);
  PhF(i) = PhF(i - 1) + 0.5 * (o.sFim0(i) + o.sFim0(i - 1)) * dt;
  PhS(i) = PhS(i - 1) + 0.5 * (o.sSim0(i) + o.sSim0(i - 1)) * dt;
end
eFa = exp(1i * PhF); eSa = exp(1i * PhS);
CFa = o.CmagmSec; CSa = o.CmagmSec; o.dfOK(:) = 0;
zRes = complex(o.xiRe - o.xieRe, o.xiIm - o.xieIm);
for s = 1:size(segs, 1)
  a = segs(s, 1); b = segs(s, 2);
  cen = [];
  c = a;
  while c < b
    cen(end + 1) = c; %#ok<AGROW>
    dw = abs(o.sFim0(c) - o.sSim0(c));
    if dw > 0, Hw = min(1.5 * 2 * pi / dw, 1.0); else, Hw = 0.1; end
    c2 = c + 1; while c2 < b && t(c2) - t(c) < Hw / 4, c2 = c2 + 1; end
    c = c2;
  end
  if cen(end) ~= b - 1, cen(end + 1) = b - 1; end %#ok<AGROW>
  nc = numel(cen); val = NaN(nc, 4);       % [a_F, a_S, C_F, C_S]; NaN – brak dopasowania
  for k = 1:nc
    f = armFit(cen(k), a, b); if isempty(f), continue; end
    [CF, CS] = describingMagnus(cmagmAtMach(P, ch.M(cen(k))), abs(complex(o.xieRe(cen(k)), o.xieIm(cen(k))) + f(3)), f(1), f(2));
    val(k, :) = [f(1), f(2), CF, CS];
  end
  for k = 1:nc
    c0 = cen(k);
    if k < nc, c1 = cen(k + 1); v1 = val(k + 1, :); else, c1 = c0 + 1; v1 = val(k, :); end
    v0 = val(k, :);
    if isnan(v0(1)) || isnan(v1(1)), continue; end
    for i = c0:min(c1, b) - 1
      if c1 > c0 && c1 <= n && t(c1) > t(c0), w = (t(i) - t(c0)) / (t(c1) - t(c0)); else, w = 0; end
      o.KFm(i) = v0(1) + w * (v1(1) - v0(1)); o.KSm(i) = v0(2) + w * (v1(2) - v0(2));
      CFa(i) = v0(3) + w * (v1(3) - v0(3)); CSa(i) = v0(4) + w * (v1(4) - v0(4)); o.dfOK(i) = 1;
    end
  end
end

% ---- przebieg 4: wariant „przy bieżącej amplitudzie”: ζ_F, ζ_S, pierwiastki, S_d,F i S_d,S, kod stanu ----
pF = []; pS = [];
for i = 1:n
  base = o.base(i); eta = o.eta(i); wg = o.wg(i); wa2 = o.wa2(i);
  zF = base * (o.CLa(i) + (md2 / P.Ix) * CFa(i)); zS = base * (o.CLa(i) + (md2 / P.Ix) * CSa(i));
  o.CmagmF(i) = CFa(i); o.CmagmS(i) = CSa(i); o.zetaF(i) = zF; o.zetaS(i) = zS;
  sF = pick(orderRoots(quadRoots(complex(eta, -wg), complex(-wa2, -wg * zF))), pF, true);
  sS = pick(orderRoots(quadRoots(complex(eta, -wg), complex(-wa2, -wg * zS))), pS, false);
  pF = sF; pS = sS;
  o.sFreA(i) = real(sF); o.sFimA(i) = imag(sF); o.sSreA(i) = real(sS); o.sSimA(i) = imag(sS);
  o.SdFA(i) = 2 * zF / eta; o.SdSA(i) = 2 * zS / eta;
  if o.invSg(i) < 1, rt = sqrt(1 - o.invSg(i)); else, rt = NaN; end
  o.nuFA(i) = o.nuF0(i); o.nuSA(i) = o.nuS0(i);
  o.sigFA(i) = -(eta / 2) * (1 - (o.SdFA(i) - 1) / rt); o.sigSA(i) = -(eta / 2) * (1 + (o.SdSA(i) - 1) / rt);
  st = 0;
  if o.invSg(i) > 1, st = bitor(st, ST.GYRO); end
  if real(sF) > sEps, st = bitor(st, ST.DYN_F); end
  if real(sS) > sEps, st = bitor(st, ST.DYN_S); end
  if o.aeAbs(i) > aeWarn, st = bitor(st, ST.WEAK); end
  if ch.alphaTot(i) > aWarn, st = bitor(st, ST.AWARN); end
  if o.sFre0(i) > 0, st = bitor(st, ST.DYN_F0); end
  if o.sSre0(i) > 0, st = bitor(st, ST.DYN_S0); end
  if ch.clamp(i), st = bitor(st, ST.MACH); end
  if (bitand(st, ST.DYN_F0) && ~bitand(st, ST.DYN_F)) || (bitand(st, ST.DYN_S0) && ~bitand(st, ST.DYN_S)), st = bitor(st, ST.LIMIT); end
  o.status(i) = st;
end

% ---- przebieg 5: ramiona K_F, K_S – inicjalizacja u wylotu i po każdej zmianie podmuchu.
%      Ewolucja w dziedzinie czasu: K_j(t) = K_j(t₀)·exp(∫s_j dt)·A_j, s_F – pierwiastek F równania
%      s² + (η_t − iω_g)s − (ω_α² + iω_g·ζ_F) = 0, s_S – pierwiastek S równania z ζ_S;
%      czynnik WKB d ln A_j = −ds_j/(2s_j + η_t − iω_g). (Postać po drodze z H* i V̇ z 6DOF dawałaby
%      fałszywe „pompowanie” ramienia F o ≈ −g·sin θ/V przy zmianie prędkości od grawitacji.) ----
K = []; pF = []; pS = [];
for i = 1:n
  xe = complex(o.xieRe(i), o.xieIm(i));
  sF = pick(orderRoots(quadRoots(complex(etaT(i), -o.wg(i)), complex(-o.wa2(i), -o.wg(i) * o.zetaF(i)))), pF, true);
  sS = pick(orderRoots(quadRoots(complex(etaT(i), -o.wg(i)), complex(-o.wa2(i), -o.wg(i) * o.zetaS(i)))), pS, false);
  pF = sF; pS = sS;
  if isTrans(i) || ~isfinite(real(xe)), K = []; continue; end
  bt = complex(etaT(i), -o.wg(i));
  xi = complex(o.xiRe(i), o.xiIm(i));
  if i == 1 || ch.reinit(i) == 1 || isempty(K)
    x0 = xi - xe; xd0 = complex(ch.dxiRe(i) - xedRe(i), ch.dxiIm(i) - xedIm(i));
    KF0 = (xd0 - sS * x0) / (sF - sS);
    K = struct('KF0', KF0, 'KS0', x0 - KF0, 'IF', 0, 'IS', 0, 'lnAF', 0, 'lnAS', 0, 't', t(i), 'sF', sF, 'sS', sS, 'bt', bt);
  else
    dt = t(i) - K.t;
    if dt > 0
      K.IF = K.IF + (sF + K.sF) * 0.5 * dt;
      K.IS = K.IS + (sS + K.sS) * 0.5 * dt;
      if wkb   % wartości środkowe: 2s̄ + η̄_t − iω̄_g
        bm = (bt + K.bt) / 2;
        K.lnAF = K.lnAF - (sF - K.sF) / (sF + K.sF + bm);
        K.lnAS = K.lnAS - (sS - K.sS) / (sS + K.sS + bm);
      end
      K.t = t(i); K.sF = sF; K.sS = sS; K.bt = bt;
    end
  end
  KF = K.KF0 * exp(K.IF + K.lnAF);
  KS = K.KS0 * exp(K.IS + K.lnAS);
  o.KFre(i) = real(KF); o.KFim(i) = imag(KF); o.KSre(i) = real(KS); o.KSim(i) = imag(KS);
  aF = abs(KF); aS = abs(KS);
  o.KFabs(i) = aF; o.KSabs(i) = aS; o.envP(i) = aF + aS; o.envM(i) = abs(aF - aS);
  o.xiLinRe(i) = real(xe + KF + KS); o.xiLinIm(i) = imag(xe + KF + KS);
  o.devE(i) = abs(xi - xe);
end

  % -------------------------------------------------------------
  function f = armFit(c, a, b)
    % Lokalne dopasowanie: ξ − ξ_e ≈ [c] + K_F·e^{iΦ_F} + K_S·e^{iΦ_S}. Wyraz stały c (przesunięcie środka ruchu
    % względem quasi-stacjonarnego ξ_e) tylko wtedy, gdy okno obejmuje ≥ 1,2 okresu składowej wolnej (T_S ≤ 3 s);
    % w krótszym oknie ramię wolne i stałe przesunięcie są nierozróżnialne.
    % Wynik [|K_F|, |K_S|, c] albo [], gdy ramion nie da się rozdzielić.
    f = [];
    dwc = abs(o.sFim0(c) - o.sSim0(c));
    if ~(o.invSg(c) < 1) || ~(dwc > 2 * pi), return; end            % dudnienia wolniejsze niż 1 Hz
    Tb = 2 * pi / dwc;
    if abs(o.sSim0(c)) > 0, TS = 2 * pi / abs(o.sSim0(c)); else, TS = Inf; end
    use3 = TS <= 3;
    if use3, Hwc = min(max(1.5 * Tb, 0.6 * TS), 1.5); else, Hwc = min(1.5 * Tb, 1.0); end
    tc = t(c);
    jj0 = c; while jj0 > a && t(jj0 - 1) >= tc - Hwc, jj0 = jj0 - 1; end
    jj1 = c; while jj1 < b - 1 && t(jj1 + 1) <= tc + Hwc, jj1 = jj1 + 1; end
    if jj1 - jj0 < 8 || t(jj1) - t(jj0) < 0.5 * Tb, return; end
    three = use3 && t(jj1) - t(jj0) >= 1.0 * TS;
    idx = jj0:jj1;
    wv = 0.5 + 0.5 * cos(pi * (t(idx) - tc) / Hwc);
    z = zRes(idx); eF = eFa(idx); eS = eSa(idx);
    Wt = sum(wv); F0 = sum(wv .* eF); S0 = sum(wv .* eS); G12 = sum(wv .* conj(eF) .* eS);
    z0 = sum(wv .* z); r1 = sum(wv .* conj(eF) .* z); r2 = sum(wv .* conj(eS) .* z);
    if three
      x = solveC([Wt, F0, S0; conj(F0), Wt, G12; conj(S0), conj(G12), Wt], [z0; r1; r2]);
      if ~isempty(x), f = [abs(x(2)), abs(x(3)), x(1)]; return; end
    end
    x = solveC([Wt, G12; conj(G12), Wt], [r1; r2]);
    if isempty(x), return; end
    f = [abs(x(1)), abs(x(2)), 0];
  end
end

% ---------------------------------------------------------------------
function x = solveC(G, r)
% układ normalny MNK (zespolony, hermitowski) – eliminacja Gaussa z wyborem elementu głównego;
% [] przy źle uwarunkowanym układzie
m = numel(r); Aug = [G, r(:)];
dmax = max(abs(diag(G)));
for k = 1:m
  [~, pk] = max(abs(Aug(k:m, k))); pk = pk + k - 1;
  if ~(abs(Aug(pk, k)) > 1e-3 * dmax), x = []; return; end
  if pk ~= k, Aug([pk k], :) = Aug([k pk], :); end
  for q = k + 1:m
    fq = Aug(q, k) / Aug(k, k);
    Aug(q, k:m + 1) = Aug(q, k:m + 1) - fq * Aug(k, k:m + 1);
  end
end
x = zeros(m, 1);
for k = m:-1:1
  s = Aug(k, m + 1) - Aug(k, k + 1:m) * x(k + 1:m);
  x(k) = s / Aug(k, k);
end
end

function s = pick(rt, prevR, wantF)
% wybór pierwiastka przez ciągłość (F – pierwszy z orderRoots, S – drugi)
if isempty(prevR)
  if wantF, s = rt(1); else, s = rt(2); end
  return;
end
if abs(rt(1) - prevR) <= abs(rt(2) - prevR), s = rt(1); else, s = rt(2); end
end
