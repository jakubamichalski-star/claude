classdef Wiatr < handle
%WIATR Wiatr W̄ = W̄₀(h) + ΔW̄·f (składowe w układzie ogniowym).
%   w = Wiatr(Wcfg, met, hOff)
%   h – wysokość nad poziomem wylotu (X₂); h_M tylko dla atmosfery.
%   W̄₀ – stały, z tabeli [h W1 W2 W3] albo z warstw METCM (met; hOff = h_M − h_MDP).
%   f – kształt podmuchu (skok / trapez / 1 − cos), w funkcji czasu od przecięcia h₁
%   albo wysokości (warstwa ze strefami przejściowymi).
%   Stan „piece”: 0 przed, 1 narastanie, 2 pełny podmuch, 3 zanikanie, 4 po.
%   Zmiany stanu wyłącznie w zdarzeniach – integrator trafia w nie dokładnie.

  properties
    on = false; mode = 'time'; dW = [0; 0; 0]; h1 = 0; branch = 'asc'; shape = 'step';
    tr = 0; td = Inf; dh = 0; HL = Inf;
    st = struct('piece', 0, 'tg', NaN, 'trigIdx', 0, 'afterApex', false, 'pulses', 0);
    met = []; hO = 0; W0c = [0; 0; 0]; tab = [];
    atApex = false; trigDirs = 1; hasTrigTime = false; triggerAt = NaN;
    bnds = struct('h', {}, 'lo', {}, 'hi', {});
  end

  methods
    function w = Wiatr(Wcfg, met, hOff)
      DEG = pi / 180;
      if nargin < 2, met = []; end
      if nargin < 3, hOff = 0; end
      w.met = met; w.hO = hOff;
      base = struct('W', [0 0 0], 'table', []);
      if isstruct(Wcfg) && isfield(Wcfg, 'base') && isstruct(Wcfg.base), base = Wcfg.base; end
      if isfield(base, 'W') && ~isempty(base.W), w.W0c = base.W(:); end
      if isfield(base, 'table') && size(base.table, 1) >= 2, w.tab = sortrows(base.table, 1); end
      G = struct();
      if isstruct(Wcfg) && isfield(Wcfg, 'gust') && isstruct(Wcfg.gust), G = Wcfg.gust; end
      w.shape = gf(G, 'shape', 'step');
      chi = gf(G, 'azimuth_deg', 0) * DEG;                  % kierunek „dokąd wieje”, od osi 1 ku osi 3
      A = gf(G, 'A', 0); W2 = gf(G, 'W2', 0);
      w.dW = [A * cos(chi); W2; A * sin(chi)];
      if strcmp(w.shape, 'step'), w.tr = 0; else, w.tr = max(0, gf(G, 'tr', 0)); end
      td = gf(G, 'td', 0); if td > 0, w.td = td; else, w.td = Inf; end   % ≤ 0: bez końca
      HL = gf(G, 'HL', 0); if HL > 0, w.HL = HL; else, w.HL = Inf; end
      if strcmp(w.shape, 'step'), w.dh = 0; else, w.dh = max(0, gf(G, 'dh', 0)); end
      if isfinite(w.HL), w.dh = min(w.dh, w.HL / 2); end
      w.h1 = gf(G, 'h1', 0);
      w.on = logical(gf(G, 'enabled', false)) && (abs(A) > 0 || abs(W2) > 0);
      if strcmp(gf(G, 'mode', 'time'), 'height'), w.mode = 'height'; else, w.mode = 'time'; end
      w.branch = gf(G, 'branch', 'asc');
      if strcmp(w.branch, 'apex') && strcmp(w.mode, 'height'), w.branch = 'both'; end   % wierzchołek tylko przy zadaniu czasem
      w.atApex = strcmp(w.branch, 'apex');
      if w.atApex, w.trigDirs = [];
      elseif strcmp(w.branch, 'asc'), w.trigDirs = 1;
      elseif strcmp(w.branch, 'desc'), w.trigDirs = -1;
      else, w.trigDirs = [1 -1];
      end
      w.triggerAt = gf(G, 'triggerAt', NaN);                 % tylko testy: wyzwolenie w zadanej chwili
      w.hasTrigTime = isfinite(w.triggerAt);
      % granice warstwy (tryb wysokościowy)
      if strcmp(w.mode, 'height')
        h1 = w.h1; dh = w.dh; HL_ = w.HL;
        if dh > 0, w.addB(h1, 0, 1); else, w.addB(h1, 0, 2); end
        if dh > 0, w.addB(h1 + dh, 1, 2); end
        if isfinite(HL_)
          if dh > 0, w.addB(h1 + HL_ - dh, 2, 3); w.addB(h1 + HL_, 3, 4);
          else, w.addB(h1 + HL_, 2, 4);
          end
        end
      end
    end

    function init(w, h0)
      if ~w.on, return; end
      if strcmp(w.mode, 'height'), w.st.piece = w.pieceFromH(h0);
      elseif ~w.hasTrigTime && ~isempty(w.trigDirs) && w.trigDirs(1) == 1 && w.h1 <= h0, w.startPulse(0);
      end
    end

    function out = base(w, h)
      % wiatr bazowy W̄₀(h) [3×1]
      if ~isempty(w.met), out = metW(w.met, h + w.hO); return; end   % warstwy METCM
      tb = w.tab;
      if isempty(tb), out = w.W0c; return; end
      n = size(tb, 1);
      if h <= tb(1, 1), out = tb(1, 2:4).'; return; end
      if h >= tb(n, 1), out = tb(n, 2:4).'; return; end
      i = find(tb(:, 1) < h, 1, 'last');
      a = (h - tb(i, 1)) / (tb(i + 1, 1) - tb(i, 1));
      out = (tb(i, 2:4) + a * (tb(i + 1, 2:4) - tb(i, 2:4))).';
    end

    function out = baseGrad(w, h)
      % dW̄₀/dh [3×1]
      if ~isempty(w.met), out = metdW(w.met, h + w.hO); return; end
      out = [0; 0; 0]; tb = w.tab;
      if isempty(tb) || h <= tb(1, 1) || h >= tb(end, 1), return; end
      i = find(tb(:, 1) < h, 1, 'last');
      out = ((tb(i + 1, 2:4) - tb(i, 2:4)) / (tb(i + 1, 1) - tb(i, 1))).';
    end

    function [out, f] = W(w, t, h)
      % całkowity wiatr w układzie ogniowym [3×1] i kształt podmuchu f
      if ~w.on && isempty(w.met) && isempty(w.tab), out = w.W0c; f = 0; return; end   % wiatr stały
      out = w.base(h);
      f = w.fval(t, h);
      if f ~= 0, out = out + f * w.dW; end
    end

    function v = dfdt(w, t, h, U2)
      % df/dt wzdłuż toru (U₂ = dX₂/dt)
      v = 0;
      if ~w.on, return; end
      pc = w.st.piece;
      if pc ~= 1 && pc ~= 3, return; end
      if pc == 1, sg = 1; else, sg = -1; end
      if strcmp(w.mode, 'time')
        if pc == 1, x = (t - w.st.tg) / w.tr; else, x = (t - w.st.tg - w.tr - w.td) / w.tr; end
        v = sg * shapeRiseDer(w.shape, x) / w.tr;
        return;
      end
      if pc == 1, x = (h - w.h1) / w.dh; else, x = (h - (w.h1 + w.HL - w.dh)) / w.dh; end
      v = sg * shapeRiseDer(w.shape, x) / w.dh * U2 * w.gate();
    end

    function b = inTransition(w)
      b = w.on && (w.st.piece == 1 || w.st.piece == 3) && (strcmp(w.mode, 'time') || w.gate() == 1);
    end

    function tb = nextTimeBreak(w, t)
      tb = Inf;
      if ~w.on || ~strcmp(w.mode, 'time'), return; end
      switch w.st.piece
        case 0
          if w.hasTrigTime && w.st.trigIdx == 0 && w.triggerAt > t - 1e-12, tb = w.triggerAt; end
        case 1, tb = w.st.tg + w.tr;
        case 2, tb = w.st.tg + w.tr + w.td;
        case 3, tb = w.st.tg + 2 * w.tr + w.td;
      end
    end

    function changed = onTimeBreak(w, t)
      % przejście w chwili znanej z góry; zwraca true, gdy zmienił się wiatr
      changed = true;
      pc = w.st.piece;
      if pc == 0, w.startPulse(t); return; end
      if pc == 1, w.st.piece = 2; return; end
      if pc == 2
        if w.tr > 0, w.st.piece = 3; else, w.st.piece = 4; end
      elseif pc == 3
        w.st.piece = 4;
      end
      if w.st.piece == 4 && w.st.trigIdx < numel(w.trigDirs) && ~w.hasTrigTime, w.st.piece = 0; end  % „obie” gałęzie: czekaj na drugie przecięcie
    end

    function r = onApex(w, t)
      % r = 'start' (początek podmuchu w wierzchołku), true (zmiana bramki gałęzi) albo false
      g0 = w.gate(); w.st.afterApex = true;
      if w.on && strcmp(w.mode, 'time') && w.atApex && ~w.hasTrigTime && w.st.piece == 0 && w.st.trigIdx == 0
        w.startPulse(t); r = 'start'; return;
      end
      r = w.on && strcmp(w.mode, 'height') && w.gate() ~= g0 && w.st.piece ~= 0;
    end

    function ev = events(w)
      % zdarzenia stanu (przecięcia wysokości) aktywne w bieżącym stanie: name, level, dir, bi
      ev = struct('name', {}, 'level', {}, 'dir', {}, 'bi', {});
      if ~w.on, return; end
      if strcmp(w.mode, 'time')
        if w.st.piece == 0 && ~w.hasTrigTime && w.st.trigIdx < numel(w.trigDirs)
          ev(1) = struct('name', 'gust-trigger', 'level', w.h1, 'dir', w.trigDirs(w.st.trigIdx + 1), 'bi', 0);
        end
      else
        for k = 1:numel(w.bnds)
          ev(end + 1) = struct('name', 'gust-layer', 'level', w.bnds(k).h, 'dir', 0, 'bi', k); %#ok<AGROW>
        end
      end
    end

    function text = fire(w, ev, t, up)
      % obsługa zdarzenia z events(); zwraca opis
      if strcmp(ev.name, 'gust-trigger')
        w.startPulse(t);
        if ev.dir > 0, gal = 'wznosząca'; else, gal = 'opadająca'; end
        text = sprintf('początek podmuchu (h = h₁ = %g m, gałąź %s)', w.h1, gal);
      else
        b = w.bnds(ev.bi);
        if up, w.st.piece = b.hi; else, w.st.piece = b.lo; end
        if w.gate(), s = ''; else, s = ' (na tej gałęzi toru podmuch nieaktywny)'; end
        text = sprintf('granica warstwy podmuchu h = %.1f m%s', b.h, s);
      end
    end
  end

  methods (Access = private)
    function addB(w, h, lo, hi)
      w.bnds(end + 1) = struct('h', h, 'lo', lo, 'hi', hi);
    end
    function pc = pieceFromH(w, h)
      pc = 0;
      for k = 1:numel(w.bnds), if h >= w.bnds(k).h, pc = w.bnds(k).hi; end, end
    end
    function g = gate(w)
      if ~strcmp(w.mode, 'height') || strcmp(w.branch, 'both'), g = 1; return; end
      if strcmp(w.branch, 'asc'), g = double(~w.st.afterApex); else, g = double(w.st.afterApex); end
    end
    function f = fval(w, t, h)
      f = 0;
      if ~w.on, return; end
      pc = w.st.piece;
      if strcmp(w.mode, 'time')
        if pc == 1, f = shapeRise(w.shape, (t - w.st.tg) / w.tr);
        elseif pc == 2, f = 1;
        elseif pc == 3, f = 1 - shapeRise(w.shape, (t - w.st.tg - w.tr - w.td) / w.tr);
        end
        return;
      end
      if pc == 1, f = shapeRise(w.shape, (h - w.h1) / w.dh);
      elseif pc == 2, f = 1;
      elseif pc == 3, f = 1 - shapeRise(w.shape, (h - (w.h1 + w.HL - w.dh)) / w.dh);
      end
      f = f * w.gate();
    end
    function startPulse(w, t)
      w.st.tg = t;
      if w.tr > 0, w.st.piece = 1; else, w.st.piece = 2; end
      w.st.trigIdx = w.st.trigIdx + 1; w.st.pulses = w.st.pulses + 1;
    end
  end
end

function v = gf(G, name, def)
if isfield(G, name) && ~isempty(G.(name)), v = G.(name); else, v = def; end
end

function v = shapeRise(shape, x)
if x <= 0, v = 0; return; end
if x >= 1, v = 1; return; end
if strcmp(shape, 'cos'), v = 0.5 * (1 - cos(pi * x));
elseif strcmp(shape, 'step'), v = 1;
else, v = x;
end
end

function v = shapeRiseDer(shape, x)
if x <= 0 || x >= 1, v = 0; return; end
if strcmp(shape, 'cos'), v = 0.5 * pi * sin(pi * x);
elseif strcmp(shape, 'step'), v = 0;
else, v = 1;
end
end
