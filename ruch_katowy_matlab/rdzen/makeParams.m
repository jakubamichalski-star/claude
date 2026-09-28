function P = makeParams(A, cfgIn)
%MAKEPARAMS Konfiguracja (scalona z domyślną) i parametry pochodne obliczeń.
%   P = makeParams(A, cfg) – A z buildAero, cfg – struktura (pola jak w defaultConfig; brakujące – domyślne).
if nargin < 2, cfgIn = struct(); end
DEG = pi / 180; OMEGA_E = 7.292115e-5;
cfg = deepMerge(defaultConfig(), cfgIn);
P.cfg = cfg; P.A = A;
if strcmp(cfg.norm, 'pd/V'), P.kp = 1; else, P.kp = 0.5; end
P.m = cfg.mass.m; P.Ix = cfg.mass.Ix * cfg.mult.Ix; P.IT = cfg.mass.IT * cfg.mult.IT;
P.d = A.d; P.S = A.S;
P.dXcg = (cfg.mass.XCG_mm - A.refCG_mm) / A.refDia_mm;   % ΔX_CG [kal.], dodatnie ku tyłowi
P.f = cfg.f; P.fd = cfg.fd;
P.fullCmqTransfer = logical(cfg.fullCmqTransfer);
P.frozen = logical(cfg.test.frozen);
if P.frozen, P.f.CD0 = 0; P.f.CD2 = 0; P.f.Cspin = 0; P.fd.CD0 = 0; P.fd.Cspin = 0; end
fdv = struct2cell(P.fd); P.anyFd = any([fdv{:}] ~= 0);
% tablice do szybkiej interpolacji współczynników (aeroCoefVec)
P.mach = A.mach(:);
P.coefTab = [A.CDo, A.CD2, A.CLa, A.CLa3, A.Cma, A.Cma3, A.Cma5, A.CNa, A.CNq, A.Cmq, A.Cmq2, A.CYpa, A.Clp, A.CPN];
P.cnpaAlpha = A.cnpaAlpha(:).'; P.cnpaData = A.cnpaData;
P.fv = [P.f.CD0, P.f.CD2, P.f.CLa, P.f.CMa, P.f.CMq, P.f.CNq, P.f.Cmagf, P.f.Cmagm, P.f.Cspin];
gravity = isfield(cfg.test, 'gravity') && logical(cfg.test.gravity);
if P.frozen && ~gravity, P.g = 0; else, P.g = cfg.earth.g; end   % test.gravity: współczynniki zamrożone, grawitacja włączona
P.hM = cfg.fire.hM;
P.met = metProfile(cfg.met, cfg.earth.az_deg, P.hM);
P.coriolis = ~P.frozen && logical(cfg.earth.coriolis);
if P.coriolis
  L = cfg.earth.lat_deg * DEG; az = cfg.earth.az_deg * DEG;
  P.LamF = [OMEGA_E * cos(L) * cos(az), OMEGA_E * sin(L), -OMEGA_E * cos(L) * sin(az)];
  P.LamI = [P.LamF(1), P.LamF(3), -P.LamF(2)];
else
  P.LamF = [0 0 0]; P.LamI = [0 0 0];
end
F = cfg.fire;
switch F.QEunit
  case 'mil6400', P.QE = F.QE * 2 * pi / 6400;
  case 'mil6000', P.QE = F.QE * 2 * pi / 6000;
  otherwise, P.QE = F.QE * DEG;
end
P.u0 = F.u0;
P.p0 = 2 * pi * P.u0 / (F.tc * F.bore_mm / 1000) * cfg.mult.p0;   % p₀ = 2π·u₀/(t_c·kaliber lufy)
end
